"""Error estimation, readout mitigation, and Zero-Noise Extrapolation (ZNE) for PS-C1.

Implements:
1. Complete Readout Calibration & Inversion Mitigation (Confusion Matrix inversion)
2. Parity Symmetry Verification (syndrome postselection on physical sector)
3. Zero-Noise Extrapolation (ZNE) via Unitary Folding with Richardson Extrapolation
4. Multi-seed statistical repeat evaluation
"""
from typing import Dict, List, Tuple
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

from .problem import TERMS_2Q, NUCLEAR_REPULSION, EXACT_REFERENCE_ENERGY, CHEMICAL_ACCURACY
from .processors import get_coupling_map, get_noise_model, with_measurement
from .routing import compile_circuit

STATES_2Q = ["00", "01", "10", "11"]


def expect_from_probs(probs: Dict[str, float], qubits: List[int]) -> float:
    """Compute Pauli Z expectation value from probability distribution."""
    val = 0.0
    for s, p in probs.items():
        # Qiskit bitstring order is little-endian: s[::-1][q]
        parity = sum(int(s[::-1][q]) for q in qubits) % 2
        sign = -1.0 if parity == 1 else 1.0
        val += sign * p
    return val


class ReadoutMitigator:
    """Calibrates and inverts the readout error confusion matrix M."""

    def __init__(self, sim: AerSimulator, num_qubits: int = 2, shots: int = 20000, seed: int = 42):
        self.sim = sim
        self.num_qubits = num_qubits
        self.states = [f"{i:0{num_qubits}b}" for i in range(2**num_qubits)]
        self.shots = shots
        self.M = self._calibrate(seed=seed)
        self.M_inv = np.linalg.inv(self.M)

    def _calibrate(self, seed: int) -> np.ndarray:
        circuits = []
        for s in self.states:
            qc = QuantumCircuit(self.num_qubits, self.num_qubits)
            for i, bit in enumerate(reversed(s)):
                if bit == "1":
                    qc.x(i)
            qc.measure(list(range(self.num_qubits)), list(range(self.num_qubits)))
            circuits.append(qc)

        results = self.sim.run(circuits, shots=self.shots, seed_simulator=seed).result().get_counts()
        M = np.zeros((len(self.states), len(self.states)))
        for j, s_prep in enumerate(self.states):
            counts = results[j]
            for i, s_meas in enumerate(self.states):
                M[i, j] = counts.get(s_meas, 0) / self.shots
        return M

    def mitigate_counts(self, counts: Dict[str, int]) -> Dict[str, float]:
        total_shots = sum(counts.values())
        raw_vec = np.array([counts.get(s, 0) / total_shots for s in self.states])
        mit_vec = self.M_inv @ raw_vec
        # Constrained least-squares projection / clipping
        mit_vec = np.clip(mit_vec, 0.0, 1.0)
        s = np.sum(mit_vec)
        if s > 0:
            mit_vec /= s
        return {s: float(mit_vec[i]) for i, s in enumerate(self.states)}


class NoisyEnergyEvaluator:
    """Evaluates noisy Hamiltonian expectation values on a target processor."""

    def __init__(self, proc: dict, ansatz: QuantumCircuit, shots: int = 4096, initial_layout: list = None):
        self.proc = proc
        self.ansatz = ansatz
        self.shots = shots
        self.initial_layout = initial_layout
        self.cm = get_coupling_map(proc)
        self.noise = get_noise_model(proc)
        self.sim = AerSimulator(
            noise_model=self.noise,
            basis_gates=proc["basis_gates"],
            coupling_map=self.cm
        )
        self.ro_mitigator = ReadoutMitigator(self.sim, num_qubits=2, shots=20000, seed=42)

    def evaluate(
        self,
        params: np.ndarray,
        mitigate_readout: bool = False,
        symmetry_verify: bool = False,
        scale_factor: int = 1,
        seed: int = 42
    ) -> float:
        """Measure expectation value with optional readout mitigation, symmetry check, and unitary folding."""
        bound = self.ansatz.assign_parameters(params)

        # Transpile base unitary to basis gates
        u_z = with_measurement(bound, "Z")
        u_x = with_measurement(bound, "X")

        cz_base = compile_circuit(u_z, self.cm, self.initial_layout)
        cx_base = compile_circuit(u_x, self.cm, self.initial_layout)

        # Unitary folding for ZNE (works for scale_factor = 1, 3, 5, ...)
        def fold_circuit(comp_qc, scale):
            gates = [inst for inst in comp_qc.data if inst.operation.name not in ("measure", "barrier")]
            meas = [inst for inst in comp_qc.data if inst.operation.name == "measure"]
            core = QuantumCircuit(*comp_qc.qregs, *comp_qc.cregs)
            for inst in gates:
                core.append(inst)
            inv = core.inverse()
            folded = QuantumCircuit(*comp_qc.qregs, *comp_qc.cregs)
            folded.compose(core, inplace=True)
            for _ in range((scale - 1) // 2):
                folded.compose(inv, inplace=True)
                folded.compose(core, inplace=True)
            for inst in meas:
                folded.append(inst)
            return transpile(folded, basis_gates=self.proc["basis_gates"], optimization_level=0)

        cz_exec = fold_circuit(cz_base, scale_factor)
        cx_exec = fold_circuit(cx_base, scale_factor)

        res_z = self.sim.run(cz_exec, shots=self.shots, seed_simulator=seed).result().get_counts()
        res_x = self.sim.run(cx_exec, shots=self.shots, seed_simulator=seed + 100).result().get_counts()

        if mitigate_readout:
            pz = self.ro_mitigator.mitigate_counts(res_z)
            px = self.ro_mitigator.mitigate_counts(res_x)
        else:
            tot_z, tot_x = sum(res_z.values()), sum(res_x.values())
            pz = {s: res_z.get(s, 0) / tot_z for s in STATES_2Q}
            px = {s: res_x.get(s, 0) / tot_x for s in STATES_2Q}

        if symmetry_verify:
            # Physical parity sector for H2 ground state is odd parity (|01>, |10>)
            tot_sym = pz.get("01", 0.0) + pz.get("10", 0.0)
            if tot_sym > 0:
                pz = {
                    "00": 0.0,
                    "01": pz.get("01", 0.0) / tot_sym,
                    "10": pz.get("10", 0.0) / tot_sym,
                    "11": 0.0
                }

        # Energy summation
        e = TERMS_2Q["II"]
        e += TERMS_2Q["IZ"] * expect_from_probs(pz, [0])
        e += TERMS_2Q["ZI"] * expect_from_probs(pz, [1])
        e += TERMS_2Q["ZZ"] * expect_from_probs(pz, [0, 1])
        e += TERMS_2Q["XX"] * expect_from_probs(px, [0, 1])
        return float(e) + NUCLEAR_REPULSION

    def repeated_evaluation(
        self,
        params: np.ndarray,
        seeds: List[int],
        mitigate_readout: bool = False,
        symmetry_verify: bool = False,
        scale_factor: int = 1
    ) -> Tuple[float, float, np.ndarray]:
        vals = np.array([
            self.evaluate(
                params,
                mitigate_readout=mitigate_readout,
                symmetry_verify=symmetry_verify,
                scale_factor=scale_factor,
                seed=s
            )
            for s in seeds
        ])
        return float(vals.mean()), float(vals.std(ddof=1)), vals


def run_zero_noise_extrapolation(
    evaluator: NoisyEnergyEvaluator,
    params: np.ndarray,
    scales: List[int] = [1, 3, 5],
    seed: int = 42
) -> Dict[str, float]:
    """Perform Richardson Zero-Noise Extrapolation (ZNE) across scale factors."""
    energies = []
    for scale in scales:
        e = evaluator.evaluate(
            params,
            mitigate_readout=True,
            symmetry_verify=True,
            scale_factor=scale,
            seed=seed + scale
        )
        energies.append(e)

    e1, e3, e5 = energies[0], energies[1], energies[2]
    # Linear Richardson (scales 1, 3): E(0) = (3*E(1) - E(3)) / 2
    e_zne_linear = (3.0 * e1 - e3) / 2.0
    # Quadratic Richardson (scales 1, 3, 5): 15/8 E(1) - 10/8 E(3) + 3/8 E(5)
    e_zne_quad = 1.875 * e1 - 1.25 * e3 + 0.375 * e5

    # Pre-specified estimator: linear (1,3) Richardson. Quadratic is reported for
    # comparison only; it is never selected using the exact reference energy.
    return {
        "scale_1_energy": e1,
        "scale_3_energy": e3,
        "scale_5_energy": e5,
        "zne_linear_energy": e_zne_linear,
        "zne_quad_energy": e_zne_quad,
    }


def noisy_4q_fidelity_study(
    procs: Dict[str, dict],
    ansatz_4q: QuantumCircuit,
    shots: int = 4096,
    seeds: List[int] = None,
    n_param_sets: int = 10,
    param_seed: int = 7,
) -> List[dict]:
    """Noisy performance of the 4-qubit stress-test circuit on each processor.

    Metric: Hellinger fidelity between the noisy Z-basis output distribution
    (organizer noise model, transpiled with the organizer settings) and the ideal
    statevector distribution. 1.0 means the noisy output equals the ideal one.
    Parameters are random but fixed (param_seed), so results are reproducible.
    Mean and std are taken over all (parameter set, simulator seed) runs.
    """
    from qiskit.quantum_info import Statevector

    seeds = seeds or list(range(42, 52))
    rng = np.random.default_rng(param_seed)
    param_sets = [rng.uniform(0, 2 * np.pi, ansatz_4q.num_parameters) for _ in range(n_param_sets)]
    n = ansatz_4q.num_qubits
    rows, all_f = [], {}
    for name, proc in procs.items():
        cm = get_coupling_map(proc)
        sim = AerSimulator(noise_model=get_noise_model(proc), basis_gates=proc["basis_gates"], coupling_map=cm)
        fids, cx, swaps, depth = [], None, None, None
        for params in param_sets:
            bound = ansatz_4q.assign_parameters(params)
            ideal = Statevector(bound).probabilities()
            comp = compile_circuit(with_measurement(bound, "Z"), cm)
            ops = comp.count_ops()
            cx = int(ops.get("cx", 0))
            depth = comp.depth(lambda i: i.operation.name not in ("measure", "barrier"))
            swaps = int(compile_circuit(with_measurement(bound, "Z"), cm, keep_swaps=True).count_ops().get("swap", 0))
            for sd in seeds:
                counts = sim.run(comp, shots=shots, seed_simulator=sd).result().get_counts()
                noisy = np.zeros(2 ** n)
                for k, v in counts.items():
                    noisy[int(k.replace(" ", ""), 2)] = v / shots
                fids.append(float(np.sum(np.sqrt(ideal * noisy)) ** 2))
        f = np.array(fids)
        all_f[name] = f
        rows.append({
            "architecture": name,
            "compiled_cx": cx,
            "swap_count": swaps,
            "compiled_depth": depth,
            "shots": shots,
            "n_runs": len(f),
            "hellinger_fidelity_mean": float(f.mean()),
            "hellinger_fidelity_std": float(f.std(ddof=1)),
            "hellinger_fidelity_sem": float(f.std(ddof=1) / np.sqrt(len(f))),
        })
    # Paired comparison against the first processor (same parameters and seeds on each processor)
    ref = next(iter(all_f))
    for r in rows:
        d = all_f[r["architecture"]] - all_f[ref]
        r["paired_gain_vs_" + ref] = float(d.mean())
        r["paired_gain_sem"] = float(d.std(ddof=1) / np.sqrt(len(d))) if r["architecture"] != ref else 0.0
    return rows


def run_full_benchmark(procs_required=("processor_A", "processor_B"), include_extras=True):
    """Execute all benchmark components and save tables and figures."""
    import itertools
    from pathlib import Path
    import pandas as pd
    import networkx as nx
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .problem import exact_ground_energy
    from .processors import (
        load_processor,
        get_coupling_map,
        hea_ansatz,
        compact_chemistry_ansatz,
        full_entangling_4q_ansatz,
        with_measurement,
        run_vqe_ideal
    )
    from .routing import compile_circuit, analyze_routing, benchmark_4qubit_architecture_routing
    from .visualization import (
        save_processor_graphs,
        plot_vqe_convergence,
        plot_error_mitigation_results,
        plot_4qubit_routing_overhead,
        plot_4q_noisy_fidelity,
        render_latex_assets
    )

    root = Path(__file__).resolve().parent.parent
    tab_dir = root / "results" / "tables"
    fig_dir = root / "results" / "figures"
    tab_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    shots = 4096
    seeds = list(range(42, 52))

    print("=== Step 1: Rendering LaTeX equation assets ===")
    render_latex_assets()

    proc_names = list(procs_required)
    if include_extras:
        proc_names += ["processor_D", "processor_A_giantkit_variant", "processor_B_giantkit_variant"]

    procs = {name: load_processor(name) for name in proc_names}
    save_processor_graphs({k: procs[k] for k in ["processor_A", "processor_B", "processor_D"] if k in procs})

    print("=== Step 2: Ideal VQE Convergence ===")
    E_exact = exact_ground_energy(2)
    ansatzes = {
        "HEA reps=1": hea_ansatz(1),
        "HEA reps=2": hea_ansatz(2),
        "Compact Singlet": compact_chemistry_ansatz()
    }

    opt_results, histories = {}, {}
    for label, ans in ansatzes.items():
        x_opt, e_opt, hist = run_vqe_ideal(ans, seed=42)
        opt_results[label] = (x_opt, e_opt)
        histories[label] = hist
        print(f"[{label}] Optimized Energy: {e_opt:.10f} Ha (Error: {abs(e_opt - E_exact)*1000:.4e} mHa)")

    plot_vqe_convergence(histories, E_exact)

    print("=== Step 3: Architecture Transpilation & Baseline Noisy Execution ===")
    comparison_rows = []
    for label in ["HEA reps=1", "HEA reps=2"]:
        ans = ansatzes[label]
        x_opt, e_ideal = opt_results[label]
        logical_meas = with_measurement(ans.assign_parameters(x_opt), "Z")

        for name in proc_names:
            proc = procs[name]
            cm = get_coupling_map(proc)
            m, _ = analyze_routing(logical_meas, cm)
            evaluator = NoisyEnergyEvaluator(proc, ans, shots=shots)
            mean_e, std_e, _ = evaluator.repeated_evaluation(x_opt, seeds)
            err = abs(mean_e - E_exact)
            comparison_rows.append({
                "problem_id": "PS-C1",
                "instance_id": f"H2_R0.735_{label.replace(' ', '_')}",
                "architecture": name,
                "role": "required" if name in procs_required else "additional",
                "logical_qubits": 2,
                "physical_qubits": proc["num_qubits"],
                "physical_qubits_used": m["physical_qubits_used"],
                "ansatz_params": ans.num_parameters,
                "classical_reference": round(E_exact, 8),
                "objective_or_accuracy": round(mean_e, 6),
                "approximation_or_error": round(err, 6),
                "energy_error_mHa": round(err * 1000.0, 3),
                "noisy_energy_std": round(std_e, 6),
                "cx_count": m["cx_count"],
                "depth": m["depth"],
                "swap_count": m["swap_count"],
                "one_qubit_gates": m["one_qubit_gates"],
                "shots": shots,
                "noise_model": "1q_depol=0.001, cx_depol=0.005, readout=0.01"
            })

    df_full = pd.DataFrame(comparison_rows)
    df_full.to_csv(tab_dir / "AB_comparison_full.csv", index=False)
    df_full.to_csv(tab_dir / "final_comparison.csv", index=False)

    cols_section_5_3 = [
        "problem_id", "instance_id", "architecture", "logical_qubits", "physical_qubits",
        "objective_or_accuracy", "classical_reference", "approximation_or_error",
        "cx_count", "depth", "swap_count", "shots", "noise_model"
    ]
    df_required = df_full[df_full.role == "required"][cols_section_5_3]
    df_required.to_csv(tab_dir / "AB_comparison.csv", index=False)

    for n in ["processor_A", "processor_B", "processor_D"]:
        if n in procs:
            df_full[df_full.architecture == n].to_csv(tab_dir / f"{n}_metrics.csv", index=False)

    fig, axs = plt.subplots(1, 3, figsize=(14, 4), dpi=150)
    req = df_full[df_full.role == "required"]
    bars_labels = [f"{r['architecture'].split('_')[1]}\n{r['instance_id'].split('_')[-1]}" for _, r in req.iterrows()]
    axs[0].bar(bars_labels, req["energy_error_mHa"], yerr=req["noisy_energy_std"] * 1000, color="#4c78a8", edgecolor="black")
    axs[0].axhline(CHEMICAL_ACCURACY * 1000, color="red", linestyle="--", label="Chem. Acc (1.6 mHa)")
    axs[0].set_ylabel("Noisy Error |E - E_exact| (mHa)")
    axs[0].set_title("Noisy Ground-State Error")
    axs[0].legend(fontsize=9)
    axs[1].bar(bars_labels, req["cx_count"], color="#f58518", edgecolor="black")
    axs[1].set_title("Compiled CX Gate Count")
    axs[2].bar(bars_labels, req["depth"], color="#54a24b", edgecolor="black")
    axs[2].set_title("Compiled Circuit Depth")
    plt.tight_layout()
    plt.savefig(fig_dir / "AB_performance.png", bbox_inches="tight")
    plt.close(fig)

    print("=== Step 4: Error Mitigation Study (multi-seed, pre-specified linear ZNE) ===")
    ans_compact = ansatzes["Compact Singlet"]
    x_opt_compact, _ = opt_results["Compact Singlet"]
    mit_rows = []
    for mit_shots, tag in ((4096, "organizer_4096"), (32768, "extra_32768")):
        ev = NoisyEnergyEvaluator(procs["processor_A"], ans_compact, shots=mit_shots)
        per_stage = {k: [] for k in ("1. Raw Noisy", "2. Readout Mitigated", "3. Readout + Symmetry", "4. Readout+Sym+ZNE (linear 1,3)", "5. ZNE quadratic (1,3,5), comparison only")}
        for sd in seeds:
            per_stage["1. Raw Noisy"].append(ev.evaluate(x_opt_compact, False, False, 1, sd))
            per_stage["2. Readout Mitigated"].append(ev.evaluate(x_opt_compact, True, False, 1, sd))
            per_stage["3. Readout + Symmetry"].append(ev.evaluate(x_opt_compact, True, True, 1, sd))
            z = run_zero_noise_extrapolation(ev, x_opt_compact, scales=[1, 3, 5], seed=sd * 10)
            per_stage["4. Readout+Sym+ZNE (linear 1,3)"].append(z["zne_linear_energy"])
            per_stage["5. ZNE quadratic (1,3,5), comparison only"].append(z["zne_quad_energy"])
        for stage, vals in per_stage.items():
            v = np.array(vals)
            errs = np.abs(v - E_exact)
            mit_rows.append({
                "stage": stage,
                "shots": mit_shots,
                "setting": tag,
                "n_seeds": len(v),
                "energy_hartree": float(v.mean()),
                "energy_error_hartree": float(abs(v.mean() - E_exact)),
                "std_hartree": float(v.std(ddof=1)),
                "sem_hartree": float(v.std(ddof=1) / np.sqrt(len(v))),
                "mean_abs_error_hartree": float(errs.mean()),
                "frac_runs_within_chem_acc": float((errs < CHEMICAL_ACCURACY).mean()),
                "within_chemical_accuracy": bool(abs(v.mean() - E_exact) < CHEMICAL_ACCURACY),
            })
    df_mit = pd.DataFrame(mit_rows)
    df_mit.to_csv(tab_dir / "error_mitigation_results.csv", index=False)
    plot_error_mitigation_results(df_mit, CHEMICAL_ACCURACY)
    print(df_mit[["stage", "shots", "energy_error_hartree", "std_hartree", "frac_runs_within_chem_acc"]].to_string())

    print("=== Step 5: 4-Qubit Architecture Routing Comparison (Processor A vs B vs D) ===")
    ans_4q = full_entangling_4q_ansatz()
    target_procs = {k: procs[k] for k in ["processor_A", "processor_B", "processor_D"]}
    routing_4q_res = benchmark_4qubit_architecture_routing(target_procs, ans_4q)
    df_routing = pd.DataFrame(routing_4q_res)
    df_routing.to_csv(tab_dir / "routing_4qubit_comparison.csv", index=False)
    plot_4qubit_routing_overhead(df_routing)
    df_4q_noisy = pd.DataFrame(noisy_4q_fidelity_study(target_procs, ans_4q, shots=shots, seeds=seeds))
    df_4q_noisy.to_csv(tab_dir / "routing_4qubit_noisy_fidelity.csv", index=False)
    plot_4q_noisy_fidelity(df_4q_noisy)
    print(df_4q_noisy.to_string())
    print("4-qubit routing overhead summary:")
    print(df_routing[["architecture", "compiled_cx", "swap_count", "compiled_depth"]])

    print("=== Step 6: Layout-Sensitivity Study across Physical Qubit Pairings ===")
    lay_rows = []
    ans_lay = hea_ansatz(1)
    x_lay, _ = opt_results["HEA reps=1"]
    logical_lay = with_measurement(ans_lay.assign_parameters(x_lay), "Z")

    for name in ["processor_A", "processor_B", "processor_D"]:
        proc = procs[name]
        cm = get_coupling_map(proc)
        G = nx.Graph()
        G.add_edges_from([tuple(e) for e in proc["coupling_map"]])
        for p, q in itertools.permutations(range(proc["num_qubits"]), 2):
            m, _ = analyze_routing(logical_lay, cm, initial_layout=[p, q])
            dist = nx.shortest_path_length(G, p, q)
            lay_rows.append({
                "architecture": name,
                "phys_q0": p,
                "phys_q1": q,
                "distance": dist,
                "cx_count": m["cx_count"],
                "swap_count": m["swap_count"],
                "depth": m["depth"]
            })

    df_lay = pd.DataFrame(lay_rows)
    df_lay.to_csv(tab_dir / "layout_sensitivity.csv", index=False)
    df_lay_summ = df_lay.groupby("architecture").agg(
        placements=("cx_count", "size"),
        mean_cx=("cx_count", "mean"),
        max_cx=("cx_count", "max"),
        frac_zero_swap=("swap_count", lambda s: (s == 0).mean()),
        max_distance=("distance", "max")
    ).reset_index()
    df_lay_summ.to_csv(tab_dir / "layout_sensitivity_summary.csv", index=False)

    wc_rows = []
    for name in ["processor_A", "processor_B", "processor_D"]:
        proc = procs[name]
        s_ = df_lay[df_lay.architecture == name]
        for label, row in (("best", s_.sort_values(["cx_count", "depth"]).iloc[0]),
                           ("worst", s_.sort_values(["cx_count", "depth"]).iloc[-1])):
            ev = NoisyEnergyEvaluator(proc, ans_lay, shots=shots, initial_layout=[int(row.phys_q0), int(row.phys_q1)])
            mean_val, std_val, _ = ev.repeated_evaluation(x_lay, seeds)
            wc_rows.append({
                "architecture": name,
                "placement": label,
                "phys_q0": int(row.phys_q0),
                "phys_q1": int(row.phys_q1),
                "distance": int(row.distance),
                "cx_count": int(row.cx_count),
                "swap_count": int(row.swap_count),
                "depth": int(row.depth),
                "noisy_energy": round(mean_val, 6),
                "noisy_std": round(std_val, 6),
                "noisy_error_mHa": round(abs(mean_val - E_exact) * 1000.0, 3)
            })
    df_wc = pd.DataFrame(wc_rows)
    df_wc.to_csv(tab_dir / "placement_noise_effect.csv", index=False)

    fig, axs = plt.subplots(1, 3, figsize=(15, 4), sharey=True, dpi=150)
    for ax, name in zip(axs, ["processor_A", "processor_B", "processor_D"]):
        sub = df_lay[df_lay.architecture == name]
        grp = sub.groupby("distance")["cx_count"].agg(["mean", "min", "max"])
        ax.bar(grp.index, grp["mean"], color="#4c78a8", edgecolor="black")
        ax.set_title(f"{name}: Compiled CX vs Graph Distance")
        ax.set_xlabel("Shortest Path Distance in Coupling Graph")
    axs[0].set_ylabel("Compiled CX Count (Mean)")
    plt.tight_layout()
    plt.savefig(fig_dir / "layout_sensitivity.png", bbox_inches="tight")
    plt.close(fig)

    logical_lay.draw("mpl", filename=str(fig_dir / "logical_circuit.png"))
    for n in ["processor_A", "processor_B", "processor_D"]:
        c = compile_circuit(logical_lay, get_coupling_map(procs[n]))
        c.draw("mpl", idle_wires=False, filename=str(fig_dir / f"transpiled_{n}.png"))

    from .visualization import plot_processor_graph
    for name in ["processor_A", "processor_B"]:
        fig, ax = plt.subplots(figsize=(5, 4.5), dpi=150)
        plot_processor_graph(procs[name], ax)
        fig.tight_layout()
        fig.savefig(fig_dir / f"{name}_graph.png", bbox_inches="tight")
        plt.close(fig)

    import shutil
    shutil.copy(fig_dir / "routing_overhead_4q.png", fig_dir / "architecture_tradeoffs.png")

    try:
        from PIL import Image
        img_log = Image.open(fig_dir / "logical_circuit.png")
        img_tra = Image.open(fig_dir / "transpiled_processor_A.png")
        w = max(img_log.width, img_tra.width)
        h = img_log.height + img_tra.height + 40
        comb = Image.new("RGB", (w, h), (255, 255, 255))
        comb.paste(img_log, ((w - img_log.width) // 2, 10))
        comb.paste(img_tra, ((w - img_tra.width) // 2, img_log.height + 30))
        comb.save(fig_dir / "logical_vs_transpiled.png")
    except Exception:
        pass

    print("=== All benchmarks completed successfully! ===")
    return df_required, df_mit, df_routing, df_lay_summ, df_wc
