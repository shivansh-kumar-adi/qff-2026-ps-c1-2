"""Regenerate README.md and report.md with every number read from results/tables/*.csv.

Usage:  python scripts/build_docs.py     (run after src.run_full_benchmark())
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
T = ROOT / "results" / "tables"
ab = pd.read_csv(T / "AB_comparison_full.csv")
mit = pd.read_csv(T / "error_mitigation_results.csv")
rt = pd.read_csv(T / "routing_4qubit_comparison.csv")
lay = pd.read_csv(T / "layout_sensitivity_summary.csv").set_index("architecture")
wc = pd.read_csv(T / "placement_noise_effect.csv")
CA = 1.6


def mm(x):  # hartree -> mHa string
    return f"{x*1000:.2f}"


# ---------- tables ----------
req = ab[ab.role == "required"]
t_ab = ["| Processor | Ansatz | Phys. qubits | Compiled CX | Depth | SWAPs | Noisy energy (Ha, mean of 10 seeds) | Error (mHa) | Std over seeds (mHa) |",
        "|---|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|"]
for _, r in req.iterrows():
    t_ab.append(f"| {r.architecture.replace('processor_', 'Processor ')} | {r.instance_id.split('_', 2)[2]} | {r.physical_qubits} | {r.cx_count} | {r.depth} | {r.swap_count} | {r.objective_or_accuracy:.6f} | {r.energy_error_mHa:.2f} | {r.noisy_energy_std*1000:.2f} |")
t_ab = "\n".join(t_ab)

t_mit = ["| Stage | Shots | Mean energy (Ha) | \\|Bias\\| (mHa) | Std over 10 seeds (mHa) | Runs within 1.6 mHa | Mean bias within 1.6 mHa? |",
         "|---|:-:|:-:|:-:|:-:|:-:|:-:|"]
for _, r in mit.iterrows():
    t_mit.append(f"| {r.stage} | {int(r.shots)} | {r.energy_hartree:.6f} | {r.energy_error_hartree*1000:.2f} | {r.std_hartree*1000:.2f} | {int(round(r.frac_runs_within_chem_acc*r.n_seeds))}/{int(r.n_seeds)} | {'yes' if r.within_chemical_accuracy else 'no'} |")
t_mit = "\n".join(t_mit)

t_rt = ["| Processor | Phys. qubits | Logical CX | Compiled CX | SWAPs | Depth |", "|---|:-:|:-:|:-:|:-:|:-:|"]
for _, r in rt.iterrows():
    t_rt.append(f"| {r.architecture.replace('processor_', 'Processor ')} | {r.physical_qubits} | {r.logical_cx} | {r.compiled_cx} | {r.swap_count} | {r.compiled_depth} |")
t_rt = "\n".join(t_rt)

t_lay = ["| Processor | Ordered placements of 2 logical qubits | Mean CX | Max CX | Zero-SWAP placements | Max distance |", "|---|:-:|:-:|:-:|:-:|:-:|"]
for a, r in lay.iterrows():
    t_lay.append(f"| {a.replace('processor_', 'Processor ')} | {int(r.placements)} | {r.mean_cx:.2f} | {int(r.max_cx)} | {r.frac_zero_swap*100:.0f}% | {int(r.max_distance)} |")
t_lay = "\n".join(t_lay)

t_wc = ["| Processor | Placement | Phys. qubits | Distance | CX | SWAPs | Depth | Error (mHa) |", "|---|---|:-:|:-:|:-:|:-:|:-:|:-:|"]
for _, r in wc.iterrows():
    t_wc.append(f"| {r.architecture.replace('processor_', 'Processor ')} | {r.placement} | ({r.phys_q0},{r.phys_q1}) | {r.distance} | {r.cx_count} | {r.swap_count} | {r.depth} | {r.noisy_error_mHa:.1f} |")
t_wc = "\n".join(t_wc)

# headline numbers
m4 = mit[mit.shots == 4096].set_index("stage")
m32 = mit[mit.shots == 32768].set_index("stage")
S_ZNE = "4. Readout+Sym+ZNE (linear 1,3)"
raw4 = m4.loc["1. Raw Noisy"].energy_error_hartree * 1000
zne4, zne4s = m4.loc[S_ZNE].energy_error_hartree * 1000, m4.loc[S_ZNE].std_hartree * 1000
zne4f = int(round(m4.loc[S_ZNE].frac_runs_within_chem_acc * 10))
zne32, zne32s = m32.loc[S_ZNE].energy_error_hartree * 1000, m32.loc[S_ZNE].std_hartree * 1000
zne32f = int(round(m32.loc[S_ZNE].frac_runs_within_chem_acc * 10))
A4, B4, D4 = (rt.set_index("architecture").loc[k] for k in ("processor_A", "processor_B", "processor_D"))

FACTS = """| Processor | Qubits | Edges | Degree (min-max) | Diameter | Independent cycles |
|---|:-:|:-:|:-:|:-:|:-:|
| A (organizer) | 5 | 4 | 1-3 | 3 | 0 (tree) |
| B (organizer) | 9 | 12 | 1-5 | 4 | 4 |
| D (team design) | 12 | 18 | 2-4 | 3 | 7 |"""

# ---------- README ----------
README = f"""# Qiskit Fall Fest 2026 - PS-C1: H2 Ground-State Energy (VQE) on Constrained Processors

**Team:** Team X - Shivansh Kumar and Aparna Bhadoriya
**Track:** Quantum Chemistry - PS-C1 (Phase 1)

We estimate the ground-state energy of H2 (STO-3G, R = 0.735 A) with VQE, run it under the organizers' fixed noise model on Processor A (5-qubit tree) and Processor B (9-qubit graph), study how coupling topology changes routing cost, and propose a 12-qubit Processor D.
**Every number below is read from `results/tables/*.csv` by `scripts/build_docs.py`.**

## 1. Headline results

- **Exact reference:** E_FCI = -1.1373060358 Ha (2-qubit parity-mapped Hamiltonian + nuclear repulsion 0.7199689944 Ha). Chemical accuracy = 1.6 mHa.
- **Ideal VQE** (statevector, COBYLA) reaches the exact energy to ~1e-7 mHa for all three ansatze.
- **Required A vs B benchmark:** the 2-qubit circuit needs only one CX, so both processors compile it identically (0 SWAPs, same CX/depth) and give the same noisy energy. The raw noisy error is about {raw4:.1f} mHa, far above chemical accuracy.
- **Error mitigation (Processor A, compact 1-parameter ansatz):** readout inversion + parity postselection + linear ZNE (scales 1 and 3, chosen in advance) lowers the mean bias to **{zne4:.2f} mHa at the organizer's 4096 shots**, but run-to-run std is **{zne4s:.1f} mHa**, so only **{zne4f}/10** individual runs are within 1.6 mHa. With 32768 shots (our extra setting) the mean bias is {zne32:.2f} mHa, std {zne32s:.2f} mHa, and {zne32f}/10 runs are within 1.6 mHa. **We therefore do not claim reliable chemical accuracy; we claim a large, measured bias reduction whose precision is shot-limited.**
- **Topology matters once the circuit is bigger.** On a 4-qubit all-pairs-style entangling circuit (a routing stress test, not a chemistry result) compiled CX / SWAPs / depth are A: {A4.compiled_cx}/{A4.swap_count}/{A4.compiled_depth}, B: {B4.compiled_cx}/{B4.swap_count}/{B4.compiled_depth}, D: {D4.compiled_cx}/{D4.swap_count}/{D4.compiled_depth}.

## 2. Problem and method

- 4 spin-orbitals -> parity mapping with Z2 tapering -> 2 qubits: H = g0 II + g1 IZ + g2 ZI + g3 ZZ + g4 XX with g = (-1.05237325, +0.39793742, -0.39793742, -0.01128010, +0.18093120) Ha (`data/processed/h2_2qubit_hamiltonian.json`).
- Ansatze: hardware-efficient RY+CX (reps 1 and 2) and a 1-parameter, 1-CX singlet ansatz (X, RY, CX).
- Noise (organizer, fixed): 1q depolarizing 0.001 on sx/x; CX depolarizing 0.005; readout flip 0.01. 4096 shots. Transpile: basis {{rz, sx, x, cx}}, `optimization_level=1`, SABRE layout/routing, `seed_transpiler=42`.
- Mitigation: (1) readout confusion-matrix inversion (20000 calibration shots); (2) parity postselection (odd parity); (3) ZNE by whole-circuit unitary folding at scales 1, 3 (primary, linear: E0 = (3E1 - E3)/2) and 5 (quadratic Richardson reported only for comparison). **The primary estimator is fixed in code (linear 1,3); it is never selected by comparing against the exact energy.**
- Statistics: mean and std over 10 seeds (42-51). Extrapolation amplifies shot noise (about 1.6x for linear, 2.3x for quadratic), which is why ZNE has a larger std than its inputs.

## 3. Processors

{FACTS}

Processor B uses the 12-edge map from the organizers' Benchmark Configuration PDF. The Giant Kit lists different maps; they are run as extra `*_giantkit_variant` rows in `AB_comparison_full.csv` and give identical results for this 2-qubit problem.

![Processor graphs](results/figures/processor_A_B_graphs.png)

## 4. Required benchmark: Processor A vs B (Phase 1 table)

{t_ab}

Both processors place the 2 logical qubits on an adjacent pair, so compiled circuits, noise streams and energies coincide. The deeper ansatz (reps 2) adds a CX and more depth and is worse (more noise, no accuracy benefit since reps 1 is already exact in the ideal case). Full table: `results/tables/AB_comparison.csv`.

![A vs B](results/figures/AB_performance.png)

## 5. Error mitigation study (Processor A)

{t_mit}

![Mitigation](results/figures/error_mitigation_comparison.png)

Reading it honestly: readout inversion and parity postselection give robust, large improvements. ZNE removes most of the remaining *bias* but multiplies the shot noise, so at 4096 shots a single run lands within 1.6 mHa only occasionally. More shots (or averaging several runs) is what makes the result usable.

## 6. Architecture analysis

**4-qubit routing stress test** (`full_entangling_4q_ansatz`: RY layer, 6 CX on pairs (0,1),(1,2),(2,3),(3,0),(0,2),(1,3), RY layer; it is *not* a UCCSD circuit):

{t_rt}

A is a tree with no cycles, so the 4-cycle and chords cannot be embedded without SWAPs; B has cycles and needs fewer; D's K4 core hosts the pattern directly.

**Placement sensitivity of the 2-qubit circuit** (all ordered physical placements, forced initial layout):

{t_lay}

{t_wc}

Placement matters: the worst placement costs several extra CX and tens of mHa. Note the best-placement rows use a forced layout and different seed handling than the Phase 1 table, so their errors differ slightly from it.

![Layout](results/figures/layout_sensitivity.png)

## 7. Processor D (team design)

- **12 qubits, 18 edges.** K4 core on qubits {{0,1,2,3}}; spokes (0,4), (1,5), (2,6), (3,7); links (4,8), (5,9), (6,10), (7,11); ring 8-9-10-11-8. Degree 2-4, diameter 3.
- **Same noise as A/B** (no favourable errors), <= 12 qubits, explicit graph, CX + 1q gates.
- **Bottleneck addressed:** the tree/sparse topologies of A and B force SWAPs for circuits whose interaction graph contains cycles or K4-like patterns. D embeds a 4-qubit clique natively: 0 SWAPs, {D4.compiled_cx} CX vs {A4.compiled_cx} on A for the stress test.
- **Limits, stated plainly:** D is tailored to the K4 stress-test circuit. Only 4 qubits are used for it (the other 8 are unused here), and the 7 independent cycles/higher degree cost hardware (crosstalk, fabrication). For the 2-qubit problem D gives no benefit (same CX/depth/energy as A and B), and across *all* 2-qubit placements its mean CX is not lower than A's or B's (table above). A fair claim is "D removes SWAPs for 4-qubit all-pairs circuits", not "D is better in general".

## 8. Repository layout

```
README.md  report.md  main.ipynb  requirements.txt
data/raw, data/processed         H2 geometry and Hamiltonians (JSON)
processors/                      A, B, D (+ Giant-Kit variants) as JSON coupling maps
src/problem.py processors.py routing.py metrics.py visualization.py
scripts/build_docs.py            regenerates README/report numbers from the CSVs
results/tables, results/figures   all generated outputs
assets/qff_phase1_latex_images   LaTeX equation PNGs (Starter Guide rendering note)
```

## 9. Reproduce

```bash
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\\Scripts\\activate
pip install -r requirements.txt                     # pinned: qiskit==2.5.2, qiskit-aer==0.17.2
python -c "from src import run_full_benchmark; run_full_benchmark()"   # ~30 s
python scripts/build_docs.py                        # refresh README/report numbers
jupyter nbconvert --execute --inplace main.ipynb
```

Results were generated with Python 3.13, Qiskit 2.5.2, Qiskit Aer 0.17.2. Other versions may change SABRE routing and sampled noise slightly.

## 10. Limitations

- Minimal STO-3G basis, one bond length (0.735 A), synthetic Aer noise (no T1/T2, crosstalk, drift).
- The 4-qubit circuit is a generic routing stress test, not a variationally optimized H2 4-qubit state.
- Readout calibration uses exact basis-state circuits (no state-prep error), and ZNE uses whole-circuit folding with only 2-3 scale factors.
- Processor D is a proposal evaluated only in simulation on one circuit family.

## 11. AI assistance disclosure

AI assistants (Claude, and earlier Gemini/Antigravity) were used to draft and refactor code, run reviews, and revise documentation. The team ran the pipeline, checked results against the CSVs and the organizer specification, and is responsible for the submission.

## License
Apache 2.0 - see [LICENSE](LICENSE).
"""

# ---------- report ----------
REPORT = f"""# PS-C1 Technical Report: H2 Ground-State Energy under Constrained Processor Topologies

**Team X** - Shivansh Kumar and Aparna Bhadoriya. Qiskit Fall Fest 2026, Phase 1.
All numbers are read from `results/tables/*.csv`.

## 1. Problem
Find the ground-state energy of H2 (STO-3G, R = 0.735 A) with VQE and quantify how noise and processor connectivity affect it. The exact total energy is E_FCI = -1.1373060358 Ha (electronic -1.8572750302 Ha + nuclear repulsion 0.7199689944 Ha). Chemical accuracy is 1.6 mHa.

## 2. Hamiltonian
Parity mapping with Z2 tapering reduces the 4 spin-orbital problem to 2 qubits: H = g0 II + g1 IZ + g2 ZI + g3 ZZ + g4 XX, g = (-1.05237325, +0.39793742, -0.39793742, -0.01128010, +0.18093120) Ha. Exact diagonalization of this 4x4 matrix reproduces E_FCI.

## 3. Ansatze and ideal VQE
HEA reps 1 (4 parameters, 1 CX), HEA reps 2 (6 parameters, 2 CX), and a singlet ansatz with 1 parameter (X, RY, CX). With COBYLA on the noiseless statevector, all three reach the exact energy to ~1e-7 mHa, so expressibility is not the bottleneck; noise is.

## 4. Required benchmark: Processor A vs Processor B
Setup: organizer noise model (1q depolarizing 0.001, CX depolarizing 0.005, readout flip 0.01), 4096 shots, basis {{rz, sx, x, cx}}, `optimization_level=1`, SABRE layout and routing, `seed_transpiler=42`, 10 seeds (42-51).

{t_ab}

Both processors offer an adjacent pair, so the 2-qubit circuit compiles identically (0 SWAPs). The A-vs-B difference is therefore invisible at this size; the raw noisy error (about {raw4:.1f} mHa) comes from gate and readout noise, not routing.

## 5. Error mitigation
Applied to the 1-parameter ansatz on Processor A. Readout inversion uses a 4x4 confusion matrix from 20000 calibration shots with simplex projection. Parity postselection keeps odd-parity outcomes. ZNE folds the whole circuit (U -> U (U^dag U)^((s-1)/2)) at scales 1, 3, 5. The **pre-specified** estimator is linear Richardson on scales 1 and 3, E0 = (3E1 - E3)/2; the quadratic (1, 3, 5) extrapolation, E0 = 15/8 E1 - 5/4 E3 + 3/8 E5, is reported only for comparison. Neither was selected using the exact energy.

{t_mit}

![Mitigation](results/figures/error_mitigation_comparison.png)

Findings. (i) Readout inversion removes most of the bias (raw {raw4:.1f} mHa to {abs(m4.loc['2. Readout Mitigated'].energy_error_hartree)*1000:.1f} mHa at 4096 shots). (ii) Postselection helps slightly. (iii) Linear ZNE reduces the mean bias to {zne4:.2f} mHa at 4096 shots, but extrapolation amplifies shot noise: the std over seeds is {zne4s:.1f} mHa and only {zne4f}/10 runs are within chemical accuracy. At 32768 shots the mean bias is {zne32:.2f} mHa with std {zne32s:.2f} mHa and {zne32f}/10 runs within 1.6 mHa. **Chemical accuracy is therefore reached on average only in the higher-shot setting, and not reliably in any single run.** Quadratic extrapolation is noisier still.

## 6. Architecture analysis
Graph facts:

{FACTS}

4-qubit routing stress test (generic ansatz with CX on (0,1),(1,2),(2,3),(3,0),(0,2),(1,3); not a chemistry-derived circuit):

{t_rt}

A tree cannot embed a 4-cycle or its chords, so SWAPs are unavoidable; each SWAP costs 3 CX, each CX carries 0.5% depolarizing error. B's cycles reduce but do not remove the overhead.

Placement sensitivity for the 2-qubit circuit (ordered physical placements, forced layout):

{t_lay}

{t_wc}

Poor placement costs up to several extra CX and tens of mHa of added error. (The "best" rows use a forced layout, so they differ slightly from the Phase 1 table.)

![Layout](results/figures/layout_sensitivity.png)

## 7. Processor D
12 qubits, 18 edges: K4 core on {{0,1,2,3}}, spokes (0,4),(1,5),(2,6),(3,7), links (4,8),(5,9),(6,10),(7,11), ring 8-9-10-11-8; degree 2-4, diameter 3; same noise as A and B. It addresses the SWAP bottleneck that A and B hit on circuits with cycles/K4 patterns: the stress test compiles to {D4.compiled_cx} CX, 0 SWAPs, depth {D4.compiled_depth}, versus {A4.compiled_cx} CX / {A4.swap_count} SWAPs / depth {A4.compiled_depth} on A and {B4.compiled_cx} / {B4.swap_count} / {B4.compiled_depth} on B.

Caveats: D is designed around that one circuit; 8 of its 12 qubits are unused in it; higher connectivity has real hardware costs not modeled here; and for the 2-qubit problem D gives no advantage and over all 2-qubit placements has no lower mean CX than A or B.

## 8. Limitations
STO-3G and a single bond length; synthetic noise only; generic (non-optimized) 4-qubit circuit; idealized readout calibration; ZNE with few scale factors; D validated in simulation on one circuit family.

## 9. Reproducibility
Qiskit 2.5.2, Qiskit Aer 0.17.2 (pinned in `requirements.txt`). `python -c "from src import run_full_benchmark; run_full_benchmark()"` then `python scripts/build_docs.py`.
"""

(ROOT / "README.md").write_text(README)
(ROOT / "report.md").write_text(REPORT)
print("docs written")
