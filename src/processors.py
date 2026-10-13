"""Processor loading, noise models, ansatz definitions, and ideal VQE execution for PS-C1.

Contains:
- Processor specification loaders and coupling map builders
- Fixed synthetic noise model definitions
- Variational ansätze (HEA reps=1, reps=2, Compact Singlet, 4Q molecular active space)
- Ideal statevector expectation value and COBYLA optimization
"""
import json
from pathlib import Path
from typing import Dict, Tuple, List, Optional
import numpy as np
from scipy.optimize import minimize
from qiskit import QuantumCircuit, ClassicalRegister
from qiskit.circuit import ParameterVector, Parameter
from qiskit.quantum_info import Statevector
from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError
from qiskit.transpiler import CouplingMap

from .problem import hamiltonian_2q, hamiltonian_4q, NUCLEAR_REPULSION, TERMS_2Q

ROOT = Path(__file__).resolve().parent.parent


def load_processor(name_or_path: str) -> dict:
    """Load processor specification JSON from processors/ directory."""
    p = Path(name_or_path)
    if not p.suffix:
        p = ROOT / "processors" / f"{name_or_path}.json"
    with p.open("r", encoding="utf-8") as f:
        return json.load(f)


def get_coupling_map(proc: dict) -> CouplingMap:
    """Construct Qiskit CouplingMap from processor JSON."""
    edges = [tuple(e) for e in proc["coupling_map"]]
    if proc.get("bidirectional", True):
        edges += [(b, a) for a, b in edges]
    return CouplingMap(edges)


# Alias for backward compatibility
coupling_map = get_coupling_map


def get_noise_model(proc: dict) -> NoiseModel:
    """Construct Aer NoiseModel matching the fixed challenge specification:
    - 1q depolarizing = 0.001 on {sx, x}
    - 2q depolarizing = 0.005 on {cx}
    - readout bit-flip = 0.01 per measured qubit
    - virtual rz treated as noise-free
    """
    n = proc["noise"]
    nm = NoiseModel(basis_gates=proc["basis_gates"])
    nm.add_all_qubit_quantum_error(
        depolarizing_error(n["one_qubit_depolarizing"], 1), n["one_qubit_gates"]
    )
    nm.add_all_qubit_quantum_error(
        depolarizing_error(n["cx_depolarizing"], 2), ["cx"]
    )
    p = n["readout_bitflip"]
    nm.add_all_qubit_readout_error(ReadoutError([[1 - p, p], [p, 1 - p]]))
    return nm


# Alias for backward compatibility
noise_model = get_noise_model


def hea_ansatz(reps: int = 1) -> QuantumCircuit:
    """Hardware-efficient RY ansatz on 2 qubits with linear CX entangler."""
    n = 2
    th = ParameterVector("θ", n * (reps + 1))
    qc = QuantumCircuit(n, name=f"HEA_reps{reps}")
    k = 0
    for _ in range(reps):
        for q in range(n):
            qc.ry(th[k], q); k += 1
        qc.cx(0, 1)
    for q in range(n):
        qc.ry(th[k], q); k += 1
    return qc


def compact_chemistry_ansatz() -> QuantumCircuit:
    """Minimal 1-parameter symmetry-preserving chemistry ansatz.
    Directly prepares the singlet ground state in span{|01>, |10>}:
    |psi(θ)> = cos(θ/2)|01> - sin(θ/2)|10>.
    Circuit cost: 1 X gate, 1 RY gate, 1 CX gate.
    """
    th = Parameter("θ")
    qc = QuantumCircuit(2, name="Compact_Singlet_Ansatz")
    qc.x(0)
    qc.ry(-th, 1)
    qc.cx(1, 0)
    return qc


def full_entangling_4q_ansatz() -> QuantumCircuit:
    """Generic 4-qubit entangling ansatz used as a routing stress test.

    NOT a chemically derived UCCSD circuit. It applies an RY layer, six CX gates
    on the pairs (0,1),(1,2),(2,3),(3,0),(0,2),(1,3) (a K4 interaction pattern
    containing a 4-cycle and two chords), then a second RY layer. It exists only
    to expose how coupling-graph topology affects SWAP/CX overhead.
    """
    qc = QuantumCircuit(4, name="Entangling_4Q_Ansatz")
    qc.x(0); qc.x(2)
    th = ParameterVector("θ", 8)
    for i in range(4):
        qc.ry(th[i], i)
    qc.cx(0, 1)
    qc.cx(1, 2)
    qc.cx(2, 3)
    qc.cx(3, 0)
    qc.cx(0, 2)
    qc.cx(1, 3)
    for i in range(4):
        qc.ry(th[i+4], i)
    return qc


def with_measurement(qc: QuantumCircuit, basis: str = "Z") -> QuantumCircuit:
    """Append basis rotation and measurement operations."""
    out = qc.copy()
    n = qc.num_qubits
    out.add_register(ClassicalRegister(n, "c"))
    if basis == "X":
        for i in range(n):
            out.h(i)
    elif basis == "Y":
        for i in range(n):
            out.sdg(i)
            out.h(i)
    out.measure(list(range(n)), list(range(n)))
    return out


def ideal_energy(ansatz: QuantumCircuit, params: np.ndarray, qubits: int = 2) -> float:
    """Compute exact expectation value via Statevector simulation."""
    bound = ansatz.assign_parameters(params)
    sv = Statevector(bound)
    op = hamiltonian_2q() if qubits == 2 else hamiltonian_4q()
    e_elec = float(np.real(sv.expectation_value(op)))
    return e_elec + NUCLEAR_REPULSION


def run_vqe_ideal(ansatz: QuantumCircuit, qubits: int = 2, seed: int = 42, maxiter: int = 500):
    """Run ideal statevector VQE with COBYLA classical optimizer."""
    rng = np.random.default_rng(seed)
    x0 = rng.uniform(-np.pi, np.pi, ansatz.num_parameters)
    if ansatz.name == "Compact_Singlet_Ansatz":
        x0 = np.array([0.22])
    history = []
    
    def objective(x):
        e = ideal_energy(ansatz, x, qubits=qubits)
        history.append(e)
        return e

    res = minimize(objective, x0, method="COBYLA", options={"maxiter": maxiter, "rhobeg": 0.5})
    return res.x, res.fun, history
