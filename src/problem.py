"""PS-C1: Molecular Ground-State Energy (H2) Problem Formulation.

Constructs 2-qubit parity-reduced and 4-qubit Jordan-Wigner molecular Hamiltonians
for H2 (STO-3G basis, R = 0.735 Angstrom), calculates exact classical ground-state energies,
and enforces the chemical accuracy threshold (~0.0016 Hartree).
"""
import json
from pathlib import Path
import numpy as np
from qiskit.quantum_info import SparsePauliOp

ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = ROOT / "data" / "raw"
DATA_PROCESSED = ROOT / "data" / "processed"

# Load raw molecular parameters
with open(DATA_RAW / "h2_geometry_sto3g.json", "r", encoding="utf-8") as f:
    RAW_CONFIG = json.load(f)

BOND_LENGTH_ANGSTROM = RAW_CONFIG["bond_length_angstrom"]
NUCLEAR_REPULSION = RAW_CONFIG["nuclear_repulsion_hartree"]
CHEMICAL_ACCURACY = RAW_CONFIG["chemical_accuracy_threshold_hartree"]
EXACT_REFERENCE_ENERGY = RAW_CONFIG["fci_exact_ground_energy_hartree"]

# Load processed 2-qubit Hamiltonian
with open(DATA_PROCESSED / "h2_2qubit_hamiltonian.json", "r", encoding="utf-8") as f:
    HAM_2Q_DATA = json.load(f)

TERMS_2Q = HAM_2Q_DATA["pauli_terms"]

# Load processed 4-qubit Hamiltonian
with open(DATA_PROCESSED / "h2_4qubit_jordan_wigner_hamiltonian.json", "r", encoding="utf-8") as f:
    HAM_4Q_DATA = json.load(f)

TERMS_4Q = HAM_4Q_DATA["pauli_terms"]


def hamiltonian_2q() -> SparsePauliOp:
    """Return the 2-qubit parity-reduced electronic Hamiltonian as a SparsePauliOp."""
    return SparsePauliOp.from_list(list(TERMS_2Q.items()))


def hamiltonian_4q() -> SparsePauliOp:
    """Return the 4-qubit Jordan-Wigner electronic Hamiltonian as a SparsePauliOp."""
    return SparsePauliOp.from_list(list(TERMS_4Q.items()))


def exact_ground_energy(qubits: int = 2, include_nuclear: bool = True) -> float:
    """Compute the exact ground-state energy by classical matrix diagonalization.
    
    Args:
        qubits: 2 for reduced parity Hamiltonian, 4 for full Jordan-Wigner.
        include_nuclear: Whether to add the nuclear repulsion constant.
        
    Returns:
        Exact ground-state energy in Hartree.
    """
    op = hamiltonian_2q() if qubits == 2 else hamiltonian_4q()
    e_elec = float(np.linalg.eigvalsh(op.to_matrix())[0])
    return e_elec + (NUCLEAR_REPULSION if include_nuclear else 0.0)


if __name__ == "__main__":
    e2 = exact_ground_energy(2)
    print(f"H2 (2-qubit parity) Exact Ground Energy: {e2:.10f} Ha (Ref: {EXACT_REFERENCE_ENERGY:.10f} Ha)")
    print(f"Chemical Accuracy Target: {CHEMICAL_ACCURACY * 1000:.2f} mHa ({CHEMICAL_ACCURACY:.5f} Ha)")
