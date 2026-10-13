"""Transpilation, routing analysis, and architecture benchmarking for PS-C1.

Implements the official transpilation baseline:
- Qiskit transpile with optimization_level=1
- SABRE-based layout and routing
- Transpiler seed 42
- Target basis gates: {rz, sx, x, cx}
"""
from typing import Dict, Tuple, Optional, List
from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap

BASIS_GATES = ["rz", "sx", "x", "cx"]
SEED_TRANSPILER = 42


def compile_circuit(
    qc: QuantumCircuit,
    cm: CouplingMap,
    initial_layout: Optional[List[int]] = None,
    keep_swaps: bool = False,
    optimization_level: int = 1
) -> QuantumCircuit:
    """Compile a quantum circuit to a target coupling map under official baseline settings."""
    basis = BASIS_GATES + (["swap"] if keep_swaps else [])
    kwargs = {
        "coupling_map": cm,
        "basis_gates": basis,
        "optimization_level": optimization_level,
        "seed_transpiler": SEED_TRANSPILER,
    }
    if initial_layout is None:
        kwargs["layout_method"] = "sabre"
        kwargs["routing_method"] = "sabre"
    else:
        kwargs["initial_layout"] = initial_layout
        kwargs["routing_method"] = "sabre"

    return transpile(qc, **kwargs)


def circuit_metrics(qc: QuantumCircuit) -> Dict[str, int]:
    """Extract operational gate counts and depth metrics from a transpiled circuit."""
    ops = qc.count_ops()
    gate_ops = {k: v for k, v in ops.items() if k not in ("measure", "barrier", "reset")}
    return {
        "depth": qc.depth(lambda inst: inst.operation.name not in ("measure", "barrier")),
        "cx_count": ops.get("cx", 0),
        "swap_count": ops.get("swap", 0),
        "one_qubit_gates": sum(v for k, v in gate_ops.items() if k not in ("cx", "swap")),
        "total_gates": sum(gate_ops.values()),
    }


def analyze_routing(
    logical_qc: QuantumCircuit,
    cm: CouplingMap,
    initial_layout: Optional[List[int]] = None
) -> Tuple[Dict[str, int], QuantumCircuit]:
    """Determine routing overhead (SWAP count) and final compiled circuit metrics."""
    # Step 1: Route keeping SWAPs as explicit gates to count routing insertions
    routed_with_swaps = compile_circuit(logical_qc, cm, initial_layout, keep_swaps=True)
    swap_count = routed_with_swaps.count_ops().get("swap", 0)

    # Step 2: Fully decompose into basis {rz, sx, x, cx}
    final_compiled = compile_circuit(logical_qc, cm, initial_layout, keep_swaps=False)
    metrics = circuit_metrics(final_compiled)
    metrics["swap_count"] = swap_count

    # Determine unique physical qubits touched by operations
    used_indices = set()
    for inst in final_compiled.data:
        for q in inst.qubits:
            used_indices.add(final_compiled.find_bit(q).index)
    metrics["physical_qubits_used"] = len(used_indices)

    return metrics, final_compiled


def benchmark_4qubit_architecture_routing(procs_dict: Dict[str, dict], ansatz_4q: QuantumCircuit) -> List[dict]:
    """Run 4-qubit molecular ansatz routing comparison across Processor A, B, and D.
    
    Shows concrete architecture-dependent routing differences:
    - Processor A (tree): limited diameter requires SWAP insertions
    - Processor B (9-qubit graph with cycles): different routing path and lower depth
    - Processor D (K4 core + ring): zero SWAPs because the ansatz pattern fits the K4 core
    """
    from .processors import get_coupling_map, with_measurement
    
    logical = with_measurement(ansatz_4q, "Z")
    results = []
    
    for name, proc in procs_dict.items():
        cm = get_coupling_map(proc)
        m, transpiled = analyze_routing(logical, cm)
        results.append({
            "architecture": name,
            "logical_qubits": 4,
            "physical_qubits": proc["num_qubits"],
            "logical_cx": int(ansatz_4q.count_ops().get("cx", 0)),
            "ansatz": "full_entangling_4q (routing stress test, not UCCSD)",
            "compiled_cx": m["cx_count"],
            "swap_count": m["swap_count"],
            "compiled_depth": m["depth"],
            "one_qubit_gates": m["one_qubit_gates"],
            "physical_qubits_used": m["physical_qubits_used"]
        })
    return results
