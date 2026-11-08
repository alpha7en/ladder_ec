from qiskit import QuantumCircuit
from typing import List

def parallel_pattern(qc: QuantumCircuit, qubits: List[int]) -> QuantumCircuit:
    """
    Apply H gates to all qubits in parallel.

    Args:
        qc: Quantum circuit to modify
        qubits: List of qubit indices

    Returns:
        Modified quantum circuit
    """
    # Parallel application of H gates
    for q in qubits:
        qc.h(q)
    return qc


def solve(n):
    qc = QuantumCircuit(n)
    parallel_pattern(qc, list(range(0, n)))
    return qc
