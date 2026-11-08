from qiskit import QuantumCircuit
from typing import List

def all_basis_vectors_superposition(qc: QuantumCircuit, qs: List[int]) -> QuantumCircuit:
    """
    Create equal superposition of all N-qubit basis vectors.

    Args:
        qc: Quantum circuit to modify
        qs: List of qubit indices

    Returns:
        Modified quantum circuit
    """
    # Apply Hadamard to all qubits
    for q in qs:
        qc.h(q)

    return qc


def solve(n):
    qc = QuantumCircuit(n)
    all_basis_vectors_superposition(qc, list(range(0, n)))
    return qc
