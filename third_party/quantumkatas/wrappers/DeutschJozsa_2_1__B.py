from qiskit import QuantumCircuit
from typing import List

def dj_state_prep(qc: QuantumCircuit, query: List[int], answer: int) -> QuantumCircuit:
    """
    Prepare initial state for Deutsch-Jozsa algorithm.

    Args:
        qc: Quantum circuit to modify
        query: List of query qubit indices
        answer: Answer qubit index

    Returns:
        Modified quantum circuit
    """
    # Apply Hadamard to all query qubits to create equal superposition
    for q in query:
        qc.h(q)

    # Prepare |-> state on answer qubit
    qc.x(answer)  # |0> -> |1>
    qc.h(answer)  # |1> -> |->


def solve(n):
    qc = QuantumCircuit(n + 1)
    dj_state_prep(qc, list(range(0, n)), n)
    return qc
