from qiskit import QuantumCircuit
from typing import List

def multi_target_oracle(qc: QuantumCircuit, query: List[int], target: int, pattern: List[bool]) -> QuantumCircuit:
    """
    Oracle that marks a specific bit pattern.

    Args:
        qc: Quantum circuit to modify
        query: List of query qubit indices
        target: Target qubit index
        pattern: Boolean list representing the pattern to mark

    Returns:
        Modified quantum circuit
    """
    # Flip qubits where pattern is False (to convert to all-ones check)
    for i, bit in enumerate(pattern):
        if not bit:
            qc.x(query[i])

    # Apply multi-controlled X
    if len(query) == 1:
        qc.cx(query[0], target)
    elif len(query) == 2:
        qc.ccx(query[0], query[1], target)
    else:
        qc.mcx(query, target)

    # Uncompute the flips
    for i, bit in enumerate(pattern):
        if not bit:
            qc.x(query[i])

    return qc


def solve(n, a):
    qc = QuantumCircuit(n + 1)
    multi_target_oracle(qc, list(range(0, n)), n, a)
    return qc
