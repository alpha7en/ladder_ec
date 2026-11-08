from qiskit import QuantumCircuit
from typing import List

def basic_marking_oracle(qc: QuantumCircuit, query: List[int], target: int) -> QuantumCircuit:
    """
    Oracle that marks the all-ones state |11...1⟩.

    Args:
        qc: Quantum circuit to modify
        query: List of query qubit indices
        target: Target qubit index to flip

    Returns:
        Modified quantum circuit
    """
    # Multi-controlled X gate flips target when all query qubits are |1⟩
    if len(query) == 1:
        qc.cx(query[0], target)
    elif len(query) == 2:
        qc.ccx(query[0], query[1], target)
    else:
        qc.mcx(query, target)

    return qc


def solve(n):
    qc = QuantumCircuit(n + 1)
    basic_marking_oracle(qc, list(range(0, n)), n)
    return qc
