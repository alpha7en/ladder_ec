from qiskit import QuantumCircuit
from typing import List

def or_oracle(qc: QuantumCircuit, query: List[int], target: int) -> QuantumCircuit:
    """
    Oracle that marks states with at least one |1⟩ (OR function).

    Args:
        qc: Quantum circuit to modify
        query: List of query qubit indices
        target: Target qubit index

    Returns:
        Modified quantum circuit
    """
    # OR = NOT(AND(NOT(inputs)))
    # Flip all inputs
    for q in query:
        qc.x(q)

    # Apply NOT of AND (flip target, then MCX, then flip back)
    qc.x(target)

    if len(query) == 1:
        qc.cx(query[0], target)
    elif len(query) == 2:
        qc.ccx(query[0], query[1], target)
    else:
        qc.mcx(query, target)

    # Uncompute input flips
    for q in query:
        qc.x(q)


def solve(n):
    qc = QuantumCircuit(n + 1)
    or_oracle(qc, list(range(0, n)), n)
    return qc
