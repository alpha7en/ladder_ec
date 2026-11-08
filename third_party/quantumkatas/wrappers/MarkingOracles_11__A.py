from qiskit import QuantumCircuit
from typing import List

def phase_oracle(qc: QuantumCircuit, query: List[int], pattern: List[bool]) -> QuantumCircuit:
    """
    Phase oracle that applies -1 phase to the marked pattern.

    Args:
        qc: Quantum circuit to modify
        query: List of query qubit indices
        pattern: Boolean list representing the pattern to mark

    Returns:
        Modified quantum circuit
    """
    # Flip qubits where pattern is False
    for i, bit in enumerate(pattern):
        if not bit:
            qc.x(query[i])

    # Apply multi-controlled Z (phase flip on |11...1⟩)
    # MCZ = H on last qubit, MCX, H on last qubit
    if len(query) == 1:
        qc.z(query[0])
    elif len(query) == 2:
        qc.cz(query[0], query[1])
    else:
        # MCZ using MCX with H sandwich on last qubit
        qc.h(query[-1])
        qc.mcx(query[:-1], query[-1])
        qc.h(query[-1])

    # Uncompute the flips
    for i, bit in enumerate(pattern):
        if not bit:
            qc.x(query[i])

    return qc


def solve(n, a):
    qc = QuantumCircuit(n)
    phase_oracle(qc, list(range(0, n)), a)
    return qc
