from qiskit import QuantumCircuit
from typing import List

def parity_oracle(qc: QuantumCircuit, query: List[int], target: int) -> QuantumCircuit:
    """
    Oracle that marks states with odd parity.

    Args:
        qc: Quantum circuit to modify
        query: List of query qubit indices
        target: Target qubit index

    Returns:
        Modified quantum circuit
    """
    # XOR all query qubits into target using CNOT gates
    for q in query:
        qc.cx(q, target)


def solve(n):
    qc = QuantumCircuit(n + 1)
    parity_oracle(qc, list(range(0, n)), n)
    return qc
