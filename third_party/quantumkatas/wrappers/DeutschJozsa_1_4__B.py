from qiskit import QuantumCircuit
from typing import List

def oracle_odd_number_of_ones(qc: QuantumCircuit, x: List[int], y: int) -> QuantumCircuit:
    """
    Oracle for parity function f(x) = 1 if x has odd number of 1s.

    Args:
        qc: Quantum circuit to modify
        x: List of input qubit indices
        y: Output qubit index

    Returns:
        Modified quantum circuit
    """
    # Apply CNOT from each input qubit to the output
    # This computes the XOR of all input bits
    for q in x:
        qc.cx(q, y)


def solve(n):
    qc = QuantumCircuit(n + 1)
    oracle_odd_number_of_ones(qc, list(range(0, n)), n)
    return qc
