from qiskit import QuantumCircuit
import numpy as np

def binary_fraction_quantum(qc: QuantumCircuit, q: int, j_register: list) -> QuantumCircuit:
    """
    Apply controlled phase based on quantum register to qubit q.

    Args:
        qc: Quantum circuit to modify
        q: Index of target qubit
        j_register: List of control qubit indices representing binary fraction

    Returns:
        Modified quantum circuit
    """
    # Apply controlled phase rotation for each qubit in the register
    for ind, control_qubit in enumerate(j_register):
        # Apply controlled phase 2π / 2^(ind+1)
        phase = 2.0 * np.pi / (2**(ind + 1))
        qc.cp(phase, control_qubit, q)


def solve(n):
    qc = QuantumCircuit(n + 1)
    binary_fraction_quantum(qc, n, list(range(0, n)))
    return qc
