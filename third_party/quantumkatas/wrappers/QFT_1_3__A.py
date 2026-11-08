from qiskit import QuantumCircuit
import numpy as np

def binary_fraction_classical(qc: QuantumCircuit, q: int, j: list) -> QuantumCircuit:
    """
    Apply phase based on classical binary fraction to qubit q.

    Args:
        qc: Quantum circuit to modify
        q: Index of qubit to apply phase
        j: List of bits representing binary fraction

    Returns:
        Modified quantum circuit
    """
    # Apply R1 rotation for each bit position in the binary fraction
    for ind, bit in enumerate(j):
        if bit == 1:
            # Apply phase 2π / 2^(ind+1) = π / 2^ind
            phase = 2.0 * np.pi / (2**(ind + 1))
            qc.p(phase, q)  # Phase gate equivalent to R1
    return qc


def solve(n, a):
    qc = QuantumCircuit(1)
    binary_fraction_classical(qc, 0, a)
    return qc
