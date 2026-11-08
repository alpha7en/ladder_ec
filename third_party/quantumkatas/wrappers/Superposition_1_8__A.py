from qiskit import QuantumCircuit
from typing import List

def ghz_state(qc: QuantumCircuit, qs: List[int]) -> QuantumCircuit:
    """
    Create GHZ state on N qubits.

    Args:
        qc: Quantum circuit to modify
        qs: List of qubit indices

    Returns:
        Modified quantum circuit
    """
    # Apply Hadamard to first qubit
    qc.h(qs[0])

    # Apply CNOT from first qubit to all others
    for i in range(1, len(qs)):
        qc.cx(qs[0], qs[i])

    return qc


def solve(n):
    qc = QuantumCircuit(n)
    ghz_state(qc, list(range(0, n)))
    return qc
