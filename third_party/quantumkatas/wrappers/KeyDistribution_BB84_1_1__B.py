from qiskit import QuantumCircuit

def diagonal_basis(qc: QuantumCircuit, qubits: list[int]) -> QuantumCircuit:
    """
    Convert qubits from computational basis to diagonal basis.

    Args:
        qc: Quantum circuit to modify
        qubits: List of qubit indices to convert

    Returns:
        Modified quantum circuit
    """
    # Apply Hadamard gate to each qubit to convert to diagonal basis
    # H|0> = |+> = (|0> + |1>)/sqrt2
    # H|1> = |-> = (|0> - |1>)/sqrt2
    for qubit in qubits:
        qc.h(qubit)


def solve(n):
    qc = QuantumCircuit(n)
    diagonal_basis(qc, list(range(0, n)))
    return qc
