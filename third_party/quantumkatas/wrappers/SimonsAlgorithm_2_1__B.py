from qiskit import QuantumCircuit

def sa_state_prep(qc: QuantumCircuit, query: list[int]) -> QuantumCircuit:
    """
    Prepare equal superposition state for Simon's algorithm.

    Args:
        qc: Quantum circuit to modify
        query: List of query register qubit indices

    Returns:
        Modified quantum circuit with superposition state
    """
    # Apply Hadamard gate to each qubit to create superposition
    for qubit in query:
        qc.h(qubit)


def solve(n):
    qc = QuantumCircuit(n)
    sa_state_prep(qc, list(range(0, n)))
    return qc
