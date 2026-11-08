from qiskit import QuantumCircuit

def oracle_alternating_bits(qc: QuantumCircuit, query_qubits: list[int], target: int) -> QuantumCircuit:
    """
    Oracle that identifies alternating bit strings.
    For 2 bits: 01, 10.
    For 3 bits: 010, 101.

    Args:
        qc: Quantum circuit to modify
        query_qubits: Input qubits
        target: Target qubit

    Returns:
        Modified quantum circuit
    """
    n = len(query_qubits)

    # We want to flip target if bits alternate.
    # a0 != a1, a1 != a2, ...
    # Simple way: Compute differences into ancillas?
    # Or strict pattern matching.
    # Two patterns: 0101... and 1010...

    # Pattern 1: 0101... (start 0)
    # X on even indices. MCX. X on even indices.

    # Pattern 2: 1010... (start 1)
    # X on odd indices. MCX. X on odd indices.

    # We can use the target qubit itself.
    # Pattern 1
    for i in range(0, n, 2):
        qc.x(query_qubits[i])
    qc.mcx(query_qubits, target)
    for i in range(0, n, 2):
        qc.x(query_qubits[i])

    # Pattern 2
    for i in range(1, n, 2):
        qc.x(query_qubits[i])
    qc.mcx(query_qubits, target)
    for i in range(1, n, 2):
        qc.x(query_qubits[i])


def solve(n):
    qc = QuantumCircuit(n + 1)
    oracle_alternating_bits(qc, list(range(0, n)), n)
    return qc
