from qiskit import QuantumCircuit

def oracle_arbitrary_pattern(qc: QuantumCircuit, query_register: list[int], target: int, pattern: list[bool]) -> QuantumCircuit:
    """
    Oracle that flips target if query register matches the given bit pattern.

    Args:
        qc: Quantum circuit to modify
        query_register: List of qubit indices for the query register
        target: Index of target qubit
        pattern: List of bools representing the target pattern (True=|1>, False=|0>)

    Returns:
        Modified quantum circuit
    """
    # Flip qubits where pattern is False (so we need |1> for all-ones check)
    for i, bit in enumerate(pattern):
        if not bit:  # If pattern bit is False (|0>)
            qc.x(query_register[i])

    # Multi-controlled X gate to flip target if all qubits are |1>
    if len(query_register) == 1:
        qc.cx(query_register[0], target)
    elif len(query_register) == 2:
        qc.ccx(query_register[0], query_register[1], target)
    else:
        qc.mcx(query_register, target)

    # Flip back the qubits we flipped earlier
    for i, bit in enumerate(pattern):
        if not bit:
            qc.x(query_register[i])


def solve(n, a):
    qc = QuantumCircuit(n + 1)
    oracle_arbitrary_pattern(qc, list(range(0, n)), n, a)
    return qc
