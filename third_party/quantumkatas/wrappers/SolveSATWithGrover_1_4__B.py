from qiskit import QuantumCircuit, QuantumRegister

def oracle_alternating_bits(qc: QuantumCircuit, query_qubits: list[int], target: int) -> QuantumCircuit:
    """
    Alternating bits oracle: marks states where adjacent bits are different.

    Args:
        qc: Quantum circuit to modify
        query_qubits: List of input qubit indices
        target: Target qubit index

    Returns:
        Modified quantum circuit
    """
    n = len(query_qubits)
    if n < 2:
        return qc

    # Create auxiliary qubits to store XOR results for each adjacent pair
    # We need to track these to properly implement the oracle
    aux_qubits = list(range(max(query_qubits + [target]) + 1,
                           max(query_qubits + [target]) + 1 + (n-1)))

    # Add auxiliary qubits to circuit
    while qc.num_qubits < max(aux_qubits) + 1:
        qc.add_register(QuantumRegister(1))

    # Compute XOR for each adjacent pair
    for i in range(n-1):
        # XOR of adjacent bits stored in auxiliary qubit
        qc.cx(query_qubits[i], aux_qubits[i])
        qc.cx(query_qubits[i+1], aux_qubits[i])

    # AND all the XOR results - target flips if ALL adjacent pairs are different
    if n == 2:
        qc.cx(aux_qubits[0], target)
    elif n == 3:
        qc.ccx(aux_qubits[0], aux_qubits[1], target)
    else:
        qc.mcx(aux_qubits[:n-1], target)

    # Uncompute auxiliary qubits
    for i in range(n-1):
        qc.cx(query_qubits[i+1], aux_qubits[i])
        qc.cx(query_qubits[i], aux_qubits[i])


def solve(n):
    qc = QuantumCircuit(n + 1)
    oracle_alternating_bits(qc, list(range(0, n)), n)
    return qc
