from qiskit import QuantumCircuit

def oracle_all_ones(qc: QuantumCircuit, query_register: list[int], target: int) -> QuantumCircuit:
    """
    Oracle that flips target if all qubits in query register are |1>.

    Args:
        qc: Quantum circuit to modify
        query_register: List of qubit indices for the query register
        target: Index of target qubit

    Returns:
        Modified quantum circuit
    """
    # Multi-controlled X gate: flip target if all control qubits are |1>
    if len(query_register) == 1:
        qc.cx(query_register[0], target)
    elif len(query_register) == 2:
        qc.ccx(query_register[0], query_register[1], target)
    else:
        # For more than 2 controls, use multi-controlled X
        qc.mcx(query_register, target)

    return qc


def solve(n):
    qc = QuantumCircuit(n + 1)
    oracle_all_ones(qc, list(range(0, n)), n)
    return qc
