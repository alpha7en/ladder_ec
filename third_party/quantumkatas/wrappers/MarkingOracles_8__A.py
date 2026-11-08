from qiskit import QuantumCircuit
from typing import List

def equality_oracle(qc: QuantumCircuit, reg_a: List[int], reg_b: List[int], target: int) -> QuantumCircuit:
    """
    Oracle that marks states where reg_a == reg_b.

    Args:
        qc: Quantum circuit to modify
        reg_a: First register qubit indices
        reg_b: Second register qubit indices
        target: Target qubit index

    Returns:
        Modified quantum circuit
    """
    n = len(reg_a)

    # XOR corresponding bits - if equal, result is 0
    for i in range(n):
        qc.cx(reg_a[i], reg_b[i])

    # Flip all reg_b qubits - now all should be 1 if equal
    for q in reg_b:
        qc.x(q)

    # Multi-controlled X - fires if all reg_b are 1 (meaning all pairs were equal)
    if len(reg_b) == 1:
        qc.cx(reg_b[0], target)
    elif len(reg_b) == 2:
        qc.ccx(reg_b[0], reg_b[1], target)
    else:
        qc.mcx(reg_b, target)

    # Uncompute
    for q in reg_b:
        qc.x(q)
    for i in range(n):
        qc.cx(reg_a[i], reg_b[i])

    return qc


def solve(n):
    qc = QuantumCircuit(2 * n + 1)
    equality_oracle(qc, list(range(0, n)), list(range(n, 2 * n)), 2 * n)
    return qc
