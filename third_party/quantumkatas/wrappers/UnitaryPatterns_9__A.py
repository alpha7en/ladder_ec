from qiskit import QuantumCircuit
from typing import List

def multi_controlled_pattern(qc: QuantumCircuit, controls: List[int], target: int) -> QuantumCircuit:
    """
    Apply multi-controlled X gate.

    Args:
        qc: Quantum circuit to modify
        controls: List of control qubit indices
        target: Target qubit index

    Returns:
        Modified quantum circuit
    """
    # Multi-controlled X gate
    if len(controls) == 1:
        qc.cx(controls[0], target)
    elif len(controls) == 2:
        qc.ccx(controls[0], controls[1], target)
    else:
        qc.mcx(controls, target)
    return qc


def solve(n):
    qc = QuantumCircuit(n + 1)
    multi_controlled_pattern(qc, list(range(0, n)), n)
    return qc
