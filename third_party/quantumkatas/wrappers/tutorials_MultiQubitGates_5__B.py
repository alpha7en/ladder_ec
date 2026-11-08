from qiskit import QuantumCircuit
from typing import List

def multi_controls(qc: QuantumCircuit, controls: List[int], target: int, control_bits: List[bool]) -> QuantumCircuit:
    """
    Apply controlled X gate with multiple controls matching a pattern.

    Args:
        qc: Quantum circuit to modify
        controls: List of control qubit indices
        target: Target qubit index
        control_bits: Pattern that controls must match

    Returns:
        Modified quantum circuit
    """
    # Use multi-controlled X with state preparation
    # First flip qubits where control_bits is False to invert the control
    for i, bit in enumerate(control_bits):
        if not bit:
            qc.x(controls[i])

    # Apply multi-controlled X
    qc.mcx(controls, target)

    # Flip back the inverted controls
    for i, bit in enumerate(control_bits):
        if not bit:
            qc.x(controls[i])


def solve(n, a):
    qc = QuantumCircuit(n + 1)
    multi_controls(qc, list(range(0, n)), n, a)
    return qc
