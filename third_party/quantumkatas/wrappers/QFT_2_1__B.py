from qiskit import QuantumCircuit
from typing import List
from qiskit.circuit.library import QFT

def prepare_equal_superposition(qc: QuantumCircuit, register: List[int]) -> QuantumCircuit:
    """
    Prepare equal superposition of all basis states using QFT.

    Args:
        qc: Quantum circuit to modify
        register: List of qubit indices

    Returns:
        Modified quantum circuit
    """
    # QFT applied to |0...0> state gives equal superposition of all states
    qft = QFT(len(register))
    qc.append(qft, register)


def solve(n):
    qc = QuantumCircuit(n)
    prepare_equal_superposition(qc, list(range(0, n)))
    return qc
