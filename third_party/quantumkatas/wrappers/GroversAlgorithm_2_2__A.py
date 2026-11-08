import numpy as np
from qiskit import QuantumCircuit

def conditional_phase_flip(qc: QuantumCircuit, register: list[int]) -> QuantumCircuit:
    """
    Apply phase flip to all states except |0...0>.

    Args:
        qc: Quantum circuit to modify
        register: List of qubit indices

    Returns:
        Modified quantum circuit
    """
    # Flip all qubits to convert |0...0> check to |1...1> check
    for qubit in register:
        qc.x(qubit)

    # Apply multi-controlled Z gate (phase flip when all qubits are |1>)
    if len(register) == 1:
        qc.z(register[0])
    elif len(register) == 2:
        qc.cz(register[0], register[1])
    else:
        # Multi-controlled Z: use MCX with auxiliary qubits or direct MCZ
        # For simplicity, use a phase gate pattern
        qc.mcx(register[:-1], register[-1])
        qc.z(register[-1])
        qc.mcx(register[:-1], register[-1])

    # Flip qubits back
    for qubit in register:
        qc.x(qubit)

    # Add global phase correction (π phase)
    qc.global_phase += np.pi

    return qc


def solve(n):
    qc = QuantumCircuit(n)
    conditional_phase_flip(qc, list(range(0, n)))
    return qc
