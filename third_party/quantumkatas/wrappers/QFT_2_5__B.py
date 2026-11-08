from qiskit import QuantumCircuit
import numpy as np

def prepare_square_wave_signal(qc: QuantumCircuit, register: list) -> QuantumCircuit:
    """
    Prepare square wave signal using QFT.

    Args:
        qc: Quantum circuit to modify
        register: List of qubit indices (n ≥ 2)

    Returns:
        Modified quantum circuit
    """
    n = len(register)
    if n < 2:
        raise ValueError("Register must have at least 2 qubits")

    # Following the reference implementation:
    # Prepare superposition exp(-iπ/4)|2ⁿ⁻²> + exp(iπ/4)|2ⁿ⁻²+2ⁿ⁻¹>
    # This corresponds to |010...0> and |110...0>

    # Set second qubit to |1>: |010...0>
    qc.x(register[1])

    # Put first qubit in superposition: |010...0> + |110...0>
    qc.h(register[0])

    # Apply phase corrections
    # T gate adds phase exp(iπ/4)
    qc.t(register[0])

    # Apply -T to |110...0> component (when first qubit is |1>)
    qc.x(register[0])  # Flip first qubit
    qc.tdg(register[0])  # Apply T† (which is -T for our purposes)
    qc.x(register[0])  # Flip back

    # Apply QFT
    for i in range(n):
        qc.h(register[i])
        for j in range(i+1, n):
            phase = 2.0 * np.pi / (2**(j-i+1))
            qc.cp(phase, register[j], register[i])

    # Reverse the register (part of QFT)
    for i in range(n // 2):
        qc.swap(register[i], register[n - 1 - i])


def solve(n):
    qc = QuantumCircuit(n)
    prepare_square_wave_signal(qc, list(range(0, n)))
    return qc
