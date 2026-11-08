from qiskit import QuantumCircuit

def prepare_equal_superposition_even_states(qc: QuantumCircuit, qubits: list) -> QuantumCircuit:
    """
    Prepare equal superposition of even states |0>, |2>, |4>...
    Even states have LSB = 0.
    So q0 (LSB) must be |0>.
    Other qubits (q1..qn) can be in |+>.
    """
    n = len(qubits)
    if n == 0: return qc

    # Apply H to all EXCEPT LSB (qubits[0])
    for i in range(1, n):
        qc.h(qubits[i])

    # LSB remains |0>.
    return qc


def solve(n):
    qc = QuantumCircuit(n)
    prepare_equal_superposition_even_states(qc, list(range(0, n)))
    return qc
