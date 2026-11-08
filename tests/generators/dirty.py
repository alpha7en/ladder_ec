from qiskit import QuantumCircuit


def ref_parity(n, a):
    """data qubits 0..n-1, ancilla n: flipped once per bit a_i = 1, so the reference is clean iff a has even parity"""
    qc = QuantumCircuit(n + 1)
    for i in range(n):
        if (a >> i) & 1:
            qc.x(n)
    return qc
