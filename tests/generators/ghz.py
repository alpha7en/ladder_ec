from qiskit import QuantumCircuit


def ghz_ref(n):
    """GHZ state: H on qubit 0, then a CNOT ladder"""
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(n - 1):
        qc.cx(i, i + 1)
    return qc


def ghz_star(n):
    """the same state, all CNOTs controlled by qubit 0"""
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(1, n):
        qc.cx(0, i)
    return qc


def ghz_phase(n):
    """GHZ with a Z on the last qubit: (|0...0> - |1...1>)/sqrt2"""
    qc = QuantumCircuit(n)
    qc.h(0)
    for i in range(n - 1):
        qc.cx(i, i + 1)
    qc.z(n - 1)
    return qc
