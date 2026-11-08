from qiskit import QuantumCircuit

def twos_complement(qc: QuantumCircuit, a: list) -> QuantumCircuit:
    n = len(a)
    # One's complement: flip all bits
    for i in range(n):
        qc.x(a[i])
    # Add 1: quantum increment (process MSB to LSB to avoid corrupting controls)
    for i in range(n-1, 0, -1):
        qc.mcx(a[:i], a[i])
    qc.x(a[0])  # Always flip LSB


def solve(n):
    qc = QuantumCircuit(n)
    twos_complement(qc, list(range(0, n)))
    return qc
