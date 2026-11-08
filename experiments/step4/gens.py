"""Step-4 experiment (Table 8): correct solutions with a constant K. Solution = reference + an extra loop from K
whose body cancels (a block and its inverse), like the mutation C6_K of the benchmark, but without an error.
Written by make_gens.py; do not edit by hand."""
import math
from qiskit import QuantumCircuit


def cx_ref(n):
    qc = QuantumCircuit(n)
    for i in range(0, n - 1):
        qc.cx(i, i + 1)
    return qc


def cx_k5(n):
    qc = QuantumCircuit(n)
    for i in range(0, n - 1):
        qc.cx(i, i + 1)
    for i in range(5, n - 1):
        qc.cx(i, i + 1)
        qc.cx(i, i + 1)
    return qc


def cx_k8(n):
    qc = QuantumCircuit(n)
    for i in range(0, n - 1):
        qc.cx(i, i + 1)
    for i in range(8, n - 1):
        qc.cx(i, i + 1)
        qc.cx(i, i + 1)
    return qc


def cx_k11(n):
    qc = QuantumCircuit(n)
    for i in range(0, n - 1):
        qc.cx(i, i + 1)
    for i in range(11, n - 1):
        qc.cx(i, i + 1)
        qc.cx(i, i + 1)
    return qc


def cx_k14(n):
    qc = QuantumCircuit(n)
    for i in range(0, n - 1):
        qc.cx(i, i + 1)
    for i in range(14, n - 1):
        qc.cx(i, i + 1)
        qc.cx(i, i + 1)
    return qc


def cx_k20(n):
    qc = QuantumCircuit(n)
    for i in range(0, n - 1):
        qc.cx(i, i + 1)
    for i in range(20, n - 1):
        qc.cx(i, i + 1)
        qc.cx(i, i + 1)
    return qc


def rzz_ref(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(0, n - 1):
        qc.rzz(math.pi / 4, i, i + 1)
    return qc


def rzz_k5(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(0, n - 1):
        qc.rzz(math.pi / 4, i, i + 1)
    for i in range(5, n - 1):
        qc.rzz(math.pi / 8, i, i + 1)
        qc.rzz(-math.pi / 8, i, i + 1)
    return qc


def rzz_k8(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(0, n - 1):
        qc.rzz(math.pi / 4, i, i + 1)
    for i in range(8, n - 1):
        qc.rzz(math.pi / 8, i, i + 1)
        qc.rzz(-math.pi / 8, i, i + 1)
    return qc


def rzz_k11(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(0, n - 1):
        qc.rzz(math.pi / 4, i, i + 1)
    for i in range(11, n - 1):
        qc.rzz(math.pi / 8, i, i + 1)
        qc.rzz(-math.pi / 8, i, i + 1)
    return qc


def rzz_k14(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(0, n - 1):
        qc.rzz(math.pi / 4, i, i + 1)
    for i in range(14, n - 1):
        qc.rzz(math.pi / 8, i, i + 1)
        qc.rzz(-math.pi / 8, i, i + 1)
    return qc


def rzz_k20(n):
    qc = QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    for i in range(0, n - 1):
        qc.rzz(math.pi / 4, i, i + 1)
    for i in range(20, n - 1):
        qc.rzz(math.pi / 8, i, i + 1)
        qc.rzz(-math.pi / 8, i, i + 1)
    return qc


def bv_ref(n, a):
    qc = QuantumCircuit(n)
    qc.x(n - 1)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        if (a >> i) & 1:
            qc.cx(i, n - 1)
    for i in range(n - 1):
        qc.h(i)
    return qc


def bv_k5(n, a):
    qc = QuantumCircuit(n)
    qc.x(n - 1)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        if (a >> i) & 1:
            qc.cx(i, n - 1)
    for i in range(n - 1):
        qc.h(i)
    for i in range(5, n - 1):
        if (a >> i) & 1:
            qc.cx(i, n - 1)
            qc.cx(i, n - 1)
    return qc


def bv_k8(n, a):
    qc = QuantumCircuit(n)
    qc.x(n - 1)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        if (a >> i) & 1:
            qc.cx(i, n - 1)
    for i in range(n - 1):
        qc.h(i)
    for i in range(8, n - 1):
        if (a >> i) & 1:
            qc.cx(i, n - 1)
            qc.cx(i, n - 1)
    return qc


def bv_k11(n, a):
    qc = QuantumCircuit(n)
    qc.x(n - 1)
    for i in range(n):
        qc.h(i)
    for i in range(n - 1):
        if (a >> i) & 1:
            qc.cx(i, n - 1)
    for i in range(n - 1):
        qc.h(i)
    for i in range(11, n - 1):
        if (a >> i) & 1:
            qc.cx(i, n - 1)
            qc.cx(i, n - 1)
    return qc

