"""Shared helpers and the task description of the own generator bank (Sec. 7.5, Table 10, row "own generators").

Conventions of the bank: little-endian qubit order as in Qiskit (qubit i is bit i of an integer); data qubits first,
ancillas last. The specification of a task is an isometry W of shape (2^d_out, 2^d_in): input on qubits
0..d_in-1, the other qubits start in |0>; output on qubits 0..d_out-1, the other qubits must return to |0>.
E1/E2: d_in = d_out = number of data qubits; E5: d_in = 0 (input |0...0>); E3: the style has extra qubits
(ancillas) or the input is restricted to a subspace.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from qiskit import QuantumCircuit


@dataclass
class Style:
    name: str
    fn: Callable[..., QuantumCircuit]
    model: str
    tags: list[str] = field(default_factory=list)
    note: str = ""


@dataclass
class Task:
    id: str
    title: str
    source: str
    equiv: str
    ptype: str
    params: list[str]
    spec: Callable[..., tuple[np.ndarray, int, int]]  # -> (W, d_in, d_out)
    grid: Callable[[], list[dict]]
    styles: list[Style]
    depth_limit: Callable[..., int] | None = None
    qcoder_range: Callable[[], list[dict]] | None = None
    note: str = ""


def mcx(qc: QuantumCircuit, controls, target, ctrl_state: str | None = None):
    """Multi-controlled X; ctrl_state is given in the order of controls (controls[0] is the first character)."""
    controls = list(controls)
    if not controls:
        qc.x(target)
        return
    if ctrl_state is not None:
        ctrl_state = ctrl_state[::-1]
    qc.mcx(controls, target, ctrl_state=ctrl_state)


def mcz_pattern(qc: QuantumCircuit, qubits, bits):
    """Phase -1 on the basis states where qubits[i] == bits[i] for all i."""
    qubits, bits = list(qubits), list(bits)
    zeros = [q for q, b in zip(qubits, bits) if b == 0]
    for q in zeros:
        qc.x(q)
    if len(qubits) == 1:
        qc.z(qubits[0])
    else:
        qc.h(qubits[-1])
        qc.mcx(qubits[:-1], qubits[-1])
        qc.h(qubits[-1])
    for q in zeros:
        qc.x(q)


def bits_of(x: int, n: int) -> list[int]:
    return [(x >> i) & 1 for i in range(n)]


def qft_textbook(qc: QuantumCircuit, qubits, swaps: bool = True, swap_by_cx: bool = False):
    """|j> -> 2^{-n/2} sum_k exp(2 pi i j k / 2^n) |k>, little-endian."""
    q = list(qubits)
    n = len(q)
    for j in reversed(range(n)):
        qc.h(q[j])
        for k in reversed(range(j)):
            qc.cp(math.pi / 2 ** (j - k), q[j], q[k])
    if swaps:
        for i in range(n // 2):
            a, b = q[i], q[n - 1 - i]
            if swap_by_cx:
                qc.cx(a, b); qc.cx(b, a); qc.cx(a, b)
            else:
                qc.swap(a, b)


def iqft_textbook(qc: QuantumCircuit, qubits):
    q = list(qubits)
    n = len(q)
    for i in range(n // 2):
        qc.swap(q[i], q[n - 1 - i])
    for j in range(n):
        for k in range(j):
            qc.cp(-math.pi / 2 ** (j - k), q[j], q[k])
        qc.h(q[j])


def dft_matrix(n: int, sign: int = +1) -> np.ndarray:
    N = 2 ** n
    j = np.arange(N)
    return np.exp(sign * 2j * np.pi * np.outer(j, j) / N) / math.sqrt(N)


def perm_matrix(n_qubits: int, f: Callable[[int], int]) -> np.ndarray:
    N = 2 ** n_qubits
    M = np.zeros((N, N), dtype=complex)
    for x in range(N):
        M[f(x), x] = 1
    return M


def state_spec(vec) -> tuple[np.ndarray, int, int]:
    v = np.asarray(vec, dtype=complex).reshape(-1, 1)
    v = v / np.linalg.norm(v)
    d_out = int(round(math.log2(v.shape[0])))
    return v, 0, d_out


def unitary_spec(M) -> tuple[np.ndarray, int, int]:
    M = np.asarray(M, dtype=complex)
    d = int(round(math.log2(M.shape[0])))
    return M, d, d
