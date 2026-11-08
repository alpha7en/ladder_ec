"""Gate entries modulo a prime (Sec. 5.4 of the paper).

Entries of the supported gates lie in Z[zeta_m, 1/2]. For a prime p = 1 (mod m) and r in F_p of order m,
zeta_m -> r is a ring homomorphism, so every construction that is polynomial in gate entries can be run in F_p
exactly. Complex conjugation corresponds to the other root r^{-1}.

Random primes are below 2^25: a product of two residues is below 2^50, so an int64 matrix product is exact as long
as the inner dimension is below 2^13. The exact plan and the exact ranks take primes upward from 2^24 without that
upper bound, so the condition D (p - 1)^2 < 2^63 is asserted before every product (automata.no_overflow, Sec. 6.2).
"""
import math
from fractions import Fraction
import numpy as np
from .automata import no_overflow


SMALL_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime(q):
    """deterministic Miller-Rabin with the first 12 primes as bases: exact for q < 3.3 * 10^24"""
    if q < 2:
        return False
    for b in SMALL_PRIMES:
        if q % b == 0:
            return q == b
    assert q < 3_317_044_064_679_887_385_961_981, "is_prime: outside the deterministic range"
    d, s = q - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for b in SMALL_PRIMES:
        x = pow(b, d, q)
        if x in (1, q - 1):
            continue
        for _ in range(s - 1):
            x = x * x % q
            if x == q - 1:
                break
        else:
            return False
    return True


def random_prime(m, rng):
    """random prime p = 1 (mod m) from [2^24, 2^25)"""
    while True:
        p = rng.randrange(2 ** 24 // m, 2 ** 25 // m) * m + 1
        if is_prime(p):
            return p


def field_order(angles):
    """smallest m (multiple of 8) such that e^{i pi t / 2} lies in Q(zeta_m) for all angles t pi"""
    m = 8
    for t in angles:
        d = 4 * Fraction(t).denominator
        m = m * d // math.gcd(m, d)
    return m


class Field:
    """F_p with a fixed root r of order m, the image of zeta_m; conj=True takes r^{-1} (complex conjugation)."""

    def __init__(self, p, m, conj=False):
        self.p, self.m = p, m
        g = 2
        while True:
            r = pow(g, (p - 1) // m, p)
            if all(pow(r, m // q, p) != 1 for q in range(2, m + 1) if m % q == 0 and is_prime(q)):
                break
            g += 1
        self.r = pow(r, p - 2, p) if conj else r

    def e(self, t):
        """e^{i pi t} for rational t"""
        u = Fraction(t) * self.m / 2
        assert u.denominator == 1, "angle outside Q(zeta_m)"
        return pow(self.r, int(u) % self.m, self.p)

    def inv(self, v):
        return pow(int(v) % self.p, self.p - 2, self.p)

    def mat(self, rows):
        return np.array([[int(v) % self.p for v in row] for row in rows], dtype=np.int64)

    def gate(self, op, t=None):
        """matrix of a gate; t is the angle divided by pi; the first listed qubit is the most significant bit"""
        e, i, h = self.e, self.e(Fraction(1, 2)), self.inv(2)
        s = self.inv(e(Fraction(1, 4)) + e(Fraction(-1, 4)))          # 1/sqrt2
        c = lambda u: [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, u[0][0], u[0][1]], [0, 0, u[1][0], u[1][1]]]
        one = {
            "h": lambda: [[s, s], [s, -s]],
            "x": lambda: [[0, 1], [1, 0]],
            "y": lambda: [[0, -i], [i, 0]],
            "z": lambda: [[1, 0], [0, -1]],
            "s": lambda: [[1, 0], [0, i]],
            "sdg": lambda: [[1, 0], [0, -i]],
            "t": lambda: [[1, 0], [0, e(Fraction(1, 4))]],
            "tdg": lambda: [[1, 0], [0, e(Fraction(-1, 4))]],
            "p": lambda: [[1, 0], [0, e(t)]],
            "rz": lambda: [[e(-t / 2), 0], [0, e(t / 2)]],
            "rx": lambda: [[(e(t / 2) + e(-t / 2)) * h, -(e(t / 2) - e(-t / 2)) * h],
                           [-(e(t / 2) - e(-t / 2)) * h, (e(t / 2) + e(-t / 2)) * h]],
            "ry": lambda: [[(e(t / 2) + e(-t / 2)) * h, -(e(t / 2) - e(-t / 2)) * h * self.inv(i)],
                           [(e(t / 2) - e(-t / 2)) * h * self.inv(i), (e(t / 2) + e(-t / 2)) * h]],
        }
        if op in one:
            return self.mat(one[op]())
        if op[0] == "c" and op[1:] in one and op != "cswap" and op[:2] != "cc":
            return self.mat(c(one[op[1:]]()))
        if op == "swap":
            return self.mat([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
        if op == "rzz":
            return self.mat(np.diag([e(-t / 2), e(t / 2), e(t / 2), e(-t / 2)]))
        if op in ("ccx", "ccz", "cswap"):
            inner = self.gate({"ccx": "cx", "ccz": "cz", "cswap": "swap"}[op])
            m = np.eye(8, dtype=np.int64)
            m[4:, 4:] = inner
            return m
        raise ValueError(op)


def apply(F, U, wires, psi):
    """apply the gate U to the given wires of psi, an array of shape (2,)*nq + (columns,)"""
    k = len(wires)
    no_overflow(F.p, 2 ** k)
    out = np.tensordot(U.reshape((2,) * 2 * k), psi, axes=(list(range(k, 2 * k)), list(wires))) % F.p
    return np.moveaxis(out, list(range(k)), list(wires))
