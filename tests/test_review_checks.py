"""Guards of the implementation: overflow bounds (Sec. 6.2), the clean-ancilla test, the E2 bound of Table 1 and
the angle restriction (A5) of the translator."""
import os

import numpy as np
import pytest

from conftest import ROOT, GEN
from ladder_ec.automata import no_overflow, minimize, WA
from ladder_ec.check import family
from ladder_ec.translator import analyze
from ladder_ec import decide as V
from ladder_ec import plans as C
from ladder_ec import instances as SN
from ladder_ec.field import is_prime


def test_overflow_boundary_passes():
    p25 = 2 ** 25 - 39                                       # a prime below 2^25
    no_overflow(p25, 2 ** 13 - 1)
    no_overflow(p25, 2 ** 13)                                # (p - 1)^2 2^13 < 2^63


@pytest.mark.parametrize("p, D", [(2 ** 31 - 1, 4), (2 ** 26 + 15, 2 ** 11 + 1)])
def test_overflow_refused(p, D):
    with pytest.raises(AssertionError):
        no_overflow(p, D)


def test_minimize_refuses_large_prime():
    W = WA(np.array([1]), np.ones((8, 1, 1), np.int64), np.array([1]))
    with pytest.raises(AssertionError):
        minimize(3037000507, W)                              # (p - 1)^2 >= 2^63


def test_instances_refuse_large_prime(monkeypatch):
    """the direct simulation reduces products of residues below 2^62: it needs p < 2^31"""
    monkeypatch.setattr(SN, "primes_1mod", lambda m, count, rng: [2 ** 31 + 11])
    ref = (os.path.join(GEN, "ghz.py"), "ghz_ref")
    with pytest.raises(AssertionError):
        SN.check(ref, ref, 3, ["E1"], primes=1)


@pytest.mark.parametrize("a, dirty", [([1, 0, 0], True), ([1, 1, 0], False), ([0, 0, 0, 0], False),
                                      ([1, 1, 1, 0], True)])
def test_reference_dirty_at(a, dirty):
    G = family(analyze(os.path.join(GEN, "dirty.py"), "ref_parity"))
    assert (V.reference_dirty_at(8, G, lambda W: W, a) is not None) == dirty


@pytest.mark.parametrize("D", [[3, 3], [5, 20], [1, 40]])
def test_e2_product_bound_not_larger(D):
    """Table 1: for E2 the bound of the product series is at most the bound of the pair series"""
    G2 = family(analyze(os.path.join(ROOT, "data", "benchmark", "gen_b3.py"), "b3_l2_r1__ref"))
    assert C.killers(16, (G2, G2), D, "product", 4) <= C.killers(16, (G2, G2), D, True, 4)


def test_angle_residue_refused():
    d = analyze(os.path.join(GEN, "angle.py"), "rz_residue")
    assert d["status"] != "ladder" and d["reason"]["code"] == "nonconst_angle"


def test_angle_exact_accepted():
    assert analyze(os.path.join(GEN, "angle.py"), "rz_exact")["status"] == "ladder"


@pytest.mark.parametrize("m", [8, 16, 32, 48, 64, 128])
def test_step4_primes_are_prime(m):
    """step 4 takes the primes p = 1 (mod m) from 2^30 upward (Proposition 2)"""
    it = SN.primes_in_order(m)
    for _ in range(20):
        p = next(it)
        assert p % m == 1 and p >= 2 ** 30 and is_prime(p)


def test_is_prime_miller_rabin():
    """deterministic Miller-Rabin against trial division, Carmichael numbers and strong pseudoprimes to small bases"""
    trial = lambda q: q > 1 and all(q % d for d in range(2, int(q ** 0.5) + 1))
    assert all(is_prime(q) == trial(q) for q in range(20000))
    assert all(is_prime(q) == trial(q) for q in range(2 ** 30 - 2000, 2 ** 30 + 2000))
    for c in (561, 1105, 1729, 2465, 2821, 6601, 8911, 3215031751, 2152302898747, 3474749660383,
              341550071728321, 3825123056546413051):
        assert not is_prime(c)
    assert is_prime(2 ** 31 - 1) and is_prime(2 ** 61 - 1)
