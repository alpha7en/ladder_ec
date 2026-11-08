"""Plans: which primes make the answer "equal" sound (Sec. 5.5 of the paper, Online Resource 1 Sec. S1).

A non-zero element of Q(zeta_m) with bounded conjugates vanishes modulo only a few primes (Lemmas 4 and 5). Hence:
exact plan (Thm. S.12): the first T primes p = 1 (mod m) by size, T = (number of primes that can fail) + 1;
probabilistic plan (Thm. S.13): when T is large (sqrt2 gates inside loops), certified segment ranks and t random primes,
error at most 2^-40.
Parameters of Sec. 6.2: T_max = EXACT_MAX = 8, eps_0 = 2^-40 (target), t_max = 200, exact ranks up to RANK_MAX = 64
primes.
"""
import itertools
import math
import random
from functools import lru_cache
import numpy as np
from .field import Field, is_prime, random_prime
from .ladders import Fixed, family_wa, segment_wa
from .automata import minimize

SQRT2 = {"h": 1, "ch": 1, "rx": 2, "ry": 2, "crx": 2, "cry": 2}     # powers of sqrt2 that clear denominators


@lru_cache(None)
def n_primes(m):
    """|P_m|: number of primes p = 1 (mod m) in [2^24, 2^25) (sieve)"""
    sieve = np.ones(2 ** 25, bool)
    sieve[:2] = False
    for q in range(2, int(2 ** 12.5) + 1):
        if sieve[q]:
            sieve[q * q::q] = False
    ps = np.flatnonzero(sieve[2 ** 24:]) + 2 ** 24
    return int(np.count_nonzero(ps % m == 1))


def sqrt2_count(seg, size):
    """kappa: an upper bound on the sqrt2 denominators of a segment on `size` qubits"""
    if isinstance(seg, Fixed):
        return sum(SQRT2.get(g.op, 0) for g in seg.gates)
    return size * sum(SQRT2.get(g.op, 0) for g in seg.body)


def support_bad(seg, e, phi):
    """number of primes >= 2^24 that can make a non-zero entry of a block (or of the fixed gates) vanish.
    An entry times sqrt2^g is an algebraic integer with all conjugates at most 2^{g/2} in absolute value
    (Lemmas 4 and 5), so it has at most floor(phi g / 2 / 24) such primes. If no entry has one, the reachable states
    of the machine mod p are exactly those over the cyclotomic field."""
    gates = seg.gates if isinstance(seg, Fixed) else seg.body
    g = sum(SQRT2.get(x.op, 0) for x in gates)
    refs = {q for x in gates for q in x.qubits} | {r for x in gates for r, _ in x.cond}
    return 4 ** len(refs) * 2 ** len(refs) * int(phi * g / 2 // 24)


def rank_certificate(fams, m, rng, target):
    """minimal segment dimensions over random primes; returns (dims per family, probability that some of them is
    below the true rank), at most 2^-target unless hopeless. Two series, their lengths fixed before drawing:
    t0 primes for the machine sizes M (exact unless all t0 are bad for some segment), then t1 primes for
    the ranks, t1 chosen from M (Lemma S.10). The rank series stops early once every rank equals its M: then the
    ranks are certain (true rank between the two)."""
    N, phi = n_primes(m), sum(1 for k in range(1, m + 1) if math.gcd(k, m) == 1)
    segs = [(f, sg) for f in fams for sg in f.segments]
    bad = [support_bad(sg, f.extra, phi) for f, sg in segs]
    err0 = lambda t: sum((b / N) ** t for b in bad)
    t0 = 3 if max(bad) >= N / 2 else next((t for t in range(1, 200) if err0(t) <= 2.0 ** -(target + 1)), 200)
    if err0(t0) >= 1:                        # stop rule (Sec. 6.2): the bound is known before the series is drawn
        return None, 1.0
    M = [0] * len(segs)
    for _ in range(t0):
        F = Field(random_prime(m, rng), m)
        M = [max(x, len(segment_wa(F, sg, f.extra).alpha)) for x, (f, sg) in zip(M, segs)]
    q = [min(1.0, (phi * (x / 2 * math.log2(x) + x * sqrt2_count(sg, 2 * x) / 2) / 24 + 1) / N)
         for x, (f, sg) in zip(M, segs)]
    err1 = lambda t: sum(qq ** t for qq in q)
    t1 = 3 if max(q) >= 0.5 else next((t for t in range(1, 200) if err1(t) <= 2.0 ** -(target + 1)), 200)
    r = [0] * len(segs)
    for _ in range(t1):
        if r == M:
            break
        F = Field(random_prime(m, rng), m)
        r = [max(x, len(minimize(F.p, segment_wa(F, sg, f.extra)).alpha)) for x, (f, sg) in zip(r, segs)]
    k = len(fams[0].segments)
    return (r[:k], r[k:]), min(1.0, err0(t0) + (0.0 if r == M else err1(t1)))


def killers(m, fams, D, pair, L, c=1):
    """at most this many primes >= 2^24 kill the shortest witness of a non-zero difference (Lemma 4),
    if the families have automata of dimensions D over Q(zeta_m); c: factor on D (2 for first_free).
    pair: False for f_G - f_R (dimension d_G + d_R), True for the pair series (d_G^2 + d_R^2, doubled sqrt2 count),
    "product" for E2: U_G U_R^+ - c I has dimension 2 d_G d_R and a single sqrt2 count (Table 1)"""
    phi = sum(1 for k in range(1, m + 1) if math.gcd(k, m) == 1)
    if pair == "product":
        size = L + 2 * c * D[0] * c * D[1]
    else:
        size = L + sum((c * d) ** (2 if pair else 1) for d in D)
    g = (2 if pair is True else 1) * sum(sqrt2_count(sg, size) for f in fams for sg in f.segments)
    return phi * (1 + g / 2) / 24


def primes_in_order(m):
    p = 2 ** 24 // m * m + 1
    while True:
        if p >= 2 ** 24 and is_prime(p):
            yield p
        p += m


HOPELESS = float("inf")    # error bound of a probabilistic plan that cannot get below 1: the answer is "no answer"
EXACT_MAX = 8      # take the exact plan when it needs at most this many primes
RANK_MAX = 64      # primes spent on an exact segment rank; above that the machine size is the bound


def exact_rank(m, build, g, M, phi):
    """the minimal dimension r over Q(zeta_m) of the series of build(F) (a segment or a whole family: on every word an
    entry of a unitary circuit, g(n) sqrt2 powers on n qubits, at most M states over Q(zeta_m)), or M when that needs
    more than RANK_MAX primes. Reduction mod p never raises the rank, and a minor of order r is killed by at most
    b(M - 1)/24 primes (Lemma S.10), so the largest rank over the first b/24 + 1 primes in order is r
    (Corollary S.11)."""
    k = int(phi * (M / 2 * math.log2(M) + M * g(2 * M) / 2) / 24)
    if k + 1 > RANK_MAX:
        return M
    return max(len(minimize(p, build(Field(p, m))).alpha) for p in itertools.islice(primes_in_order(m), k + 1))


def plan(fams, m, pair, L, c, target, seed, exact_max=EXACT_MAX):
    """(primes for the zero test, error bound of "equal", dimension bounds D of the families over Q(zeta_m)).
    Exact (Thm. S.12): the machine sizes M over Q(zeta_m) are seen on any of bad + 1 primes, the products
    of the segment ranks (exact_rank, at most M) bound the family ranks, which exact_rank then pins down when cheap,
    and more primes than killers() leave the witness non-zero modulo one of them.
    Random (Thm. S.13): certified segment ranks from a first series, then t random primes. None when the exact plan
    does not apply and there is no prime p = 1 (mod m) in [2^24, 2^25) to draw from.
    exact_max: the exact plan is taken when it needs at most this many primes (T_max of the paper)."""
    phi = sum(1 for k in range(1, m + 1) if math.gcd(k, m) == 1)
    bad = max(support_bad(sg, f.extra, phi) for f in fams for sg in f.segments)
    if bad < exact_max:
        first = list(itertools.islice(primes_in_order(m), bad + 1))
        M = lambda f, sg: max(len(segment_wa(Field(p, m), sg, f.extra).alpha) for p in first)
        D = [math.prod(exact_rank(m, lambda F: segment_wa(F, sg, f.extra), lambda n: sqrt2_count(sg, n), M(f, sg), phi)
                       for sg in f.segments) for f in fams]
        D = [exact_rank(m, lambda F: family_wa(F, f), lambda n: sum(sqrt2_count(sg, n) for sg in f.segments), d, phi)
             for f, d in zip(fams, D)]
        T = max(bad, int(killers(m, fams, D, pair, L, c))) + 1
        if T <= exact_max:
            return list(itertools.islice(primes_in_order(m), T)), 0.0, D
    if n_primes(m) == 0:
        return None
    rng = random.Random(seed)
    dims, err = rank_certificate(fams, m, rng, target + 1)
    if err >= 1:                             # stop rule: eps >= 1 whatever follows, no third series
        return [], HOPELESS, None
    D = [math.prod(d) for d in dims]
    q0 = min(1.0, (killers(m, fams, D, pair, L, c) + 1) / n_primes(m))
    t = 1 if q0 >= 1 or q0 < 2.0 ** -(target + 1) else min(200, math.ceil((target + 1) / -math.log2(q0)))
    if err + q0 ** t >= 1:                   # the same before the zero-test series: its length t is fixed already
        return [], HOPELESS, D
    return [random_prime(m, rng) for _ in range(t)], err + q0 ** t, D
