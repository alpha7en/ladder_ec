"""Steps 2-3 of Algorithm 1: equivalence of two ladder programs for all n >= n_min and all bit strings a.

Relation keys (paper symbols in brackets, Definition 2):
E1  [equiv]        U_G = U_R;
E2  [simeq]        U_G = c U_R with |c| = 1 (c may depend on n and a);
E5  [equiv_0]      U_G|0> = U_R|0>;  E5p [simeq_0]  the same up to a phase;
E3  [equiv_anc]    ancillas (the last e qubits) start in |0> and must return to |0>: R(psi (x) 0) = (V psi) (x) 0
                   and G(psi (x) 0) = (V psi) (x) 0, the ancilla counts of G and R may differ;
E3p [simeq_anc]    the same up to a phase.
E2: for unitaries U_G = c U_R iff U_G U_R^+ = c I, so E2 is E1 for the series of U_G U_R^+ against
c(a) delta(x, y), c = <0|U_G U_R^+|0> (Table 1). This product cancels what G and R share: its automaton is small
even when G and R are not (an extra loop from qubit 20 on a Clifford+T block: 26 states against 578 for the pair
series f conj f of the solution). E5p and E3p use the pair series f(a,x,y) conj f(a,x',y') (Lemma S.8).

A difference series that is non-zero modulo p is non-zero over C, so "not equal" comes with an exact
counterexample (n, a, x, y). "Equal" means: zero modulo every prime of a plan (plans.py, Thm. 2).
If few primes suffice, the plan is the first T primes p = 1 (mod m) and "equal" is certain (Thm. S.12);
otherwise t random primes with error at most 2^-40 (Thm. S.13). Sizes n < n_min are outside the programs and
are checked on the code itself (step 4, instances.py).
"""
import itertools
import math
import secrets
import numpy as np
from .field import Field, field_order
from .plans import EXACT_MAX, HOPELESS, plan, primes_in_order, sqrt2_count
from .ladders import Fixed, family_wa, segment_wa
from .automata import WA, _add, no_overflow, compose, difference, pair, letters_subset, nonzero_word, minimize

NO_EPS = "error bound >= 1: the probabilistic plan cannot certify (stop rule)"
Y0 = [0, 2, 4, 6]                                            # letters (a, x, y = 0)
NO_PRIMES = "the exact plan does not apply and there is no prime p = 1 (mod m) in [2^24, 2^25) for the random one"
Y0_PAIR = [16 * a + 8 * x + 2 * x2 for a in (0, 1) for x in (0, 1) for x2 in (0, 1)]


def angles(fam):
    return [g.angle for s in fam.segments for g in (s.gates if isinstance(s, Fixed) else s.body) if g.angle is not None]


def dagger(W):
    """letters (a, x, y) -> (a, y, x); with W built over conj roots, the series of U^+"""
    return letters_subset(W, [4 * a + 2 * y + x for a in (0, 1) for x in (0, 1) for y in (0, 1)])


def scalar(H):
    """the series H(a, 0, 0) delta(x, y): letters (a, x, x) carry the matrices of (a, 0, 0), letters x != y zero"""
    return WA(H.alpha, np.stack([H.A[4 * a] * (x == y) for a in (0, 1) for x in (0, 1) for y in (0, 1)]), H.beta)


def times_dagger(F, Fc, G, R):
    """automaton of U_G U_R^+ = G_k ... G_1 R_1^+ ... R_m^+ (segments in circuit order), built from the middle out:
    G_1 R_1^+ first, then the next G on the left and the next R^+ on the right, so that segments shared by G and R
    cancel before the rest is added (Sec. 6.2). The automaton of U_G alone is never built (an extra loop after a
    star with two fixed qubits on a Clifford+T block: U_G has 550 states, U_G U_R^+ has 310)"""
    seg = lambda F, f, i: minimize(F.p, segment_wa(F, f.segments[i], f.extra))
    W = None
    for i in range(max(len(G.segments), len(R.segments))):
        if i < len(G.segments):
            W = seg(F, G, i) if W is None else minimize(F.p, compose(F.p, seg(F, G, i), W))
        if i < len(R.segments):
            W = minimize(F.p, compose(F.p, W, dagger(seg(Fc, R, i))))
    return W


def series(F, Fc, G, R, kind):
    p = F.p
    if kind == "E2":
        H = times_dagger(F, Fc, G, R)
        return difference(p, H, scalar(H))
    WG, WR = family_wa(F, G), family_wa(F, R)
    if kind == "E1":
        return difference(p, WG, WR)
    if kind == "E5":
        return difference(p, letters_subset(WG, Y0), letters_subset(WR, Y0))
    PG, PR = minimize(p, pair(p, WG, family_wa(Fc, G))), minimize(p, pair(p, WR, family_wa(Fc, R)))
    return difference(p, minimize(p, letters_subset(PG, Y0_PAIR)), minimize(p, letters_subset(PR, Y0_PAIR)))


def data(p, W, e, z):
    """the series u -> f(u (0, z_1, 0) ... (0, z_e, 0)) on the data qubits: ancilla inputs 0, outputs z"""
    beta = W.beta
    no_overflow(p, len(beta))
    for b in reversed(z):
        beta = W.A[2 * b] @ beta % p
    return minimize(p, WA(W.alpha, W.A, beta))


def first_free(W):
    """the series on a restricted input (qubit 0 arbitrary, the rest |0>): y free on the first
    letter, y = 0 on the others. Two copies of the states, the second one entered after the first letter."""
    D = len(W.alpha)
    A = np.zeros((len(W.A), 2 * D, 2 * D), np.int64)
    A[:, :D, D:] = W.A
    A[0::2, D:, D:] = W.A[0::2]                              # letters with y = 0
    return WA(np.concatenate([W.alpha, 0 * W.alpha]), A, np.concatenate([W.beta, W.beta]))


def reference_dirty_at(m, R, cut, a):
    """exact test that the reference returns its ancillas at this (n, a): the series f_R^(z), z != 0, summed over no letters but
    restricted to the bits a (n = len(a)), is zero over Q(zeta_m). The span of alpha A_u over words u with these
    bits has dimension <= d_R, so one prime costs O(n |Sigma| d^3). An entry times sqrt2^g, g the sqrt2 count of R on
    n + e qubits, is killed by at most phi g / 2 / 24 primes >= 2^24 (Lemmas 4 and 5), so one more prime in order
    decides (Sec. 6.2).
    Returns the first z with a non-zero entry, or None."""
    phi = sum(1 for k in range(1, m + 1) if math.gcd(k, m) == 1)
    g = sum(sqrt2_count(sg, len(a) + R.extra) for sg in R.segments)
    for p in itertools.islice(primes_in_order(m), int(phi * (1 + g / 2) / 24) + 1):
        W = cut(family_wa(Field(p, m), R))
        for z in itertools.product((0, 1), repeat=R.extra):
            if not any(z):
                continue
            Dz = data(p, W, R.extra, z)
            no_overflow(p, len(Dz.alpha))
            rows, piv = np.zeros((0, len(Dz.alpha)), np.int64), []
            rows, piv, _ = _add(p, rows, piv, Dz.alpha % p)
            for ai in a:
                nxt, npiv = np.zeros((0, len(Dz.alpha)), np.int64), []
                for v in rows:
                    for l in range(4 * ai, 4 * ai + 4):
                        nxt, npiv, _ = _add(p, nxt, npiv, v @ Dz.A[l] % p)
                rows, piv = nxt, npiv
                if not len(rows):
                    break
            if len(rows) and np.any(rows @ Dz.beta % p):
                return z
    return None


def decide_e3(G, R, kind, target, seed, free_first=False, exact_max=EXACT_MAX):
    """E3 / E3p: the reference must return its ancillas modulo every prime of the plan, then the solution must return them and
    agree on the data part.
    free_first: only inputs with the data qubits 1, 2, ... in |0> (first_free)"""
    m = field_order(angles(G) + angles(R))
    L = max(G.n_min, R.n_min)
    c = 2 if free_first else 1
    cut = first_free if free_first else (lambda W: W)

    def garbage(p, W, e):                    # a word on which some ancilla output z != 0 has a non-zero value
        for z in itertools.product((0, 1), repeat=e):
            w = nonzero_word(p, data(p, W, e, z), L) if any(z) else None
            if w is not None:
                return z, w

    def reference(p):                        # step 3(a): the reference must return its ancillas
        found = garbage(p, cut(family_wa(Field(p, m), R)), R.extra)
        if found:
            return "undefined", dict(reason="the reference leaves garbage in its ancillas", z=found[0])

    def attempt(p):                          # step 3(b): the solution returns its ancillas and agrees on the data part
        F, Fc = Field(p, m), Field(p, m, conj=True)
        WG, WR = cut(family_wa(F, G)), cut(family_wa(F, R))
        found = garbage(p, WG, G.extra)
        if found:
            return "not equal", dict(n=len(found[1]), word=decode(found[1], "E1"), ancillas_out=found[0])
        DG, DR = data(p, WG, G.extra, (0,) * G.extra), data(p, WR, R.extra, (0,) * R.extra)
        if kind == "E3p":
            DG = minimize(p, pair(p, DG, data(p, cut(family_wa(Fc, G)), G.extra, (0,) * G.extra)))
            DR = minimize(p, pair(p, DR, data(p, cut(family_wa(Fc, R)), R.extra, (0,) * R.extra)))
        w = nonzero_word(p, difference(p, DG, DR), L)
        if w is not None:
            return "not equal", dict(n=len(w), word=decode(w, "E1" if kind == "E3" else "E2"))

    # an ancilla that the reference leaves non-zero is seen exactly modulo any prime, so p0 goes before the plan;
    # "not equal" needs a reference that returns its ancillas, so the solution is checked only after the reference
    # passed every prime of the plan
    p0 = next(primes_in_order(m))
    found = reference(p0)
    if found:
        return found
    planned = plan((G, R), m, kind == "E3p", L + max(G.extra, R.extra), c, target, seed, exact_max)
    if planned is None:
        return "no answer", dict(reason=NO_PRIMES)
    primes, err, dims = planned
    if err == HOPELESS:
        return "no answer", dict(reason=NO_EPS)
    NG, NR = dims
    for p in primes:
        found = p != p0 and reference(p)
        if found:
            return found
    for p in [p0] + [p for p in primes if p != p0]:
        found = attempt(p)
        if found and err > 0:                # probabilistic plan: the reference must return its ancillas at this (n, a) (Sec. 6.2)
            a = [l["a"] for l in found[1]["word"]]
            dirty = reference_dirty_at(m, R, cut, a)
            if dirty:
                return "undefined", dict(reason="the reference leaves garbage in its ancillas at the witness",
                                         n=len(a), a=a, z=dirty)
        if found:
            return found
    # words are on the data qubits; the largest series is the difference, of dimension <= d (Table 1, Lemma 3)
    d = c * (NG + NR) if kind == "E3" else c * c * (NG ** 2 + NR ** 2)
    return "equal", dict(primes=len(primes), error=err, n_from=L, n_cut=L + d)


def decode(word, kind):
    """counterexample: bits a, inputs and outputs per qubit"""
    if kind in ("E1", "E2"):
        return [dict(a=l >> 2, x=l >> 1 & 1, y=l & 1) for l in word]
    if kind == "E5":
        return [dict(a=l >> 1, x=l & 1) for l in word]
    return [dict(a=l >> 2, x=l >> 1 & 1, x2=l & 1) for l in word]


def decide(G, R, kind, seed=None, target=40, free_first=False, exact_max=EXACT_MAX):
    """seed: of the probabilistic plan; None draws a fresh one (the error bound is a probability over this draw).
    exact_max: largest exact plan taken (more primes than that: the probabilistic plan)"""
    if seed is None:
        seed = secrets.randbits(32)
    if kind in ("E3", "E3p"):
        return decide_e3(G, R, kind, target, seed, free_first, exact_max)
    assert G.extra == R.extra, "different ancilla counts: use E3"
    m = field_order(angles(G) + angles(R))
    L = max(G.n_min, R.n_min) + G.extra

    def attempt(p):
        w = nonzero_word(p, series(Field(p, m), Field(p, m, conj=True), G, R, kind), L)
        return None if w is None else ("not equal", dict(n=len(w) - G.extra, word=decode(w, kind)))

    p0 = next(primes_in_order(m))            # a witness modulo any prime is exact: one prime before the plan
    found = attempt(p0)
    if found:
        return found
    planned = plan((G, R), m, "product" if kind == "E2" else kind == "E5p", L, 1, target, seed, exact_max)
    if planned is None:
        return "no answer", dict(reason=NO_PRIMES)
    primes, err, dims = planned
    if err == HOPELESS:                      # eps >= 1 is "no answer" (Table 2); only p0 was tried
        return "no answer", dict(reason=NO_EPS)
    DG, DR = dims
    for p in primes:
        found = p != p0 and attempt(p)
        if found:
            return found
    # a non-zero difference has a witness of length <= L + d (Lemma 3): sizes up to n_cut decide (Corollary 1)
    d = {"E1": DG + DR, "E5": DG + DR, "E2": 2 * DG * DR, "E5p": DG ** 2 + DR ** 2}[kind]
    return "equal", dict(primes=len(primes), error=err, n_from=L - G.extra, n_cut=L + d - G.extra)
