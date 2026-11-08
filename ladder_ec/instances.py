"""Direct simulation modulo p: the quick search and step 4 of Algorithm 1 (Sec. 6.2, Proposition 2).

The automata speak only about n >= n_0 = max(n_tr(R), n_tr(G)) (the translator's programs are valid from there
on, Thm. 3), so the sizes 1 <= n < n_0 are checked on the code itself: run both generators, build the circuits
exactly mod p (angles rationalised to multiples of pi with denominator <= 2^12, assumption (A5) of Sec. 6.1),
compare per relation for all bit strings a.
"NOT equal" is certain: reduction mod p maps zero to zero, so a non-zero residue is a non-zero difference.
exact=True (step 4): whole matrices on all inputs, and the first T primes p = 1 (mod m) from p_1 = 2^30 on.
A non-zero difference of entries (or of 2 x 2 cross products, for "up to phase") times sqrt2^g, g the sqrt2 count
of both circuits, is an algebraic integer with conjugates at most 2^{1 + g/2}, so at most phi(m)(1 + g/2)/30 of
these primes kill it (Lemmas 4 and 5); T is one more than that (Proposition 2). exact=False (quick search): two
random probe vectors and random primes from [2^30, 2^31); it only looks for counterexamples.
Domain: if the reference generator fails at n, that n is outside the task and skipped; if only the solution
fails, that is a counterexample (the solution crashes where the reference is defined).
Products are reduced mod p before every addition; p < 2^31 keeps them below 2^62 (asserted in check).
"""
import importlib.util, inspect, itertools, math, os, random, sys
from fractions import Fraction
import numpy as np

from .field import is_prime


class Unsupported(Exception):
    """the circuit uses something the simulation cannot represent exactly (a gate or an angle)"""


# ---------------------------------------------------------------- arithmetic mod p
class Ring:
    """F_p with a fixed root r of order m (r = image of zeta_m); 8 | m so that zeta_8 = r^(m/8)."""
    def __init__(self, p, root, m):
        self.p, self.z, self.m = p, root % p, m
        self.half = pow(2, p - 2, p)
        self.z8 = pow(self.z, m // 8, p)

    def arr(self, a):
        return (np.array(a, dtype=object) % self.p).astype(np.int64)

    def inv(self, x):
        return pow(int(x), self.p - 2, self.p)

    def zpow(self, k):
        return pow(self.z8, k % 8, self.p)

    def e(self, t):
        """e^{i pi t} = zeta_m^{t m / 2}"""
        u = Fraction(t) * self.m / 2
        if u.denominator != 1:
            raise Unsupported("angle outside Q(zeta_m)")
        return pow(self.z, int(u) % self.m, self.p)

    def eye(self, n):
        return np.eye(n, dtype=np.int64)


def gate(R, op, angle):
    """exact gate matrix mod p; first listed qubit = first tensor factor."""
    p, z, h = R.p, R.z8, R.half
    s2i = (z - pow(z, 3, p)) * h % p          # 1/sqrt2 = (zeta - zeta^3)/2
    i = pow(z, 2, p)
    m = lambda rows: R.arr(rows)
    if angle is not None:
        if "pi" not in angle:
            raise Unsupported("angle without an exact pi-multiple")
        if Fraction(angle.get("rational", "0")) != 0:
            raise Unsupported("angle with a rational part")
        q = Fraction(angle["pi"])
    e = R.e
    def ctrl(U):
        k = U.shape[0]
        M = R.eye(2 * k)
        M[k:, k:] = U
        return M
    one = {
        "h": lambda: m([[s2i, s2i], [s2i, -s2i]]),
        "x": lambda: m([[0, 1], [1, 0]]),
        "y": lambda: m([[0, -i], [i, 0]]),
        "z": lambda: m([[1, 0], [0, -1]]),
        "s": lambda: m([[1, 0], [0, i]]),
        "sdg": lambda: m([[1, 0], [0, -i]]),
        "t": lambda: m([[1, 0], [0, z]]),
        "tdg": lambda: m([[1, 0], [0, R.zpow(-1)]]),
        "p": lambda: m([[1, 0], [0, e(q)]]),
        "rz": lambda: m([[e(-q / 2), 0], [0, e(q / 2)]]),
        "rx": lambda: m([[(e(q / 2) + e(-q / 2)) * h, -i * (e(q / 2) - e(-q / 2)) * h * R.inv(i)],
                         [-i * (e(q / 2) - e(-q / 2)) * h * R.inv(i), (e(q / 2) + e(-q / 2)) * h]]),
        "ry": lambda: m([[(e(q / 2) + e(-q / 2)) * h, -(e(q / 2) - e(-q / 2)) * h * R.inv(i)],
                         [(e(q / 2) - e(-q / 2)) * h * R.inv(i), (e(q / 2) + e(-q / 2)) * h]]),
    }
    if op in one:
        return one[op]()
    two = {
        "cx": lambda: ctrl(gate(R, "x", None)), "cz": lambda: ctrl(gate(R, "z", None)),
        "cy": lambda: ctrl(gate(R, "y", None)), "ch": lambda: ctrl(gate(R, "h", None)),
        "swap": lambda: R.eye(4)[[0, 2, 1, 3]],
        "cp": lambda: ctrl(gate(R, "p", angle)), "crz": lambda: ctrl(gate(R, "rz", angle)),
        "rzz": lambda: m(np.diag([e(-q / 2), e(q / 2), e(q / 2), e(-q / 2)])),
        "ccx": lambda: ctrl(ctrl(gate(R, "x", None))), "ccz": lambda: ctrl(ctrl(gate(R, "z", None))),
        "cswap": lambda: ctrl(R.eye(4)[[0, 2, 1, 3]]),
    }
    if op in two:
        return two[op]()
    raise Unsupported("gate " + op)


def primes_1mod(m, count, rng):
    """`count` distinct random primes p = 1 (mod m) from [2^30, 2^31) (quick search)"""
    out = []
    while len(out) < count:
        c = rng.randrange(2**30, 2**31) // m * m + 1
        if c < 2**30 or c in out:
            continue
        if is_prime(c):
            out.append(c)
    return out


def root_of_unity(p, m):
    """r of exact multiplicative order m mod p (p = 1 mod m)."""
    fac = [q for q in range(2, m + 1) if m % q == 0 and all(q % d for d in range(2, q))]
    for g in range(2, p):
        r = pow(g, (p - 1) // m, p)
        if all(pow(r, m // q, p) != 1 for q in fac):
            return r


HALF_ANGLE = {"rz", "rx", "ry", "crz", "rzz"}           # gates whose matrix uses the half angle

SQRT2 = {"h": 1, "ch": 1, "rx": 2, "ry": 2, "crx": 2, "cry": 2}     # powers of sqrt2 that clear denominators

SKIP = {"barrier", "id", "delay"}
_mods = {}


def load_fn(path, name):
    if path not in _mods:
        spec = importlib.util.spec_from_file_location("gen_" + str(len(_mods)), path)
        mod = importlib.util.module_from_spec(spec)
        sys.path.insert(0, os.path.dirname(path))
        spec.loader.exec_module(mod)
        _mods[path] = mod
    return getattr(_mods[path], name)


def has_bits(fn):
    return len(inspect.signature(fn).parameters) >= 2


def angle_of(x):
    q = Fraction(float(x) / math.pi).limit_denominator(1 << 12)
    if abs(float(x) - math.pi * q) > 1e-12:
        raise Unsupported(f"angle {x} is not a rational multiple of pi")
    return {"pi": str(q)}


def circuit_ops(qc):
    ops = []
    for ins in qc.data:
        op = ins.operation
        if op.name in SKIP:
            continue
        qs = [qc.find_bit(q).index for q in ins.qubits]
        ang = angle_of(op.params[0]) if op.params else None
        ops.append((op.name, ang, qs))
    return ops, qc.num_qubits


def field_of(ops_list):
    """smallest m (multiple of 8) whose roots of unity carry all angles; rz/rx/ry/crz/rzz use the half angle."""
    m = 8
    for ops in ops_list:
        for name, ang, _ in ops:
            if ang is None:
                continue
            q = Fraction(ang["pi"]) / (2 if name in HALF_ANGLE else 1)
            d = (q / 2).denominator                 # e^{i pi q} = zeta_m^{q m / 2}
            m = m * d // math.gcd(m, d)
    return m


USED_P = [0]                                   # largest prime used here (quick search, step 4)
_row_cache = {}


def apply_cols(R, ops, k, V):
    """U V for the circuit unitary U (qubit 0 = first tensor factor) and columns V (2^k x c), gate by gate:
    O(2^k c) per gate.
    The gate axes are moved to the front as views (no transpose copy); each output slice is a sum of input slices
    over the non-zero entries of the gate row, entrywise in int64 (each product < p^2 < 2^62, reduced before adding);
    an entry 1 is a plain copy, so permutation and diagonal gates cost no multiplication."""
    p, c = R.p, V.shape[1]
    X = (V % p).reshape([2] * k + [c])
    for name, ang, qs in ops:
        key = (R.p, name, None if ang is None else ang["pi"])
        if key not in _row_cache:
            G = gate(R, name, ang)
            m = G.shape[0].bit_length() - 1
            rows = [[(tuple((j >> (m - 1 - b)) & 1 for b in range(m)), int(G[i, j])) for j in np.flatnonzero(G[i])]
                    for i in range(G.shape[0])]
            _row_cache[key] = [(tuple((i >> (m - 1 - b)) & 1 for b in range(m)), r) for i, r in enumerate(rows)]
        Y = np.moveaxis(X, qs, range(len(qs)))
        Zf = np.empty_like(X)
        Z = np.moveaxis(Zf, qs, range(len(qs)))
        for i, terms in _row_cache[key]:
            acc = None
            for j, g in terms:
                t = Y[j] if g == 1 else Y[j] * g % p
                acc = t.copy() if acc is None else (acc + t) % p
            Z[i] = 0 if acc is None else acc
        X = Zf
    return np.ascontiguousarray(X.reshape(2**k, c))


def _prop(R, A, B):
    """A = c B for some c (B non-zero), exact mod p."""
    nz = np.flatnonzero(B)
    if len(nz) == 0:
        return not np.any(A)
    k = nz[0]
    c = int(A.flat[k]) * R.inv(int(B.flat[k])) % R.p
    return bool(np.array_equal(A, B * c % R.p))


def primes_in_order(m):
    p = 2 ** 30 // m * m + 1
    while True:
        if p >= 2 ** 30 and is_prime(p):
            yield p
        p += m


def inputs(size, rng, exact, p):
    """all basis inputs (exact) or two random probe vectors: a non-zero M has M r = 0 with probability <= 1/p"""
    if exact:
        return np.eye(size, dtype=np.int64)
    return np.stack([rng.integers(0, p, size=size, dtype=np.int64) for _ in range(2)], axis=1)


def data_block(R, ops, k, d, P):
    """images of the data inputs P (columns) with the ancillas, the last k - d qubits (lowest-order index bits),
    in |0>: the block at ancilla outputs 0, and whether all other outputs are 0 (the circuit returns its ancillas on P)"""
    e = k - d
    V = np.zeros((2**k, P.shape[1]), dtype=np.int64)
    V[::2**e, :] = P
    Y = apply_cols(R, ops, k, V)
    return Y[::2**e, :], not np.any(np.delete(Y, np.s_[::2**e], axis=0))


def relations(R, cG, cR, kinds, free_first, rng, n, exact=False):
    """dict kind -> bool, or None for E3 / E3phase when the reference does not return its ancillas (the relation is undefined).
    Whole unitaries are compared on the columns of inputs(); |0...0> and the restricted inputs are exact columns.
    Ancillas are the last qubits (lowest-order index bits): qubits n, n + 1, ..."""
    (opsG, kG), (opsR, kR) = cG, cR
    p = R.p
    res = {}
    if kG == kR:
        k = kG
        P = inputs(2**k, np.random.default_rng(rng.getrandbits(64)), exact, p)
        c = P.shape[1]
        e0 = np.zeros((2**k, 1), dtype=np.int64); e0[0] = 1
        e1 = np.zeros((2**k, 1), dtype=np.int64); e1[2**(k - 1)] = 1
        V = np.concatenate([e0, P, e1], axis=1)                 # last column: the restricted input of C12
        G, Rr = apply_cols(R, opsG, k, V), apply_cols(R, opsR, k, V)
        res["E1"] = bool(np.array_equal(G[:, 1:c + 1], Rr[:, 1:c + 1]))
        res["E2"] = _prop(R, G[:, 1:c + 1], Rr[:, 1:c + 1])
        res["E5strict"] = bool(np.array_equal(G[:, 0], Rr[:, 0]))
        res["E5phase"] = _prop(R, G[:, 0], Rr[:, 0])
        if free_first:
            res["E3"] = bool(np.array_equal(G[:, [0, c + 1]], Rr[:, [0, c + 1]]))
            res["E3phase"] = _prop(R, G[:, [0, c + 1]], Rr[:, [0, c + 1]])
        if free_first or not {"E3", "E3phase"} & set(kinds):
            return {kk: v for kk, v in res.items() if kk in kinds}
    P = inputs(2**n, np.random.default_rng(rng.getrandbits(64)), exact, p)   # data qubits; any ancilla counts
    (BG, cleanG), (BR, cleanR) = data_block(R, opsG, kG, n, P), data_block(R, opsR, kR, n, P)
    res["E3"] = cleanG and bool(np.array_equal(BG, BR)) if cleanR else None
    res["E3phase"] = cleanG and _prop(R, BG, BR) if cleanR else None
    return {kk: v for kk, v in res.items() if kk in kinds}


def bit_strings(n, max_a, rng):
    """all a in [0, 2^n) or, when there are more than max_a, 0, 2^n - 1 and random others (no repeats)."""
    if max_a is None or 2**n <= max_a:
        return range(2**n)
    a = {0, 2**n - 1}
    while len(a) < max_a:
        a.add(rng.randrange(2**n))
    return sorted(a)


def check(ref, sol, n_hi, kinds, primes=3, seed=7, free_first=False, n_lo=1, max_a=None, exact=False):
    """ref, sol: (path, function). Sizes n_lo..n_hi-1. Returns dict kind -> ("equal"|"NOT equal"|"undefined", info),
    "undefined" for E3 / E3phase at the first (n, a) where the reference does not return its ancillas (checked before the solution).
    max_a: test at most max_a bit strings per n (sampled); then "equal" only means no counterexample was found.
    exact: all bit strings, all inputs, deterministic primes (module docstring); then "equal" is a proof."""
    fR, fG = load_fn(*ref), load_fn(*sol)
    bits = has_bits(fR)
    if bits != has_bits(fG):
        raise Unsupported("reference and solution take different parameters")
    cases = []
    arng = random.Random(seed + 1)
    for n in range(n_lo, n_hi):
        for a in (bit_strings(n, max_a, arng) if bits else [None]):
            args = (n, a) if bits else (n,)
            try:
                cR = circuit_ops(fR(*args))
            except Unsupported:
                raise
            except Exception:
                continue                               # outside the task
            try:
                cG = circuit_ops(fG(*args))
            except Unsupported:
                raise
            except Exception as ex:
                cases.append((n, a, None, cR, repr(ex)))
                continue
            cases.append((n, a, cG, cR, None))
    out = {k: ("equal", None) for k in kinds}
    if not cases:
        return out, []
    m = field_of([c[2][0] for c in cases if c[2]] + [c[3][0] for c in cases])
    rng = random.Random(seed)
    if exact:
        phi = sum(1 for j in range(1, m + 1) if math.gcd(j, m) == 1)
        g = max(sum(SQRT2.get(name, 0) for c in (cG, cR) if c for name, _, _ in c[0]) for _, _, cG, cR, _ in cases)
        ps = list(itertools.islice(primes_in_order(m), int(phi * (1 + g / 2) / 30) + 1))
    else:
        ps = primes_1mod(m, primes, rng)
    assert max(ps) < 2 ** 31, "apply_cols needs p < 2^31 (products < 2^62), got %d" % max(ps)   # Sec. 6.2
    USED_P[0] = max(USED_P[0], max(ps))
    Rs = [Ring(p, root_of_unity(p, m), m) for p in ps]
    ns = sorted({c[0] for c in cases})
    anc = [k for k in kinds if k in ("E3", "E3phase")] if not free_first else []
    for n, a, cG, cR, err in cases:                    # in increasing n; stop once every kind is decided
        if all(v[0] != "equal" for v in out.values()):
            break
        if err is not None:
            for R in Rs if any(out[k][0] == "equal" for k in anc) else []:
                P = inputs(2**n, np.random.default_rng(rng.getrandbits(64)), exact, R.p)
                if not data_block(R, cR[0], cR[1], n, P)[1]:
                    for k in anc:
                        if out[k][0] == "equal":
                            out[k] = ("undefined", dict(n=n, a=a, prime=R.p))
            for k in kinds:
                if out[k][0] == "equal":
                    out[k] = ("NOT equal", dict(n=n, a=a, solution_fails=err))
            continue
        for R in Rs:
            r = relations(R, cG, cR, kinds, free_first, rng, n, exact)
            for k in kinds:
                if k in r and r[k] is None and out[k][0] == "equal":
                    out[k] = ("undefined", dict(n=n, a=a, prime=R.p))
                elif k in r and not r[k] and out[k][0] == "equal":
                    out[k] = ("NOT equal", dict(n=n, a=a, prime=R.p))
    return out, ns
