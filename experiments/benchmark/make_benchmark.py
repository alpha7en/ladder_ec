"""Build the benchmark (Sec. 7.1): references, correct solutions and mutants for every task.

Tasks are grouped in classes B1-B17 (ladders built from small blocks) and K1-K10 (patterns of QuantumKatas-like
tasks). Tasks outside the class of ladder programs (a reference branching on n, nested loops, gates of growing
arity, ...) are written as ready code (Raw) and give one reference-reference pair, on which a refusal is expected.
A task is a small tree (Loop, Call, Gate, If, Raw); a mutation kind is one edit of this tree (Table 3).
The tree is printed as plain Qiskit code, which the translator reads.
The answer expected by construction goes to pairs.jsonl; ground_truth.py rechecks it by direct computation.
Usage: python -I experiments/benchmark/make_benchmark.py [OUT_DIR]   (default data/benchmark) -> gen_*.py, pairs.jsonl
"""
import json
import os
import random
import sys
from dataclasses import dataclass, field, replace

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HERE = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data", "benchmark")
REL = "data/benchmark"

# --- task tree. An index is a pair (base, offset): c constant, n is n + k, i is i + k, m is n - 1 - i + k


@dataclass
class Gate:
    op: str
    qs: list
    angle: str = None


@dataclass
class Call:
    block: str
    qs: list


@dataclass
class If:
    bit: tuple
    then: list
    other: list = field(default_factory=list)


@dataclass
class Loop:
    start: tuple
    stop: tuple
    step: int
    body: list
    rev: bool = False


@dataclass
class Raw:
    lines: list


def idx(r, var="i"):
    b, k = r
    if b == "c":
        return str(k)
    if b == "m":
        return "n - %d - %s" % (1 - k, var) if k < 1 else "n + %d - %s" % (k - 1, var) if k > 1 else "n - " + var
    head = {"n": "n", "i": var}[b]
    return head if k == 0 else "%s %s %d" % (head, "+" if k > 0 else "-", abs(k))


def lines(body, var="i", ind="    "):
    out = []
    for s in body:
        if isinstance(s, Gate):
            args = ([s.angle] if s.angle else []) + [idx(q, var) for q in s.qs]
            out.append("%sqc.%s(%s)" % (ind, s.op, ", ".join(args)))
        elif isinstance(s, Call):
            out.append("%su_%s(qc, %s)" % (ind, s.block, ", ".join(idx(q, var) for q in s.qs)))
        elif isinstance(s, If):
            e = idx(s.bit, var)
            out.append("%sif (a >> %s) & 1:" % (ind, "(%s)" % e if " " in e else e))
            out += lines(s.then, var, ind + "    ") or [ind + "    pass"]
            if s.other:
                out += [ind + "else:"] + lines(s.other, var, ind + "    ")
        elif isinstance(s, Loop):
            r = "range(%s, %s%s)" % (idx(s.start), idx(s.stop), "" if s.step == 1 else ", %d" % s.step)
            out.append("%sfor %s in %s:" % (ind, var, "reversed(%s)" % r if s.rev else r))
            out += lines(s.body, var, ind + "    ")
        else:
            out += [ind + x for x in s.lines]
    return out


def function(name, body, bits, var="i", helper=False, extra=0, qubits=None):
    args = "n, a" if bits else "n"
    new = "    qc = QuantumCircuit(%s)" % (qubits or ("n + %d" % extra if extra else "n"))
    if helper:
        return ["def %s(%s):" % (name, args), new, "    %s_layers(qc, %s)" % (name, args),
                "    return qc", "", "", "def %s_layers(qc, %s):" % (name, args)] + lines(body, var)
    return ["def %s(%s):" % (name, args), new] + lines(body, var) + ["    return qc"]


# --- blocks U. A block gate is (name, formal qubits, angle); 0 -> i, 1 -> j

def entangling(gs):
    """the operator Schmidt rank of the block is above 1: the block does not split into single-qubit gates"""
    one = {"h": np.array([[1, 1], [1, -1]]) / np.sqrt(2), "s": np.diag([1, 1j]), "t": np.diag([1, np.exp(1j * np.pi / 4)])}
    cx = {(0, 1): np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]]),   # qubit 0 least significant
          (1, 0): np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])}
    U = np.eye(4)
    for g, q, _ in gs:
        G = cx[q] if g == "cx" else np.kron(np.eye(2), one[g]) if q == (0,) else np.kron(one[g], np.eye(2))
        U = G @ U
    M = U.reshape(2, 2, 2, 2).transpose(0, 2, 1, 3).reshape(4, 4)
    return np.linalg.matrix_rank(M, tol=1e-9) > 1


def random_block(seed):
    """6 Clifford+T gates, at least two CX, no adjacent cancelling pairs, an entangling block"""
    rng = random.Random(seed)
    pool = [(g, (q,), None) for g in "hst" for q in (0, 1)] + [("cx", (0, 1), None), ("cx", (1, 0), None)]
    while True:
        gs = [rng.choice(pool) for _ in range(6)]
        pairs = any(x == y and x[0] in ("h", "cx") for x, y in zip(gs, gs[1:]))
        if sum(g[0] == "cx" for g in gs) >= 2 and not pairs and entangling(gs):
            return gs


def random_wide(seed, w, size, multi):
    """size Clifford+T gates on w qubits, at least multi two-qubit gates, no adjacent cancelling pairs"""
    rng = random.Random(seed)
    pool = [(g, (q,), None) for g in "hst" for q in range(w)] + \
           [("cx", (c, t), None) for c in range(w) for t in range(w) if c != t]
    while True:
        gs = [rng.choice(pool) for _ in range(size)]
        pairs = any(x == y and x[0] in ("h", "cx") for x, y in zip(gs, gs[1:]))
        if sum(g[0] == "cx" for g in gs) >= multi and not pairs and {q for g in gs for q in g[1]} == set(range(w)):
            return gs


BLOCKS = {
    "cx": [("cx", (0, 1), None)],
    "cz": [("cz", (0, 1), None)],
    "swap": [("swap", (0, 1), None)],
    "ch": [("ch", (0, 1), None)],
    "iswap": [("s", (0,), None), ("s", (1,), None), ("h", (0,), None), ("cx", (0, 1), None), ("cx", (1, 0), None),
              ("h", (1,), None)],
    "htcx": [("cx", (0, 1), None), ("h", (0,), None), ("t", (1,), None)],
    "r1": random_block(1), "r2": random_block(2), "r3": random_block(3),
    "ccx": [("ccx", (0, 1, 2), None)], "cswap": [("cswap", (0, 1, 2), None)],
    "q1": random_wide(1, 3, 8, 3), "q2": random_wide(2, 3, 8, 3),
    "w4": random_wide(1, 4, 10, 4),
    "maj": [("cx", (2, 1), None), ("cx", (2, 0), None), ("ccx", (0, 1, 2), None)],     # adder blocks (K4)
    "uma": [("ccx", (0, 1, 2), None), ("cx", (2, 0), None), ("cx", (0, 1), None)],
}
BLOCKS2 = ["cx", "cz", "swap", "ch", "iswap", "htcx", "r1", "r2", "r3"]   # two-qubit blocks
BLOCKS3 = ["ccx", "cswap", "q1", "q2"]                                     # three-qubit blocks
FEW = ["cx", "r1"]                         # blocks of the non-default variants (a fixed part of the list)
FEW3 = ["ccx", "q1"]
INV = {"s": "sdg", "sdg": "s", "t": "tdg", "tdg": "t"}


def inverse(gs):
    return [(INV.get(g, g), q, a and "-" + a) for g, q, a in reversed(gs)]


def decompose(gs):
    """an equal decomposition: CX by a reversed CX, CZ by CX, SWAP by three CX, S by two T"""
    out = []
    for g, q, a in gs:
        if g == "cx":
            h = [("h", (q[0],), None), ("h", (q[1],), None)]
            out += h + [("cx", (q[1], q[0]), None)] + h
        elif g == "cz":
            out += [("h", (q[1],), None), ("cx", q, None), ("h", (q[1],), None)]
        elif g == "swap":
            out += [("cx", q, None), ("cx", q[::-1], None), ("cx", q, None)]
        elif g in ("s", "sdg"):
            out += [("t" if g == "s" else "tdg", q, None)] * 2
        else:
            out.append((g, q, a))
    return out


def helper_text(name, gs):
    w = max(max(q) for g, q, a in gs) + 1 if gs else 2
    body = ["    qc.%s(%s)" % (g, ", ".join(([a] if a else []) + ["ijkl"[k] for k in q])) for g, q, a in gs]
    return ["def u_%s(qc, %s):" % (name, ", ".join("ijkl"[:w]))] + (body or ["    pass"]) + ["", ""]


# --- classes and variants. Each function: variant name -> (body for block b, whether bits a are needed)

def ladder(b, down=False, step=1, start=0):
    if down:
        return Loop(("n", -2), ("c", -1), -step, [Call(b, [("i", 0), ("i", 1)])])
    return Loop(("c", start), ("n", -1), step, [Call(b, [("i", 0), ("i", 1)])])


def star_inv(b):
    return Loop(("n", -1), ("c", 0), -1, [Call(b + "_inv", [("c", 0), ("i", 0)])])


def layer(g, step=1):
    return Loop(("c", 0), ("n", 0), step, [Gate(g, [("i", 0)])])


def bit_if(cond, b):
    call = [Call(b, [("i", 0), ("i", 1)])]
    if cond == "a_i":
        return [If(("i", 0), call)]
    if cond == "a_i1":
        return [If(("i", 1), call)]
    if cond == "not_a_i":
        return [If(("i", 0), [], call)]
    return [If(("i", 0), [If(("i", 1), [], call)], [If(("i", 1), call)])]   # a_i xor a_{i+1}


def mirror(v, b):
    """frame V, middle, V^{-1}"""
    V, mid = v.split("_")
    if V == "ladder":
        outer = [ladder(b)], [ladder(b + "_inv", down=True)]
    elif V == "star":
        outer = [Loop(("c", 1), ("n", 0), 1, [Call(b, [("c", 0), ("i", 0)])])], [star_inv(b)]
    else:
        outer = [layer(V[0])], [layer(V[0])]
    return outer[0] + {"last": [Gate("t", [("n", -1)])], "fixed": [Gate("t", [("c", 0)])], "none": []}[mid] + outer[1]


def window(b, w, step):
    """B2: a block on w neighbouring qubits"""
    return [Loop(("c", 0), ("n", 1 - w), step, [Call(b, [("i", k) for k in range(w)])])]


def brick(b, L):
    """B3: a layer on even pairs, then on odd pairs; L times"""
    one = [Loop(("c", 0), ("n", -1), 2, [Call(b, [("i", 0), ("i", 1)])]),
           Loop(("c", 1), ("n", -1), 2, [Call(b, [("i", 0), ("i", 1)])])]
    return one * L


def ring(b, where):
    close = [Call(b, [("n", -1), ("c", 0)])]
    return close + [ladder(b)] if where == "start" else [ladder(b)] + close


def swaps(v):
    """B13: permutations by SWAP"""
    shift1 = Loop(("c", 0), ("n", -1), 1, [Gate("swap", [("i", 0), ("i", 1)])])
    return {"shift1": [shift1], "shift2": [shift1, shift1],
            "alt": [Loop(("c", 0), ("n", -1), 2, [Gate("swap", [("i", 0), ("i", 1)])])]}[v]


def qft_band(k, angle):
    """K1: H on all qubits, then for each distance d <= k a phase cp on the pairs (i + d, i)"""
    body = [layer("h")]
    for d in range(1, k + 1):
        body.append(Loop(("c", 0), ("n", -d), 1, [Gate("cp", [("i", d), ("i", 0)], angle(d))]))
    return body


def oracle(v):
    """K5: a gate on qubit i if the bit a_i (or a condition on neighbouring bits) is 1"""
    what, cond = v.split("_", 1)
    g = {"x": [Gate("x", [("i", 0)])], "z": [Gate("z", [("i", 0)])], "cz": [Gate("cz", [("i", 0), ("i", 1)])],
         "anc": [Gate("cx", [("i", 0), ("n", 0)])]}[what]
    stop = ("n", -1) if what == "cz" or cond == "and" else ("n", 0)
    body = {"ai": [If(("i", 0), g)], "and": [If(("i", 0), [If(("i", 1), g)])], "not": [If(("i", 0), [], g)]}[cond]
    loop = [Loop(("c", 0), stop, 1, body)]
    if what == "anc":                               # phase oracle: ancilla in |->, returned to |0> at the end
        prep = [Gate("x", [("n", 0)]), Gate("h", [("n", 0)])]
        return prep + loop + prep[::-1], True, 1, "E3"
    return loop, True


def ansatz(v):
    """K7: a layer of rotations, then entangling by a chain or a ring; L times"""
    ent, L = v.split("_")
    one = [Loop(("c", 0), ("n", 0), 1, [Gate("ry", [("i", 0)], "math.pi / 2"), Gate("rz", [("i", 0)], "math.pi / 4")]),
           Loop(("c", 0), ("n", -1), 1, [Gate("cx", [("i", 0), ("i", 1)])])]
    if ent == "ring":
        one = one + [Gate("cx", [("n", -1), ("c", 0)])]
    return one * int(L[1:]), False


def trotter(v):
    """K8: Trotter step: exp(-i theta ZZ) on the edges of a chain (CX, RZ, CX), then an RX field; T steps"""
    edges, T = v.split("_")
    zz = [Gate("cx", [("i", 0), ("i", 1)]), Gate("rz", [("i", 1)], "math.pi / 4"), Gate("cx", [("i", 0), ("i", 1)])]
    if edges == "xx":
        hh = [Gate("h", [("i", 0)]), Gate("h", [("i", 1)])]
        zz = hh + zz + hh
    step = [Loop(("c", 0), ("n", -1), 1, zz)]
    if edges == "ring":
        step += [Gate("cx", [("n", -1), ("c", 0)]), Gate("rz", [("c", 0)], "math.pi / 4"), Gate("cx", [("n", -1), ("c", 0)])]
    step += [Loop(("c", 0), ("n", 0), 1, [Gate("rx", [("i", 0)], "math.pi / 2")])]
    return step * int(T[1:]), False


def pattern(v):
    """K9: matrix patterns (Katas UnitaryPatterns): which matrix entries are non-zero"""
    return {"diag": [layer("t")], "anti": [layer("x")],
            "blockdiag": [Loop(("c", 0), ("n", -1), 1, [Gate("h", [("i", 0)])])],
            "quarters": [Gate("h", [("n", -1)])]}[v], False


def ghz(v):
    """K6: GHZ state by copying qubit 0 with a ladder or a star; task relation: state equality up to phase"""
    copy = {"ladder": Loop(("c", 0), ("n", -1), 1, [Gate("cx", [("i", 0), ("i", 1)])]),
            "star": Loop(("c", 1), ("n", 0), 1, [Gate("cx", [("c", 0), ("i", 0)])])}[v]
    return [Gate("h", [("c", 0)]), copy], False, 0, "E5phase"


def prep(v):
    """K10: state preparation from |0...0>; task relation: state equality up to phase.
    Variants with a suffix are correct solutions specific to E5: another operator, the same state"""
    ladder_cx = Loop(("c", 0), ("n", -1), 1, [Gate("cx", [("i", 0), ("i", 1)])])
    star_cx = Loop(("c", 1), ("n", 0), 1, [Gate("cx", [("c", 0), ("i", 0)])])
    bit = lambda g, ang=None: [Loop(("c", 0), ("n", 0), 1, [If(("i", 0), [Gate(g, [("i", 0)], ang)])])]   # noqa: E731
    cat = Loop(("c", 1), ("n", 0), 1, [If(("i", 0), [Gate("cx", [("c", 0), ("i", 0)])])])
    body = {"plus": [layer("h")],
            "plus_ry": [Loop(("c", 0), ("n", 0), 1, [Gate("ry", [("i", 0)], "math.pi / 2")])],
            "minus": [Gate("x", [("c", 0)]), Gate("h", [("c", 0)]), ladder_cx],
            "minus_zfirst": [Gate("h", [("c", 0)]), Gate("z", [("c", 0)]), ladder_cx],
            "minus_zlast": [Gate("h", [("c", 0)]), ladder_cx, Gate("z", [("n", -1)])],
            "minus_star": [Gate("x", [("c", 0)]), Gate("h", [("c", 0)]), star_cx],
            "basis": bit("x"), "basis_ry": bit("ry", "math.pi"),
            "cat": [Gate("h", [("c", 0)]), cat], "cat_ry": [Gate("ry", [("c", 0)], "math.pi / 2"), cat]}[v]
    return body, v.split("_")[0] in ("basis", "cat"), 0, "E5phase"


# Tasks outside the class: ready reference code and the number of qubits of the circuit (if not n)
OUTSIDE = {
    "B8": {"p2": ["for i in range(n - 1):", "    if i % 2 == 0:", "        u_cx(qc, i, i + 1)", "    else:",
                  "        u_cz(qc, i, i + 1)"],
           "p3": ["for i in range(n - 1):", "    if i % 3 == 0:", "        u_cx(qc, i, i + 1)", "    else:",
                  "        u_cz(qc, i, i + 1)"]},
    "B10": {g: ["for i in range(n):", "    for j in range(i + 1, n):", "        qc.%s(i, j)" % g]
            for g in ("cz", "cx", "swap")},
    "B11": {"direct": ["qc.h(0)", "k = 1", "while k < n:", "    for j in range(0, min(k, n - k)):",
                       "        qc.cx(j, j + k)", "    k *= 2"],
            "reverse": ["k = 1", "while 2 * k < n:", "    k *= 2", "while k >= 1:", "    for j in range(0, min(k, n - k)):",
                        "        qc.cx(j, j + k)", "    k //= 2"]},
    "B12": {"mcx": ["qc.mcx(list(range(n - 1)), n - 1)"],
            "vchain": ["qc.ccx(0, 1, n)", "for i in range(1, n - 2):", "    qc.ccx(i + 1, n + i - 1, n + i)",
                       "qc.cx(2 * n - 3, n - 1)", "for i in reversed(range(1, n - 2)):",
                       "    qc.ccx(i + 1, n + i - 1, n + i)", "qc.ccx(0, 1, n)"],
            "gray": ["qc.mcx(list(range(n - 1)), n - 1, mode='noancilla')"]},
    "B13": {"reverse": ["for i in range(n // 2):", "    qc.swap(i, n - 1 - i)"]},
    "B16": {"inc_prefix": ["for i in reversed(range(1, n)):", "    qc.mcx(list(range(i)), i)", "qc.x(0)"],
            "dec_prefix": ["qc.x(0)", "for i in range(1, n):", "    qc.mcx(list(range(i)), i)"],
            "inc_suffix": ["for i in range(n - 1):", "    qc.mcx(list(range(i + 1, n)), i)", "qc.x(n - 1)"]},
    "B17": {"bitwise": ["for i in range(n):", "    qc.cx(i, n + i)"],
            "allpairs": ["for i in range(n):", "    for j in range(n):", "        qc.cz(i, n + j)"],
            "shift": ["for i in range(n - 1):", "    qc.cx(i, n + i + 1)"]},
    "K1": {"full": ["for i in range(n):", "    qc.h(i)", "    for j in range(i + 1, n):",
                    "        qc.cp(math.pi / 2 ** (j - i), j, i)"],
           "k2_swaps": ["for i in range(n):", "    qc.h(i)", "    for d in (1, 2):", "        if i + d < n:",
                        "            qc.cp(math.pi / 2 ** d, i + d, i)", "for i in range(n // 2):",
                        "    qc.swap(i, n - 1 - i)"]},
    "K2": {"r1": ["qc.h(range(n))", "for i in range(n):", "    if not (a >> i) & 1:", "        qc.x(i)",
                  "qc.mcp(math.pi, list(range(n - 1)), n - 1)", "for i in range(n):", "    if not (a >> i) & 1:",
                  "        qc.x(i)", "qc.h(range(n))", "qc.x(range(n))", "qc.mcp(math.pi, list(range(n - 1)), n - 1)",
                  "qc.x(range(n))", "qc.h(range(n))"]},
    "K3": {"powers": ["qc.x(n)", "for j in range(n):", "    qc.h(j)", "    qc.cp(math.pi / 4 * 2 ** j, j, n)"],
           "linear": ["qc.x(n)", "for j in range(n):", "    qc.h(j)", "    qc.cp(math.pi / 4 * (j + 1), j, n)"]},
    "K4": {"cdkm": ["for i in range(n):", "    u_maj(qc, 2 * i, 2 * i + 1, 2 * i + 2)",
                    "for i in reversed(range(n)):", "    u_uma(qc, 2 * i, 2 * i + 1, 2 * i + 2)"]},
    "K6": {"tree": ["qc.h(0)", "k = 1", "while k < n:", "    for j in range(0, min(k, n - k)):",
                    "        qc.cx(j, j + k)", "    k *= 2"],
           "w": ["qc.x(0)", "for i in range(n - 1):", "    qc.ry(2 * math.acos(math.sqrt(1 / (n - i))), i + 1)",
                 "    qc.cx(i + 1, i)", "    qc.cx(i, i + 1)"]},
    "K7": {"all_L1": ["qc.ry(math.pi / 2, range(n))", "for i in range(n):", "    for j in range(i + 1, n):",
                      "        qc.cx(i, j)"]},
}
OUTSIDE_QUBITS = {"B12": {"vchain": "2 * n - 2"}, "B17": "2 * n", "K3": "n + 1", "K4": "2 * n + 1"}
OUTSIDE_BITS = {"K2"}
OUTSIDE_BLOCKS = {"B8": ["cx", "cz"], "K4": ["maj", "uma"]}

CLASSES = {
    "B1": (["up1", "down1", "up2", "up3"], lambda v, b: (
        [ladder(b, down=v.startswith("down"), step=int(v[-1]))], False)),
    "B2": (["w3_s1", "w4_s1", "w3_s3", "w4_s4"], lambda v, b: (window(b, int(v[1]), int(v[-1])), False)),
    "B3": (["L1", "L2", "L3"], lambda v, b: (brick(b, int(v[1])), False)),
    "B4": (["first", "last", "two"], lambda v, b: ([{
        "first": Loop(("c", 1), ("n", 0), 1, [Call(b, [("c", 0), ("i", 0)])]),
        "last": Loop(("c", 0), ("n", -1), 1, [Call(b, [("n", -1), ("i", 0)])]),
        "two": Loop(("c", 2), ("n", 0), 1, [Call(b, [("c", 0), ("i", 0)]), Call(b, [("c", 1), ("i", 0)])])}[v]], False)),
    "B5": (["end", "start"], lambda v, b: (ring(b, v), False)),
    "B6": (["ladder_last", "star_last", "hlayer_last", "xlayer_last", "ladder_fixed", "ladder_none"],
           lambda v, b: (mirror(v, b), False)),
    "B7": (["a_i", "a_i1", "xor", "not_a_i"], lambda v, b: (
        [Loop(("c", 0), ("n", -1), 1, bit_if(v, b))], True)),
    "B9": (["K1", "K2", "K5", "K9"], lambda v, b: ([ladder(b, start=int(v[1:]))], False)),
    "B13": (["shift1", "shift2", "alt"], lambda v, b: (swaps(v), False)),
    "B15": (["all_h", "alt_h", "all_x", "all_t"], lambda v, b: ([layer(v.split("_")[1], 2 if v == "alt_h" else 1)], False)),
    "K1": (["k2", "k1", "k3", "const_k2"], lambda v, b: (qft_band(int(v[-1]), (lambda d: "math.pi / 4") if v.startswith("const")
                                                                  else (lambda d: "math.pi / %d" % 2 ** d)), False)),
    "K5": (["x_ai", "z_ai", "anc_ai", "cz_ai", "x_and", "x_not", "anc_and", "anc_not"], lambda v, b: oracle(v)),
    "K6": (["ladder", "star"], lambda v, b: ghz(v)),
    "K7": (["chain_L1", "ring_L1", "chain_L2"], lambda v, b: ansatz(v)),
    "K8": (["chain_T1", "ring_T1", "chain_T2", "xx_T1"], lambda v, b: trotter(v)),
    "K9": (["diag", "anti", "blockdiag", "quarters"], lambda v, b: pattern(v)),
}
for c, vs in OUTSIDE.items():                      # variants outside the class are added to their class
    old = CLASSES.get(c, ([], None))
    CLASSES[c] = (old[0] + ["out_" + v for v in vs], old[1])
CLASSES["K10"] = (["plus", "minus", "basis", "cat"], lambda v, b: prep(v))   # added last: earlier pair numbers do not change
CLASS_BLOCKS = {"B1": BLOCKS2, "B3": BLOCKS2, "B4": BLOCKS2, "B5": BLOCKS2, "B6": BLOCKS2, "B7": BLOCKS2, "B9": BLOCKS2,
                "B2": BLOCKS3}


def blocks_of(cls, v, variants):
    if v.startswith("out_") or cls not in CLASS_BLOCKS or (cls == "B6" and "layer" in v):
        return [None]
    if cls == "B2" and v.startswith("w4"):
        return ["w4"]
    return CLASS_BLOCKS[cls] if v == variants[0] else (FEW3 if cls == "B2" else FEW)


# --- mutation kinds (Table 3): one edit of the tree.
# Each edit returns (new body, new blocks {name: gates}) or None if it does not apply to the task.

def main_loop(body):
    return next((k for k, s in enumerate(body) if isinstance(s, Loop)), None)


def with_loop(prog, k, **kw):
    return prog[:k] + [replace(prog[k], **kw)] + prog[k + 1:]


def first_block(body):
    for s in body:
        if isinstance(s, Call):
            return s.block
        b = first_block(getattr(s, "body", []) + getattr(s, "then", []) + getattr(s, "other", []))
        if b:
            return b


def block_edit(fn):
    def edit(body, blocks):
        b = first_block(body)
        if b is None or b.endswith("_inv"):
            return None
        gs = fn(blocks[b])
        if gs is None:
            return None
        nb = b + "_m"
        return rename(body, b, nb), {nb: gs}
    return edit


def rename(body, b, nb):
    def go(s):
        if isinstance(s, Call):
            return replace(s, block=nb) if s.block == b else s
        if isinstance(s, Loop):
            return replace(s, body=[go(t) for t in s.body])
        if isinstance(s, If):
            return replace(s, then=[go(t) for t in s.then], other=[go(t) for t in s.other])
        return s
    return [go(s) for s in body]


def first_of(gs, ops, new):
    for k, (g, q, a) in enumerate(gs):
        if g in ops:
            return gs[:k] + [new(g, q, a)] + gs[k + 1:]


def map_body(body, f):
    """apply f to every gate and call of the main loop"""
    k = main_loop(body)
    if k is None:
        return None
    def go(s):
        if isinstance(s, If):
            return replace(s, then=[go(t) for t in s.then], other=[go(t) for t in s.other])
        return f(s)
    return with_loop(body, k, body=[go(s) for s in body[k].body])


def changed(body, nb):
    return (nb, {}) if nb is not None and nb != body else None


def swap_operands(s):
    return replace(s, qs=[s.qs[1], s.qs[0]] + s.qs[2:]) if isinstance(s, (Call, Gate)) and len(s.qs) >= 2 else s


def fold(s):
    m = lambda r: ("m", -r[1]) if r[0] == "i" else r
    return replace(s, qs=[m(q) for q in s.qs]) if isinstance(s, (Call, Gate)) else s


def loop_edit(fn):
    def edit(body, blocks):
        k = main_loop(body)
        if k is None:
            return None
        nb = fn(body, k, body[k])
        return None if nb is None else (nb, {})
    return edit


def shift(r, d):
    return (r[0], r[1] + d)


def c1(body, k, L):
    return with_loop(body, k, start=shift(L.start, 1 if L.step > 0 else -1))


def c2(body, k, L):
    return with_loop(body, k, stop=shift(L.stop, -1 if L.step > 0 else 1)) if abs(L.step) == 1 else None


def c7(body, k, L):
    return with_loop(body, k, stop=shift(L.stop, -1)) if L.step == 2 else None


def c3(body, k, L):
    return with_loop(body, k, step=L.step + (1 if L.step > 0 else -1))


def c4(body, k, L):
    if L.step == 1:
        return with_loop(body, k, rev=True)
    if L.step == -1:
        return with_loop(body, k, start=shift(L.stop, 1), stop=shift(L.start, 1), step=1)


def c5(body, k, L):
    loops = [j for j, s in enumerate(body) if isinstance(s, Loop)]
    return body[:loops[-1]] + body[loops[-1] + 1:] if len(loops) >= 2 else None


def c6(K):
    def f(body, k, L):
        if L.step > 0 and L.start[0] == "c" and L.start[1] < K:
            return body[:k + 1] + [replace(L, start=("c", K))] + body[k + 1:]
    return f


def r4(body, k, L):
    """the first iteration moved out of the loop: at n where the loop is empty, the solution fails"""
    if L.step != 1 or L.start != ("c", 0) or L.rev:
        return None
    def at0(s):
        put = lambda r: ("c", r[1]) if r[0] == "i" else r
        if isinstance(s, If):
            return replace(s, bit=put(s.bit), then=[at0(t) for t in s.then], other=[at0(t) for t in s.other])
        return replace(s, qs=[put(q) for q in s.qs])
    return body[:k] + [at0(s) for s in L.body] + [replace(L, start=("c", 1))] + body[k + 1:]


def bit_edit(fn):
    def edit(body, blocks):
        k = main_loop(body)
        if k is None or not any(isinstance(s, If) for s in body[k].body):
            return None
        return with_loop(body, k, body=[fn(s) if isinstance(s, If) else s for s in body[k].body]), {}
    return edit


def cp_loops(body):
    return [k for k, s in enumerate(body) if isinstance(s, Loop) and any(isinstance(t, Gate) and t.op == "cp" for t in s.body)]


def u1(body, blocks):
    """the angle formula pi / X replaced by pi / 2X in all rotations inside loops (cp in K1, ry, rz, rx in K7, K8)"""
    half = lambda t: replace(t, angle="math.pi / %d" % (2 * int(t.angle.split("/ ")[1]))) \
        if isinstance(t, Gate) and t.angle and t.angle.startswith("math.pi / ") else t
    return changed(body, [replace(s, body=[half(t) for t in s.body]) if isinstance(s, Loop) else s for s in body])


def u2(body, blocks):
    """an extra ladder of phases at distance k + 1"""
    ks = cp_loops(body)
    if not ks:
        return None
    k, last = len(ks), body[ks[-1]].body[0]
    angle = last.angle if len(ks) > 1 and body[ks[0]].body[0].angle == last.angle else "math.pi / %d" % 2 ** (k + 1)
    extra = Loop(("c", 0), ("n", -(k + 1)), 1, [Gate("cp", [("i", k + 1), ("i", 0)], angle)])
    return body[:ks[-1] + 1] + [extra] + body[ks[-1] + 1:], {}


def a1(body, blocks):
    """ancilla not uncomputed: the gates after the last loop removed"""
    last = max(k for k, s in enumerate(body) if isinstance(s, Loop))
    return changed(body, body[:last + 1])


def a2(body, blocks):
    """ancilla in |+> instead of |->: the X gates on the ancilla removed"""
    return changed(body, [s for s in body if not (isinstance(s, Gate) and s.op == "x" and s.qs == [("n", 0)])])


def f3(body, blocks):
    """lost amplitude sign: an X or Z on qubit 0 before the first loop removed"""
    k = main_loop(body)
    j = next((j for j, s in enumerate(body[:k]) if isinstance(s, Gate) and s.op in ("x", "z") and s.qs == [("c", 0)]), None)
    return None if j is None else (body[:j] + body[j + 1:], {})


def tail(text):
    return lambda body, blocks: (body + [Raw(text)], {})


def r3(N0, cls, v, b):
    """hard-coded up to N0: for every n <= N0 the explicit gate list of the reference, then an empty circuit"""
    def edit(body, blocks):
        if any(isinstance(s, Raw) for s in body):
            return None
        out = []
        for m in range(1, N0 + 1):
            gs = expand(body, blocks, m, 0)
            out += ["if n == %d:" % m] + ["    qc.%s(%s)" % (g, ", ".join(([a] if a else []) + [str(x) for x in q]))
                                         for g, q, a in gs] + (["    pass"] if not gs else [])
        return [Raw(out)], {}
    return edit


# kind code -> (edit, which tasks, answer of the method expected by construction, where the difference shows)
# which tasks: "block" the default variant only, all blocks; "all" all tasks; "once" one task of the class;
# "anc" tasks with an ancilla
ERRORS = {
    "G1": (block_edit(lambda gs: gs + [("s", (1,), None)]), "block", "not equal", "all n"),
    "G2": (block_edit(lambda gs: gs[:-1] if len(gs) > 1 else None), "block", "not equal", "all n"),
    "G3": (block_edit(lambda gs: first_of(gs, ("cx", "cz", "ch", "swap", "ccx", "cswap"), lambda g, q, a: (
        {"cx": "cz", "ccx": "cswap", "cswap": "ccx"}.get(g, "cx"), q, a))), "block", "not equal", "all n"),
    "G4": (lambda body, blocks: changed(body, map_body(body, swap_operands)), "all", "not equal", "all n"),
    "G5": (block_edit(lambda gs: first_of(gs, ("t",), lambda g, q, a: ("p", q, "math.pi / 8"))), "block", "not equal", "all n"),
    "G6": (block_edit(lambda gs: first_of(gs, ("z", "s", "t"), lambda g, q, a: (
        "rz", q, {"z": "math.pi", "s": "math.pi / 2", "t": "math.pi / 4"}[g]))), "block", "equal", "global phase only"),
    "G7": (lambda body, blocks: ([Raw(["qc.mcx(list(range(0, n - 1)), n - 1)"])], {}), "once", "refuse", "all n"),
    "C1": (loop_edit(c1), "all", "not equal", "all n"),
    "C2": (loop_edit(c2), "all", "not equal", "all n"),
    "C3": (loop_edit(c3), "all", "not equal", "n >= 3"),
    "C4": (loop_edit(c4), "all", "not equal", "n >= 3"),
    "C5": (loop_edit(c5), "all", "not equal", "all n"),
    "C6_K3": (loop_edit(c6(3)), "all", "not equal", "n > 3"),
    "C6_K11": (loop_edit(c6(11)), "all", "not equal", "n > 11"),
    "C6_K20": (loop_edit(c6(20)), "all", "not equal", "n > 20"),
    "C7": (loop_edit(c7), "all", "not equal", "odd or even n"),
    "I3": (lambda body, blocks: changed(body, map_body(body, fold)), "all", "refuse", "n >= 2"),
    "Bb1": (bit_edit(lambda s: replace(s, bit=shift(s.bit, 1))), "all", "not equal", "some a"),
    "Bb2": (bit_edit(lambda s: replace(s, then=s.other, other=s.then)), "all", "not equal", "all a"),
    "Bb3": (bit_edit(lambda s: replace(s, bit=("m", -s.bit[1]))), "all", "refuse", "asymmetric a"),
    "Bb4": (tail(["if a == (1 << n) - 1:", "    qc.z(0)"]), "block", "refuse", "a = 1...1"),
    "R1_N3": (tail(["if n > 3:", "    qc.z(0)"]), "once", "refuse", "n > 3"),
    "R1_N11": (tail(["if n > 11:", "    qc.z(0)"]), "once", "refuse", "n > 11"),
    "R1_N20": (tail(["if n > 20:", "    qc.z(0)"]), "once", "refuse", "n > 20"),
    "R2": (tail(["if n % 3 == 1:", "    qc.z(0)"]), "once", "refuse", "n = 1 mod 3"),
    "R4": (loop_edit(r4), "all", "not equal", "n = 1"),
    "F1": (lambda body, blocks: ([Gate("s", [("c", 0)])] + body, {}), "block", "not equal", "superposition"),
    "U1": (u1, "all", "not equal", "all n"),
    "U2": (u2, "all", "not equal", "n > k + 1"),
    "A1": (a1, "anc", "not equal", "all n"),
    "A2": (a2, "anc", "not equal", "all n"),
    "F2": (lambda body, blocks: ([Raw(["qc.mcp(math.pi, list(range(1, n)), 0)"])] + body, {}), "once", "refuse",
           "superposition"),
    "F3": (f3, "all", "not equal", "all n"),
}
for N0 in (3, 11, 20):
    ERRORS["R3_N%d" % N0] = (None, "once", "refuse", "n > %d" % N0)
ONCE_ONLY = {"G7": {"B1", "B9"}, "Bb4": {"B7"}, "F3": {"K10"}}    # kinds that are natural only for some classes


def expand(body, blocks, n, a, var=None):
    """gate list of the reference at given n, a (for R3)"""
    out = []
    ev = lambda r: {"c": r[1], "n": n + r[1], "i": var + r[1] if var is not None else None,
                    "m": n - 1 - var + r[1] if var is not None else None}[r[0]]
    for s in body:
        if isinstance(s, Gate):
            out.append((s.op, [ev(q) for q in s.qs], s.angle))
        elif isinstance(s, Call):
            qs = [ev(q) for q in s.qs]
            out += [(g, [qs[k] for k in q], ang) for g, q, ang in blocks[s.block]]
        elif isinstance(s, If):
            out += expand(s.then if (a >> ev(s.bit)) & 1 else s.other, blocks, n, a, var)
        elif isinstance(s, Loop):
            r = range(ev(s.start), ev(s.stop), s.step)
            for v in (reversed(r) if s.rev else r):
                out += expand(s.body, blocks, n, a, v)
        else:
            raise ValueError("Raw in expand")
    return out


# --- correct solutions (styles): the same task written differently

def style_shift(body, blocks):
    """loop counter shifted by 1: for j in range(lo + 1, hi + 1) and indices j - 1"""
    def go(s):
        if isinstance(s, Loop):
            return replace(s, start=shift(s.start, 1), stop=shift(s.stop, 1), body=[sub(t) for t in s.body])
        return s
    def sub(s):
        d = lambda r: shift(r, -1) if r[0] == "i" else shift(r, 1) if r[0] == "m" else r
        if isinstance(s, If):
            return replace(s, bit=d(s.bit), then=[sub(t) for t in s.then], other=[sub(t) for t in s.other])
        return replace(s, qs=[d(q) for q in s.qs])
    if not any(isinstance(s, Loop) for s in body):
        return None
    return [go(s) for s in body], {}


def style_dec(body, blocks):
    used = sorted({b for b in blocks if rename(body, b, b + "_dec") != body})
    changed = [b for b in used if decompose(blocks[b]) != blocks[b]]
    if not changed:
        return None
    for b in changed:
        body = rename(body, b, b + "_dec")
    return body, {b + "_dec": decompose(blocks[b]) for b in changed}


def style_vector(body, blocks):
    """a layer of single-qubit gates in one line: qc.h(range(...))"""
    k = main_loop(body)
    if k is None or len(body[k].body) != 1 or not isinstance(body[k].body[0], Gate) or len(body[k].body[0].qs) != 1:
        return None
    L, g = body[k], body[k].body[0]
    if not (L.start[0] == "c" and L.stop[0] == "n" and 1 + L.stop[1] > L.start[1]):
        return None                              # Qiskit refuses an empty range: at n = 1 the call would fail
    r = "range(%s, %s%s)" % (idx(L.start), idx(L.stop), "" if L.step == 1 else ", %d" % L.step)
    return body[:k] + [Raw(["qc.%s(%s)" % (g.op, r)])] + body[k + 1:], {}


def style_commute(body, blocks):
    """single-qubit gates on different qubits commute: traversal in reverse order"""
    k = main_loop(body)
    if k is None or not all(isinstance(s, Gate) and len(s.qs) == 1 for s in body[k].body):
        return None
    return with_loop(body, k, rev=not body[k].rev), {}


STYLES = {"shift": style_shift, "helper": None, "dec": style_dec, "vector": style_vector, "commute": style_commute}
# correct solutions specific to a class: (class, variant) -> [(style name, variant whose body is the solution)]
CLASS_STYLES = {("K5", "anc_ai"): [("noanc", "z_ai")], ("K5", "anc_and"): [("noanc", "z_and")],
                ("K5", "anc_not"): [("noanc", "z_not")], ("K6", "ladder"): [("star", "star")],
                ("K10", "plus"): [("ry", "plus_ry")], ("K10", "basis"): [("ry", "basis_ry")],
                ("K10", "cat"): [("ry", "cat_ry")],
                ("K10", "minus"): [("zfirst", "minus_zfirst"), ("zlast", "minus_zlast"), ("star", "minus_star")]}


# --- assembly

def starts(body):
    """constants at the loop starts: the ground-truth computation must look past them (ground_truth.py)"""
    out = set()
    for s in body:
        if isinstance(s, Loop):
            out |= {s.start[1]} if s.start[0] == "c" else set()
            out |= starts(s.body)
    return out


def loop_edges(body):
    """sizes where a loop with ends c and n + s changes its number of iterations: c - s + 1 + k|t|, k = 0, 1, 2.
    Needed for edits visible only above the loop threshold (C6 with a window of width w shows from K + w, not
    from K + 1; C3 from the second iteration)."""
    out = set()
    for s in body:
        if isinstance(s, Loop):
            ends = {s.start[0], s.stop[0]}
            if ends == {"c", "n"}:
                c = s.start[1] if s.start[0] == "c" else s.stop[1]
                d = s.stop[1] if s.stop[0] == "n" else s.start[1]
                out |= {c - d + 1 + k * abs(s.step) for k in range(3)}
            out |= loop_edges(s.body)
        for b in (getattr(s, "then", []), getattr(s, "other", [])):
            out |= loop_edges(b)
    return out


def task_blocks(body):
    names = set()
    def go(s):
        if isinstance(s, Call):
            names.add(s.block)
        for t in getattr(s, "body", []) + getattr(s, "then", []) + getattr(s, "other", []):
            go(t)
    for s in body:
        go(s)
    out = {}
    for nm in names:
        base = nm[:-4] if nm.endswith("_inv") else nm
        out[nm] = inverse(BLOCKS[base]) if nm.endswith("_inv") else BLOCKS[nm]
    return out


def build_task(cls, v, b, make):
    """body, bits, number of ancillas, task relation, number of qubits of the circuit (a string if not n + ancillas)"""
    if v.startswith("out_"):
        q = OUTSIDE_QUBITS.get(cls)
        return [Raw(OUTSIDE[cls][v[4:]])], cls in OUTSIDE_BITS, 0, "E5phase" if cls == "K6" else "E2", \
            q.get(v[4:]) if isinstance(q, dict) else q
    got = make(v, b)
    body, bits, extra, label = got if len(got) == 4 else got + (0, "E2")
    return body, bits, extra, label, None


def main():
    pairs = []
    for cls, (variants, make) in CLASSES.items():
        src, helpers = [], {b: BLOCKS[b] for b in OUTSIDE_BLOCKS.get(cls, [])}
        for v in variants:
            for b in blocks_of(cls, v, variants):
                body, bits, anc, label, qubits = build_task(cls, v, b, make)
                blocks = task_blocks(body)
                helpers.update(blocks)
                tid = "%s_%s%s" % (cls, v, "_" + b if b else "")
                ref = tid.lower() + "__ref"
                src += function(ref, body, bits, extra=anc, qubits=qubits) + ["", ""]
                default = v == variants[0]
                first = default and b == blocks_of(cls, v, variants)[0]
                todo = []                    # (kind, (body, new blocks), expected answer, where visible, helper, counter, ancillas)
                if v.startswith("out_"):     # outside the class: only the reference against itself, a refusal expected
                    todo.append(("style:self", (body, {}), "refuse", None, False, "i", anc))
                else:
                    for st, fn in STYLES.items():
                        got = (body, {}) if st == "helper" else fn(body, blocks)
                        if got:
                            todo.append(("style:" + st, got, "equal", None, st == "helper", "j" if st == "shift" else "i", anc))
                    for st, w in CLASS_STYLES.get((cls, v), []):
                        other = build_task(cls, w, b, make)
                        todo.append(("style:" + st, (other[0], {}), "equal", None, False, "i", other[2]))
                for code, (fn, where, expect, when) in ERRORS.items():
                    if v.startswith("out_") or where == "block" and not default or where == "once" and not first:
                        continue
                    if where == "anc" and not anc or code in ONCE_ONLY and cls not in ONCE_ONLY[code]:
                        continue
                    if code.startswith("R3"):
                        if bits:                               # no R3 for tasks with bits a: a table over all a
                            continue
                        fn = r3(int(code[4:]), cls, v, b)
                    got = fn(body, blocks)
                    if got and got[1]:                         # an own name of the edited block for every kind
                        (old, gs), = got[1].items()
                        new = old[:-2] + "_" + code.lower()
                        got = rename(got[0], old, new), {new: gs}
                    if got:
                        todo.append((code, got, expect, when, False, "i", anc))
                for kind, (nb, extra), expect, when, helper, var, e in todo:
                    edge = {c + d for c in starts(body) | starts(nb) for d in (2, 3) if c < 10}   # loop start
                    edge |= {int(x[1:]) + d for x in kind.split("_")[1:] for d in (1, 2)}  # threshold K or N0 of the kind
                    edge |= {e + d for e in loop_edges(body) | loop_edges(nb) for d in (0, 1) if e + d <= 24}
                    check_n = sorted(c for c in edge if c > 7)
                    name = "%s__%s" % (tid.lower(), kind.replace("style:", "st_"))
                    helpers.update(extra)
                    src += function(name, nb, bits, var=var, helper=helper, extra=e, qubits=qubits) + ["", ""]
                    pairs.append({
                        "id": "%s/%s" % (tid, kind), "task": tid, "cls": cls, "variant": v, "block": b, "label": label,
                        "ref": {"file": "%s/gen_%s.py" % (REL, cls.lower()), "function": ref},
                        "sol": {"file": "%s/gen_%s.py" % (REL, cls.lower()), "function": name},
                        "origin": "benchmark", "kind": kind, "check_n": check_n,
                        "expect": {label: True if kind.startswith("style") else expect == "equal",
                                   "method": expect, "when": when}})
        head = ['"""Benchmark, class %s. Generated by experiments/benchmark/make_benchmark.py; do not edit by hand."""' % cls,
                "import math", "", "from qiskit import QuantumCircuit", "", ""]
        text = head + sum((helper_text(k, helpers[k]) for k in sorted(helpers)), []) + src
        with open(os.path.join(HERE, "gen_%s.py" % cls.lower()), "w") as f:
            f.write("\n".join(text).rstrip() + "\n")
    with open(os.path.join(HERE, "pairs.jsonl"), "w") as f:
        for p in pairs:
            f.write(json.dumps(p, ensure_ascii=False) + "\n")
    print(len(pairs), "pairs")


if __name__ == "__main__":
    main()
