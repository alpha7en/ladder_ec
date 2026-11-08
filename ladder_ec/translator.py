"""Translator: Python with Qiskit -> ladder program (Sec. 6.1 of the paper, Online Resource 1 Sec. S2).

Input: a source file and the name of a function `solve(n, a) -> QuantumCircuit` (or `solve(n)`).
Output: a JSON description of the ladder program (schema "ladder/v1") with the threshold n_tr(P) (field
"n_min"), or a refusal with a reason code (Table S2).

The generator code is NOT executed: only its AST is read. The translator accepts the constructs listed below and
refuses everything else. For an accepted function the program repeats the Python semantics for all
n >= n_tr(P) and all bit strings a (Thm. 3); the proof obligations of Table S1 are collected as affine forms
that must be nonnegative and turned into thresholds.

Accepted fragment (Sec. S2.1):
- `qc = QuantumCircuit(n + e)`, exactly one circuit, and `return qc` at the end;
- names bound to affine expressions in n (and the loop counter), `list(range(lo, hi))` and slices of it;
  rebinding only in the same context (Sec. S2.1);
- fixed-arity gates (GATES below) with indices c, n + c, i + c and constant angles q*pi;
- single-qubit gates on `range(...)` (expanded into a loop);
- `for i in range(...)` / `reversed(range(...))`: either bounded (both bounds of the same kind, expanded in
  full) or a ladder from a constant to n + c with step +-s; inside a ladder only bounded loops;
- `if <bit>:` / `else:`, where the bit is `a[idx]`, `(a >> idx) & 1` or `a & (1 << idx)` (Sec. S2.3);
- calls of helper functions of the module (and of modules next to it, `from m import f`) are inlined;
  recursion is refused.
Everything else, including any branching on n, `%`, `//`, the conditional expression, `while`, nested ladders,
multi-controlled gates of growing arity and compose/append/inverse, is refused.

In the reason codes and in the JSON, "fixed" means a qubit or a bit a[c] with a constant position (k or n + k)
used inside a loop body: the sets S and E of Definition 3.
"""
from __future__ import annotations

import ast
import json
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

SCHEMA = "ladder/v1"

# QuantumCircuit method -> (number of qubits, number of angles)
GATES = {
    "h": (1, 0), "x": (1, 0), "y": (1, 0), "z": (1, 0), "s": (1, 0), "sdg": (1, 0),
    "t": (1, 0), "tdg": (1, 0), "rz": (1, 1), "p": (1, 1), "rx": (1, 1), "ry": (1, 1),
    "cx": (2, 0), "cz": (2, 0), "cy": (2, 0), "ch": (2, 0), "swap": (2, 0),
    "cp": (2, 1), "crz": (2, 1), "rzz": (2, 1),
    "ccx": (3, 0), "ccz": (3, 0), "cswap": (3, 0),
}
CORE_GATES = {"h", "x", "z", "s", "t", "cx", "cz", "swap", "rz", "p"}
GROWING = {"mcx", "mcp", "mcrx", "mcry", "mcrz", "mct", "mcmt", "mcu1"}
LIBRARY = {"append", "compose", "inverse", "repeat", "power", "control", "unitary",
           "initialize", "prepare_state", "decompose", "tensor", "to_gate", "to_instruction",
           "unitary", "isometry", "diagonal", "ucrz", "ucry", "ucrx", "uc"}
IGNORED = {"barrier"}
BUILTINS_USED = {"range", "list", "reversed", "len"}
MAX_UNROLL = 64
ANGLE_FUNCS = {"pi", "acos", "asin", "atan", "atan2", "sqrt", "cos", "sin", "tan", "exp", "log", "arccos",
               "arcsin", "arctan"}


class Refuse(Exception):
    def __init__(self, code: str, msg: str, node=None):
        super().__init__(f"{code}: {msg}")
        self.code, self.msg = code, msg
        self.line = getattr(node, "lineno", None)
        self.file = None
        self.via: list[str] = []


# ---------------------------------------------------------------- values

@dataclass(frozen=True)
class Lin:
    """c + kn*n + ki*i, where n is the size and i the symbolic ladder counter."""
    c: int
    kn: int = 0
    ki: int = 0

    def __add__(self, o): return Lin(self.c + o.c, self.kn + o.kn, self.ki + o.ki)
    def __neg__(self): return Lin(-self.c, -self.kn, -self.ki)
    def __sub__(self, o): return self + (-o)
    def scale(self, k): return Lin(self.c * k, self.kn * k, self.ki * k)
    @property
    def const(self): return self.kn == 0 and self.ki == 0


@dataclass(frozen=True)
class Arr:          # list(range(start, start + length))
    start: Lin
    length: Lin


@dataclass(frozen=True)
class Rng:          # range(start, stop, step), constant step
    start: Lin
    stop: Lin
    step: int


@dataclass(frozen=True)
class ListVal:
    items: tuple


@dataclass(frozen=True)
class QC:
    pass


@dataclass(frozen=True)
class BitParam:
    name: str


@dataclass(frozen=True)
class BitShift:     # (a >> idx), waiting for "& 1"
    param: str
    idx: Lin


@dataclass(frozen=True)
class Bit:          # condition: bit param[idx] == value
    param: str
    idx: Lin
    access: str     # "index" (a[idx]) or "shift" ((a >> idx) & 1)
    value: int

    def neg(self): return Bit(self.param, self.idx, self.access, 1 - self.value)


@dataclass(frozen=True)
class Poison:       # a loop counter after the loop has ended
    why: str


@dataclass(frozen=True)
class Thunk:        # argument of a helper function: evaluated in the caller's scope
    node: ast.AST
    scope: "Scope"


@dataclass(frozen=True)
class FloatConst:   # float constant of the module (usable only in angles)
    value: float


@dataclass(frozen=True)
class NoneVal:
    pass


@dataclass(frozen=True)
class QList:        # list of qubits of unknown content and of length `length`
    length: Lin


@dataclass(frozen=True)
class AngleVal:     # a name bound to a constant angle
    angle: "Angle"


@dataclass(frozen=True)
class Angle:        # pi_coef * pi + rat
    pi: Fraction
    rat: Fraction

    def to_json(self):
        d = {"pi": str(self.pi), "rational": str(self.rat),
             "value": float(self.pi) * 3.141592653589793 + float(self.rat)}
        d["exact_pi_multiple"] = self.rat == 0
        return d


# ---------------------------------------------------------------- modules

class Module:
    _cache: dict[Path, "Module"] = {}

    def __init__(self, path: Path):
        self.path = path
        self.tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        self.funcs: dict[str, ast.FunctionDef] = {}
        self.consts: dict[str, int] = {}
        self.fconsts: dict[str, float] = {}
        self.imports: dict[str, tuple] = {}
        self.shadowed: set[str] = set()
        assigned: dict[str, list] = {}
        for node in self.tree.body:
            if isinstance(node, ast.FunctionDef):
                self.funcs[node.name] = node
            elif isinstance(node, ast.Import):
                for a in node.names:
                    self.imports[a.asname or a.name] = ("module", a.name)
            elif isinstance(node, ast.ImportFrom):
                for a in node.names:
                    self.imports[a.asname or a.name] = ("from", node.module, a.name, node.level)
            elif isinstance(node, (ast.Assign, ast.AnnAssign, ast.AugAssign)):
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                for t in targets:
                    for nm in _target_names(t):
                        assigned.setdefault(nm, []).append(node)
            elif isinstance(node, ast.ClassDef):
                assigned.setdefault(node.name, []).append(node)
        globals_decl = {nm for n in ast.walk(self.tree) if isinstance(n, ast.Global) for nm in n.names}
        for nm, nodes in assigned.items():
            if nm in BUILTINS_USED or nm == "QuantumCircuit":
                self.shadowed.add(nm)
            if len(nodes) == 1 and isinstance(nodes[0], ast.Assign) and nm not in globals_decl \
                    and len(nodes[0].targets) == 1:
                t, val = nodes[0].targets[0], nodes[0].value
                if isinstance(t, ast.Tuple) and isinstance(val, ast.Tuple) and len(t.elts) == len(val.elts):
                    val = val.elts[[e.id if isinstance(e, ast.Name) else None for e in t.elts].index(nm)]
                elif not isinstance(t, ast.Name):
                    continue
                try:
                    self.consts[nm] = _const_int(val, self.consts)
                except Refuse:
                    if isinstance(val, ast.Constant) and type(val.value) is float:
                        self.fconsts[nm] = val.value
        for nm in list(self.funcs):
            if nm in BUILTINS_USED or nm == "QuantumCircuit":
                self.shadowed.add(nm)
        for nm, imp in self.imports.items():
            if nm in BUILTINS_USED:
                self.shadowed.add(nm)

    @classmethod
    def load(cls, path: Path) -> "Module":
        path = path.resolve()
        if path not in cls._cache:
            cls._cache[path] = Module(path)
        return cls._cache[path]

    def resolve_func(self, name):
        if name in self.funcs:
            return self, self.funcs[name]
        imp = self.imports.get(name)
        if imp and imp[0] == "from" and imp[1] and imp[3] == 0:
            p = self.path.parent / (imp[1].replace(".", "/") + ".py")
            if p.exists():
                m = Module.load(p)
                return m.resolve_func(imp[2])
        return None

    def is_qiskit_qc(self, name):
        imp = self.imports.get(name)
        return (imp is not None and imp[0] == "from" and imp[1] in ("qiskit", "qiskit.circuit")
                and imp[2] == "QuantumCircuit" and "QuantumCircuit" not in self.shadowed)

    def is_math_module(self, name):
        imp = self.imports.get(name)
        return imp is not None and imp[0] == "module" and imp[1] in ("math", "numpy")


def _target_names(t):
    if isinstance(t, ast.Name):
        return [t.id]
    if isinstance(t, (ast.Tuple, ast.List)):
        return [x for e in t.elts for x in _target_names(e)]
    return []


def _const_int(node, consts):
    """integer constant of the module: literals, + - * // % **, references to other constants"""
    if isinstance(node, ast.Constant) and type(node.value) is int:
        return node.value
    if isinstance(node, ast.Name) and node.id in consts:
        return consts[node.id]
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -_const_int(node.operand, consts)
    if isinstance(node, ast.BinOp):
        a, b = _const_int(node.left, consts), _const_int(node.right, consts)
        ops = {ast.Add: lambda: a + b, ast.Sub: lambda: a - b, ast.Mult: lambda: a * b,
               ast.FloorDiv: lambda: a // b, ast.Mod: lambda: a % b, ast.Pow: lambda: a ** b}
        f = ops.get(type(node.op))
        if f and not (isinstance(node.op, ast.Pow) and (b < 0 or b > 64)):
            return f()
    raise Refuse("not_const", "not an integer constant", node)


# ---------------------------------------------------------------- scopes

LOOP_VAR_CTX = -1


class Scope:
    """Names of a function. Rebinding is allowed only in the context where the name was bound (the same function
    level or the same loop iteration): values are computed at once, so straight-line rebinding is exact. An
    accumulator across iterations is refused (Lemma S.16)."""

    def __init__(self, module: Module, params: dict, size_param: str | None, ctx: int = 0,
                 base_conds: tuple = (), bounds: list | None = None):
        self.module = module
        self.bounds = [] if bounds is None else bounds     # (Lin, node): proof obligations, forms that must be >= 0
        self.base_conds = base_conds
        self.vars: dict[str, object] = dict(params)
        self.ctx_of: dict[str, int] = {k: ctx for k in params}
        self.size_param = size_param

    def lookup(self, name, node):
        if name in self.vars:
            return self.vars[name]
        if name in self.module.consts:
            return Lin(self.module.consts[name])
        if name in self.module.fconsts:
            return FloatConst(self.module.fconsts[name])
        raise Refuse("unknown_name", f"name `{name}` is not in the accepted fragment", node)

    def bind(self, name, value, node, ctx):
        old = self.vars.get(name)
        if isinstance(old, QC):
            raise Refuse("alias_reassigned", f"the circuit `{name}` is rebound", node)
        if name in self.vars and not isinstance(old, Poison) and self.ctx_of.get(name) != ctx:
            what = "loop counter" if self.ctx_of.get(name) == LOOP_VAR_CTX else "name"
            raise Refuse("alias_reassigned",
                         f"{what} `{name}` is rebound inside a loop (accumulator across iterations)", node)
        self.vars[name] = value
        self.ctx_of[name] = ctx


# ---------------------------------------------------------------- expressions

def _deref(v):
    while isinstance(v, Thunk):
        v = eval_value(v.node, v.scope)
    return v


def eval_value(node, scope: Scope):
    """value of an expression in an index context: Lin, Arr, Rng, ListVal, QC, BitParam, BitShift (Table S1)"""
    if isinstance(node, ast.Constant):
        if type(node.value) is int:
            return Lin(node.value)
        if node.value is None:
            return NoneVal()
        raise Refuse("nonint_constant", f"constant {node.value!r} where an index is expected", node)
    if isinstance(node, ast.Name):
        v = scope.lookup(node.id, node)
        if isinstance(v, Thunk):
            v = _deref(v)
        if isinstance(v, Poison):
            raise Refuse("loop_var_after_loop", v.why, node)
        if isinstance(v, (FloatConst, AngleVal)):
            raise Refuse("nonint_constant", f"`{node.id}` is a float where an index is expected", node)
        return v
    if isinstance(node, ast.IfExp):
        raise Refuse("conditional_expression", "conditional expression `x if c else y`", node)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        v = eval_value(node.operand, scope)
        if isinstance(v, Lin):
            return -v
        raise Refuse("nonaffine_expression", "unary minus on a non-index value", node)
    if isinstance(node, ast.BinOp):
        a, b = eval_value(node.left, scope), eval_value(node.right, scope)
        if isinstance(node.op, ast.RShift) and isinstance(a, BitParam) and isinstance(b, Lin):
            return BitShift(a.name, b)
        for v in (a, b):
            if isinstance(v, (BitParam, BitShift)):
                raise Refuse("bit_param_misuse", "the bit-string parameter is used other than as a bit", node)
        lists = (ListVal, Arr, QList)
        if isinstance(node.op, ast.Add) and isinstance(a, lists) and isinstance(b, lists):
            if isinstance(a, ListVal) and isinstance(b, ListVal):
                return ListVal(a.items + b.items)
            return QList(_len(a) + _len(b))
        if not (isinstance(a, Lin) and isinstance(b, Lin)):
            raise Refuse("nonaffine_expression", "arithmetic on values that are not integer indices", node)
        if isinstance(node.op, ast.Add):
            return a + b
        if isinstance(node.op, ast.Sub):
            return a - b
        if isinstance(node.op, ast.Mult):
            if a.const:
                return b.scale(a.c)
            if b.const:
                return a.scale(b.c)
            raise Refuse("nonaffine_expression", "product of two variables", node)
        if a.const and b.const:
            try:
                return Lin(_const_int(ast.BinOp(ast.Constant(a.c), node.op, ast.Constant(b.c)), {}))
            except Refuse:
                pass
        opname = type(node.op).__name__
        raise Refuse("nonaffine_expression", f"operation {opname} on n or on the loop counter", node)
    if isinstance(node, ast.Subscript):
        base = eval_value(node.value, scope)
        if isinstance(node.slice, ast.Slice):
            sl = node.slice
            if isinstance(base, Arr) and sl.step is None:
                lo = eval_value(sl.lower, scope) if sl.lower is not None else Lin(0)
                hi = eval_value(sl.upper, scope) if sl.upper is not None else base.length
                if isinstance(lo, Lin) and isinstance(hi, Lin) and not (lo.const and lo.c < 0) \
                        and not (hi.const and hi.c < 0):
                    scope.bounds += [(lo, node), (hi - lo, node), (base.length - hi, node)]   # no clipping
                    return Arr(base.start + lo, hi - lo)
            raise Refuse("unsupported_expression", "list slice outside the accepted fragment", node)
        idx = eval_value(node.slice, scope)
        if isinstance(base, BitParam):
            raise Refuse("bit_param_misuse", "bit of the parameter outside an `if` condition", node)
        if not isinstance(idx, Lin):
            raise Refuse("nonaffine_expression", "list index is not affine", node)
        if isinstance(base, Arr):
            if idx.const and idx.c < 0:
                idx = base.length + idx
            scope.bounds += [(idx, node), (base.length - idx - Lin(1), node)]   # index inside the list
            return base.start + idx
        if isinstance(base, ListVal) and idx.const:
            return base.items[idx.c]
        raise Refuse("unsupported_expression", "indexing of something that is not a list of qubits", node)
    if isinstance(node, (ast.List, ast.Tuple)):
        return ListVal(tuple(eval_value(e, scope) for e in node.elts))
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
        fn = node.func.id
        if node.keywords:
            raise Refuse("unsupported_expression", f"keyword arguments in `{fn}`", node)
        if fn in BUILTINS_USED and fn not in scope.module.shadowed and fn not in scope.vars:
            args = [eval_value(a, scope) for a in node.args]
            if fn == "range":
                for a in args:
                    if isinstance(a, (BitParam, BitShift)):
                        raise Refuse("bit_param_misuse",
                                     f"the parameter is used as a repetition count: `{ast.unparse(node)}`", node)
                if not all(isinstance(a, Lin) for a in args) or not 1 <= len(args) <= 3:
                    raise Refuse("range_not_ladder", "range arguments are not affine indices", node)
                if len(args) == 1:
                    args = [Lin(0)] + args
                step = args[2] if len(args) == 3 else Lin(1)
                if not step.const or step.c == 0:
                    raise Refuse("range_not_ladder", "range step is not a non-zero constant", node)
                return Rng(args[0], args[1], step.c)
            if fn == "list" and len(args) == 1:
                a = args[0]
                if isinstance(a, Rng) and a.step == 1:
                    scope.bounds.append((a.stop - a.start, node))
                    return Arr(a.start, a.stop - a.start)
                if isinstance(a, (Arr, ListVal, QList)):
                    return a
            if fn == "reversed" and len(args) == 1:
                a = args[0]
                if isinstance(a, Arr):
                    a = Rng(a.start, a.start + a.length, 1)
                if isinstance(a, Rng) and abs(a.step) == 1:
                    s = a.step
                    return Rng(a.stop - Lin(s), a.start - Lin(s), -s)
                raise Refuse("range_not_ladder", "reversed of a range with step other than 1", node)
            if fn == "len" and len(args) == 1:
                a = args[0]
                if isinstance(a, Arr):
                    return a.length
                if isinstance(a, Rng) and a.step == 1:
                    scope.bounds.append((a.stop - a.start, node))
                    return a.stop - a.start
                if isinstance(a, (ListVal, QList)):
                    return _len(a)
            raise Refuse("unknown_call", f"`{fn}` with unsupported arguments", node)
        if fn == "QuantumCircuit":
            raise Refuse("unsupported_statement", "QuantumCircuit outside an assignment", node)
        raise Refuse("unknown_call", f"call of `{fn}` in an expression", node)
    if isinstance(node, ast.Call):
        raise Refuse("unknown_call", f"call of `{ast.unparse(node.func)}` outside the accepted fragment", node)
    if isinstance(node, ast.Attribute):
        raise Refuse("nonint_constant", f"attribute `{ast.unparse(node)}` where an index is expected", node)
    raise Refuse("unsupported_expression", f"expression `{ast.unparse(node)}`", node)


def _len(v):
    if isinstance(v, ListVal):
        return Lin(len(v.items))
    return v.length


def _kinds(node, scope):
    """which kinds of variables an expression uses: size (n), loop (counter), param (bits), growing (a list of
    qubits whose length depends on the counter)"""
    kinds = set()
    for sub in ast.walk(node):
        if isinstance(sub, ast.Name) and isinstance(sub.ctx, ast.Load):
            try:
                v = scope.lookup(sub.id, sub)
                v = _deref(v) if isinstance(v, Thunk) else v
            except Refuse:
                continue
            if isinstance(v, Lin):
                kinds |= {"size"} if v.kn else set()
                kinds |= {"loop"} if v.ki else set()
            elif isinstance(v, (Arr, QList, ListVal, Rng)):
                ln = _len(v) if not isinstance(v, Rng) else v.stop - v.start
                kinds |= {"growing"} if ln.ki else set()
                kinds |= {"size", "size_list"} if ln.kn else set()
            elif isinstance(v, (BitParam, BitShift)):
                kinds.add("param")
    return kinds


def eval_angle(node, scope: Scope) -> Angle:
    if isinstance(node, ast.Constant) and type(node.value) in (int, float):
        return Angle(Fraction(0), Fraction(repr(node.value)))
    if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.attr == "pi" \
            and scope.module.is_math_module(node.value.id) and node.value.id not in scope.vars:
        return Angle(Fraction(1), Fraction(0))
    if isinstance(node, ast.Name):
        v = scope.lookup(node.id, node)
        if isinstance(v, Thunk):
            return eval_angle(v.node, v.scope)
        if isinstance(v, Lin) and v.const:
            return Angle(Fraction(0), Fraction(v.c))
        if isinstance(v, AngleVal):
            return v.angle
        if isinstance(v, FloatConst):
            return Angle(Fraction(0), Fraction(repr(v.value)))
        k = _kinds(node, scope)
        dep = "the loop counter" if "loop" in k else ("n" if "size" in k else ("the parameter" if "param" in k else ""))
        dep = "" if node.id == scope.size_param else dep
        raise Refuse("nonconst_angle", f"angle depends on {dep} `{node.id}`".replace("  ", " "), node)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        a = eval_angle(node.operand, scope)
        return Angle(-a.pi, -a.rat) if isinstance(node.op, ast.USub) else a
    if isinstance(node, ast.BinOp):
        a, b = eval_angle(node.left, scope), eval_angle(node.right, scope)
        if isinstance(node.op, ast.Add):
            return Angle(a.pi + b.pi, a.rat + b.rat)
        if isinstance(node.op, ast.Sub):
            return Angle(a.pi - b.pi, a.rat - b.rat)
        if isinstance(node.op, ast.Mult):
            if a.pi == 0:
                return Angle(b.pi * a.rat, b.rat * a.rat)
            if b.pi == 0:
                return Angle(a.pi * b.rat, a.rat * b.rat)
        if isinstance(node.op, ast.Div) and b.pi == 0 and b.rat != 0:
            return Angle(a.pi / b.rat, a.rat / b.rat)
        if isinstance(node.op, ast.Pow) and a.pi == 0 and b.pi == 0 and b.rat.denominator == 1 \
                and abs(b.rat) <= 64:
            return Angle(Fraction(0), a.rat ** int(b.rat))
        raise Refuse("nonconst_angle", f"angle `{ast.unparse(node)}` does not reduce to q*pi + r", node)
    k = _kinds(node, scope)
    dep = ", ".join(x for x, t in (("n", "size"), ("the loop counter", "loop"), ("the parameter", "param")) if t in k)
    raise Refuse("nonconst_angle", f"angle `{ast.unparse(node)}` " +
                 (f"depends on {dep}" if dep else "is not a constant"), node)


def eval_cond(node, scope: Scope) -> Bit | None:
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
        b = eval_cond(node.operand, scope)
        return b.neg() if b else None
    if isinstance(node, ast.Compare) and len(node.ops) == 1 and isinstance(node.ops[0], (ast.Eq, ast.NotEq)):
        l, r = node.left, node.comparators[0]
        if isinstance(l, ast.Constant):
            l, r = r, l
        if isinstance(r, ast.Constant) and type(r.value) is int and r.value in (0, 1):
            if r.value == 1 and isinstance(l, ast.BinOp) and any(isinstance(x, ast.BinOp) and
                                                                 isinstance(x.op, ast.LShift) for x in (l.left, l.right)):
                return None                             # a & (1 << k) is 0 or 2^k, not 0 or 1
            b = _bit_expr(l, scope)
            if b is None:
                return None
            want = r.value if isinstance(node.ops[0], ast.Eq) else 1 - r.value
            return b if want == 1 else b.neg()
        return None
    return _bit_expr(node, scope)


def _bit_expr(node, scope):
    try:
        if isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name):
            base = eval_value(node.value, scope)
            if isinstance(base, BitParam):
                idx = eval_value(node.slice, scope)
                if isinstance(idx, Lin):
                    return Bit(base.name, idx, "index", 1)
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitAnd):
            for x, y in ((node.left, node.right), (node.right, node.left)):
                if isinstance(y, ast.Constant) and y.value == 1:
                    v = eval_value(x, scope)
                    if isinstance(v, BitShift):
                        return Bit(v.param, v.idx, "shift", 1)
                if isinstance(y, ast.BinOp) and isinstance(y.op, ast.LShift) \
                        and isinstance(y.left, ast.Constant) and y.left.value == 1:
                    v = eval_value(x, scope)
                    idx = eval_value(y.right, scope)
                    if isinstance(v, BitParam) and isinstance(idx, Lin):
                        return Bit(v.name, idx, "shift", 1)
    except Refuse:
        return None
    return None


def classify_branch(node, scope: Scope):
    """why an `if` condition is not a bit of the parameter: what it branches on"""
    kinds = _kinds(node, scope)
    text = ast.unparse(node)
    if "growing" in kinds:
        return "growing_arity", f"condition on a list of qubits whose length grows with the counter: `if {text}`"
    if "size_list" in kinds:
        return "growing_arity", f"condition on a list of qubits whose length depends on n: `if {text}`"
    if "size" in kinds:
        return "branch_on_size", f"branch on n: `if {text}`"
    if "loop" in kinds:
        return "branch_on_loop_index", f"branch on the loop counter: `if {text}`"
    if "param" in kinds:
        return "branch_on_param", f"branch on the parameter other than on one bit: `if {text}`"
    return "branch_unsupported", f"condition outside the accepted fragment: `if {text}`"


def static_truth(node, scope: Scope):
    """True/False when the condition can be evaluated during translation (all values constant), else None"""
    try:
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
            t = static_truth(node.operand, scope)
            return None if t is None else not t
        if isinstance(node, ast.BoolOp):
            vals = [static_truth(v, scope) for v in node.values]
            if isinstance(node.op, ast.And):
                return False if False in vals else (None if None in vals else True)
            return True if True in vals else (None if None in vals else False)
        if isinstance(node, ast.Compare) and len(node.ops) == 1:
            op, r = node.ops[0], node.comparators[0]
            if isinstance(op, (ast.Is, ast.IsNot)) and isinstance(r, ast.Constant) and r.value is None:
                v = eval_value(node.left, scope)
                res = isinstance(v, NoneVal)
                return res if isinstance(op, ast.Is) else not res
            a, b = eval_value(node.left, scope), eval_value(r, scope)
            if isinstance(a, Lin) and isinstance(b, Lin) and a.const and b.const:
                f = {ast.Eq: a.c == b.c, ast.NotEq: a.c != b.c, ast.Lt: a.c < b.c, ast.LtE: a.c <= b.c,
                     ast.Gt: a.c > b.c, ast.GtE: a.c >= b.c}
                return f.get(type(op))
            return None
        v = eval_value(node, scope)
        if isinstance(v, NoneVal):
            return False
        if isinstance(v, Lin) and v.const:
            return v.c != 0
        if isinstance(v, (ListVal, Arr, QList)):
            ln = _len(v)
            return ln.c != 0 if ln.const else None
    except Refuse:
        return None
    return None


# ---------------------------------------------------------------- gate collection

@dataclass
class Gate:
    op: str
    qubits: tuple       # Lin
    angle: Angle | None
    cond: tuple         # Bit
    line: int


class Collector:
    def __init__(self):
        self.segments: list[dict] = []   # {"kind": "fixed", "gates": [...]} | {"kind": "loop", ...}
        self.body: list[Gate] | None = None  # body of the current ladder

    def emit(self, g: Gate):
        if self.body is not None:
            self.body.append(g)
            return
        if not self.segments or self.segments[-1]["kind"] != "fixed":
            self.segments.append({"kind": "fixed", "gates": []})
        self.segments[-1]["gates"].append(g)


class Analyzer:
    def __init__(self, module: Module, fn: ast.FunctionDef):
        self.module, self.fn = module, fn
        self.col = Collector()
        self.nq: Lin | None = None
        self.qc_name = None
        self.stack = []
        self.bit_params: list[str] = []
        self.ctx = [0]
        self.next_ctx = 1
        self.bounds = []          # proof obligations outside ladders; inside ladders with the counter range
        self.loop_bounds = []

    # ---- entry
    def run(self):
        fn = self.fn
        a = fn.args
        if fn.decorator_list:
            raise Refuse("unsupported_statement", "decorated function", fn)
        if a.vararg or a.kwarg or a.kwonlyargs or a.posonlyargs or a.defaults:
            raise Refuse("extra_params", "default parameters, *args or **kwargs", fn)
        names = [x.arg for x in a.args]
        if not names:
            raise Refuse("no_size_param", "the function has no parameters", fn)
        size, bits = names[0], names[1:]
        self.size_param, self.bit_params = size, bits
        params = {size: Lin(0, 1, 0)}
        params.update({b: BitParam(b) for b in bits})
        scope = Scope(self.module, params, size, bounds=self.bounds)
        body = list(fn.body)
        if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                and isinstance(body[0].value.value, str):
            body = body[1:]
        if not body or not isinstance(body[-1], ast.Return):
            raise Refuse("early_return", "the function does not end with `return qc`", fn)
        ret = body[-1]
        self.block(body[:-1], scope, conds=())
        if not (isinstance(ret.value, ast.Name) and isinstance(scope.vars.get(ret.value.id), QC)):
            rv = ret.value
            if isinstance(rv, ast.Call) and isinstance(rv.func, ast.Attribute) and rv.func.attr in LIBRARY:
                raise Refuse("library_gate", f"returns `{ast.unparse(rv)}`", ret)
            raise Refuse("unsupported_statement", "returns something other than the built circuit", ret)
        if self.nq is None:
            raise Refuse("unsupported_statement", "no circuit is created", fn)
        return self.finish()

    # ---- statements
    def block(self, stmts, scope, conds):
        for st in stmts:
            self.stmt(st, scope, conds)

    def stmt(self, st, scope: Scope, conds):
        if isinstance(st, ast.Pass):
            return
        if isinstance(st, ast.Expr):
            if isinstance(st.value, ast.Constant) and isinstance(st.value.value, str):
                return
            if isinstance(st.value, ast.Call):
                return self.call(st.value, scope, conds)
            raise Refuse("unsupported_statement", f"expression statement `{ast.unparse(st)}`", st)
        if isinstance(st, ast.Assign):
            if len(st.targets) != 1:
                raise Refuse("unsupported_statement", "chained assignment", st)
            t, v = st.targets[0], st.value
            if isinstance(t, ast.Tuple) and isinstance(v, ast.Tuple) and len(t.elts) == len(v.elts) \
                    and all(isinstance(e, ast.Name) for e in t.elts):
                pairs = [(e.id, x) for e, x in zip(t.elts, v.elts)]
            elif isinstance(t, ast.Name):
                pairs = [(t.id, v)]
            else:
                raise Refuse("unsupported_statement", f"assignment to `{ast.unparse(t)}`", st)
            if len(pairs) == 1 and isinstance(v, ast.Call) and isinstance(v.func, ast.Name) \
                    and v.func.id == "QuantumCircuit":
                return self.make_qc(pairs[0][0], v, scope, conds, st)
            if conds != scope.base_conds:
                raise Refuse("alias_under_condition", f"assignment to `{ast.unparse(t)}` under a bit condition", st)
            vals = [self.assigned_value(x, scope) for _, x in pairs]   # as in Python: the right-hand side first
            for (name, _), val in zip(pairs, vals):
                scope.bind(name, val, st, self.ctx[-1])
            return
        if isinstance(st, ast.AugAssign):
            raise Refuse("alias_reassigned", f"`{ast.unparse(st)}`", st)
        if isinstance(st, ast.For):
            return self.for_loop(st, scope, conds)
        if isinstance(st, ast.While):
            raise Refuse("while_loop", "while loop", st)
        if isinstance(st, ast.If):
            t = static_truth(st.test, scope)
            if t is not None:
                return self.block(st.body if t else st.orelse, scope, conds)
            b = eval_cond(st.test, scope)
            if b is None:
                code, msg = classify_branch(st.test, scope)
                raise Refuse(code, msg, st)
            self.block(st.body, scope, conds + (b,))
            self.block(st.orelse, scope, conds + (b.neg(),))
            return
        if isinstance(st, ast.Return):
            raise Refuse("early_return", "return before the end of the function", st)
        if isinstance(st, ast.FunctionDef):
            raise Refuse("unsupported_statement", "nested function definition", st)
        raise Refuse("unsupported_statement", f"statement {type(st).__name__}", st)

    def assigned_value(self, node, scope):
        try:
            return eval_value(node, scope)
        except Refuse as rv:
            looks_angle = any(isinstance(x, ast.Constant) and type(x.value) is float or
                              isinstance(x, ast.Attribute) and isinstance(x.value, ast.Name)
                              and scope.module.is_math_module(x.value.id) and x.attr in ANGLE_FUNCS
                              for x in ast.walk(node))
            if not looks_angle:
                raise
            return AngleVal(eval_angle(node, scope))

    def make_qc(self, name, call, scope, conds, st):
        if self.qc_name is not None or conds or self.col.body is not None or self.stack:
            raise Refuse("unsupported_statement", "a second circuit, or a circuit created inside a loop or helper", st)
        if not self.module.is_qiskit_qc("QuantumCircuit") or "QuantumCircuit" in scope.vars:
            raise Refuse("unknown_call", "QuantumCircuit is not imported from qiskit", st)
        if len(call.args) != 1 or call.keywords:
            raise Refuse("unsupported_statement", "QuantumCircuit not called with one qubit count", st)
        nq = eval_value(call.args[0], scope)
        if not isinstance(nq, Lin) or nq.ki:
            raise Refuse("nonaffine_expression", "the qubit count is not affine in n", st)
        if nq.kn == 0:
            raise Refuse("no_size_dependence",
                         "the qubit count does not depend on the first parameter (fixed size)", st)
        if nq.kn != 1:
            raise Refuse("two_registers", f"qubit count {nq.kn}n{nq.c:+d}: several registers of length n", st)
        self.nq, self.qc_name = nq, name
        scope.bind(name, QC(), st, self.ctx[-1])

    # ---- calls
    def call(self, call: ast.Call, scope: Scope, conds):
        f = call.func
        if isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name):
            try:
                obj = _deref(scope.lookup(f.value.id, f.value))
            except Refuse:
                obj = None
            if isinstance(obj, QC):
                return self.gate(f.attr, call, scope, conds)
            raise Refuse("unknown_call", f"method `{ast.unparse(f)}` of something other than the circuit", call)
        if isinstance(f, ast.Name):
            if f.id in scope.vars:
                raise Refuse("unknown_call", f"call of a local value `{f.id}`", call)
            res = scope.module.resolve_func(f.id)
            if res is None:
                raise Refuse("unknown_call", f"call of `{f.id}` outside the accepted fragment", call)
            return self.inline(res[0], res[1], call, scope, conds)
        raise Refuse("unknown_call", f"call of `{ast.unparse(f)}`", call)

    def inline(self, module, fdef, call, scope, conds):
        key = (module.path, fdef.name)
        if key in self.stack:
            raise Refuse("recursion", f"recursive call of `{fdef.name}`", call)
        a = fdef.args
        if a.vararg or a.kwarg or a.kwonlyargs or a.posonlyargs or fdef.decorator_list:
            raise Refuse("unsupported_statement", f"signature of `{fdef.name}` outside the accepted fragment", call)
        names = [x.arg for x in a.args]
        binding = {}
        if len(call.args) > len(names):
            raise Refuse("unsupported_statement", f"too many arguments for `{fdef.name}`", call)
        for nm, arg in zip(names, call.args):
            binding[nm] = Thunk(arg, scope)
        for kw in call.keywords:
            if kw.arg is None or kw.arg not in names or kw.arg in binding:
                raise Refuse("unsupported_statement", f"arguments of `{fdef.name}`", call)
            binding[kw.arg] = Thunk(kw.value, scope)
        defaults = dict(zip(names[len(names) - len(a.defaults):], a.defaults))
        for nm in names:
            if nm not in binding:
                if nm not in defaults:
                    raise Refuse("unsupported_statement", f"missing argument `{nm}`", call)
                binding[nm] = Thunk(defaults[nm], Scope(module, {}, None, bounds=scope.bounds))
        body = list(fdef.body)
        if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
            body = body[1:]
        # the call is a statement, so a final `return` is a no-op: bare, or the name of a parameter (`return qc`
        # of a helper that changes the circuit in place), whose lookup cannot fail and whose value is dropped
        if body and isinstance(body[-1], ast.Return) and (body[-1].value is None or isinstance(body[-1].value, ast.Name)
                                                          and body[-1].value.id in names):
            body = body[:-1]
        c = self.new_ctx()
        inner = Scope(module, binding, None, ctx=c, base_conds=tuple(conds), bounds=scope.bounds)
        self.stack.append(key)
        self.ctx.append(c)
        try:
            self.block(body, inner, conds)
        except Refuse as r:
            if r.file is None:
                r.file = module.path.name
            r.via.append(f"{fdef.name} <- {scope.module.path.name}:{call.lineno}")
            raise
        finally:
            self.stack.pop()
            self.ctx.pop()

    def new_ctx(self):
        self.next_ctx += 1
        return self.next_ctx

    def gate(self, meth, call, scope, conds):
        if meth in IGNORED:
            return
        if meth in GROWING:
            raise Refuse("growing_arity", f"`{meth}`: multi-controlled gate of growing arity", call)
        if meth in LIBRARY:
            raise Refuse("library_gate", f"`{meth}`: library block or circuit transformation", call)
        if meth not in GATES:
            raise Refuse("unknown_gate", f"gate `{meth}` outside the gate set", call)
        if call.keywords:
            raise Refuse("unsupported_statement", f"keyword arguments of `{meth}`", call)
        nqb, nang = GATES[meth]
        if len(call.args) != nqb + nang:
            raise Refuse("unsupported_statement", f"`{meth}`: unexpected number of arguments", call)
        angle = eval_angle(call.args[0], scope) if nang else None
        qargs = [eval_value(x, scope) for x in call.args[nang:]]
        if nqb == 1 and isinstance(qargs[0], (Rng, Arr, ListVal)):
            return self.broadcast(meth, angle, qargs[0], call, scope, conds)
        for q in qargs:
            if not isinstance(q, Lin):
                raise Refuse("broadcast_multi", f"`{meth}` on a list of qubits", call)
        self.check_conds(conds, call)
        self.col.emit(Gate(meth, tuple(qargs), angle, tuple(conds), call.lineno))

    def broadcast(self, meth, angle, q, call, scope, conds):
        if isinstance(q, ListVal):
            for item in q.items:
                if not isinstance(item, Lin):
                    raise Refuse("broadcast_multi", "nested list of qubits", call)
                self.col.emit(Gate(meth, (item,), angle, tuple(conds), call.lineno))
            return
        rng = q if isinstance(q, Rng) else Rng(q.start, q.start + q.length, 1)
        self.iterate(rng, "<broadcast>", call, conds,
                     lambda sc, cs: self.col.emit(Gate(meth, (sc.vars["<broadcast>"],), angle,
                                                       tuple(cs), call.lineno)),
                     scope)

    def check_conds(self, conds, node):
        for b in conds:
            if b.idx.ki == -1:
                raise Refuse("fold_index", "bit index decreases with i (n-1-i): a fold, not a ladder", node)
            if b.idx.kn and b.idx.ki:
                raise Refuse("two_registers", "bit index depends on both n and i", node)
            if b.idx.ki not in (0, 1) or b.idx.kn not in (0, 1):
                raise Refuse("nonaffine_index", "bit index not of the form c, n+c, i+c", node)

    # ---- loops
    def for_loop(self, st: ast.For, scope: Scope, conds):
        if st.orelse:
            raise Refuse("unsupported_statement", "for ... else", st)
        if not isinstance(st.target, ast.Name):
            raise Refuse("unsupported_statement", "unpacking in a for header", st)
        it = eval_value(st.iter, scope)
        if isinstance(it, Arr):
            it = Rng(it.start, it.start + it.length, 1)
        if isinstance(it, ListVal):
            raise Refuse("range_not_ladder", "loop over an explicit list", st)
        if not isinstance(it, Rng):
            raise Refuse("range_not_ladder", "loop not over a range", st)
        var = st.target.id
        self.iterate(it, var, st, conds, lambda sc, cs: self.block(st.body, sc, cs), scope)

    def iterate(self, rng: Rng, var, node, conds, body_fn, scope):
        for v in (rng.start, rng.stop):
            if v.ki:
                raise Refuse("nested_loop", "loop bounds depend on an outer counter", node)
            if v.kn not in (0, 1):
                raise Refuse("range_not_ladder", "loop bound not of the form c or n+c", node)
        if var in scope.vars and not isinstance(scope.vars[var], Poison):
            raise Refuse("alias_reassigned", f"counter `{var}` shadows a name", node)
        before = set(scope.vars)
        local = set()

        def run(value):
            scope.vars[var] = value
            scope.ctx_of[var] = LOOP_VAR_CTX
            self.ctx.append(self.new_ctx())
            try:
                body_fn(scope, conds)
            finally:
                self.ctx.pop()
                for k in set(scope.vars) - before:
                    local.add(k)
                    del scope.vars[k]

        def poison():
            for k in local | {var}:
                scope.vars[k] = Poison(f"`{k}` from the loop body is used after the loop")
                scope.ctx_of[k] = 0
        if rng.start.kn == rng.stop.kn:
            vals = list(range(rng.start.c, rng.stop.c, rng.step))
            if len(vals) > MAX_UNROLL:
                raise Refuse("range_not_ladder", "bounded loop is too long", node)
            for v in vals:
                run(Lin(v, rng.start.kn, 0))
            poison()
            return
        up = rng.step > 0 and rng.start.kn == 0 and rng.stop.kn == 1
        down = rng.step < 0 and rng.start.kn == 1 and rng.stop.kn == 0
        if not (up or down):
            raise Refuse("range_not_ladder", "range is empty for large n", node)
        if self.col.body is not None:
            raise Refuse("nested_loop", "nested ladder (a loop up to n inside a loop up to n)", node)
        self.col.body = []
        k = len(self.bounds)
        run(Lin(0, 0, 1))
        poison()
        self.loop_bounds += [(v, nd, rng) for v, nd in self.bounds[k:]]
        del self.bounds[k:]
        body, self.col.body = self.col.body, None
        self.col.segments.append({"kind": "loop", "range": rng, "body": body,
                                  "line": getattr(node, "lineno", None)})

    # ---- result
    def qref(self, v: Lin, node_line):
        if v.ki == -1:
            raise Refuse("fold_index", "qubit index decreases with i (n-1-i): a fold, not a ladder",
                         _L(node_line))
        if v.kn and v.ki:
            raise Refuse("two_registers", "qubit index depends on both n and i (a second register)",
                         _L(node_line))
        if v.ki not in (0, 1):
            raise Refuse("stride_index", f"qubit index {v.ki}*i+c: the window stride is inside the index",
                         _L(node_line))
        if v.kn not in (0, 1):
            raise Refuse("two_registers", f"qubit index {v.kn}n+c", _L(node_line))
        return {"ref": "i" if v.ki else ("n" if v.kn else "const"), "off": v.c}

    def finish(self):
        e = self.nq.c
        nmin = [1]
        segs_out, loops_summary = [], []
        consts, ends = [], []                           # positions c and n + c of the current segment

        def p3():                                       # (P3): in a segment, positions c (qubits and bits) lie left of n + c
            if consts and ends:
                nmin.append(max(consts) - min(ends) + 1)
            consts.clear()
            ends.clear()

        def chk_fixed_index(v, line):
            r = self.qref(v, line)
            if r["ref"] == "const":
                if v.c < 0:
                    raise Refuse("index_out_of_range", f"negative qubit index {v.c}", _L(line))
                nmin.append(v.c - e + 1)
                consts.append(v.c)
            elif r["ref"] == "n":
                if v.c > e - 1:
                    raise Refuse("index_out_of_range", f"index n{v.c:+d} with {self.nq.c:+d} extra qubits",
                                 _L(line))
                nmin.append(-v.c)
                ends.append(v.c)
            return r

        def pair_fixed(a, b, line):
            if a.kn == b.kn:
                if a.c == b.c:
                    raise Refuse("qubit_collision", "gate applied twice to one qubit", _L(line))
            else:
                c, o = (a.c, b.c) if a.kn == 0 else (b.c, a.c)
                nmin.append(c - o + 1)

        def gate_json(g: Gate, refs):
            d = {"op": g.op, "qubits": refs, "line": g.line}
            if g.angle is not None and g.angle.rat != 0:     # (A3): the program keeps q*pi only, so r != 0 is refused (Sec. S2.2)
                raise Refuse("nonconst_angle", f"angle {g.angle.pi}*pi + {g.angle.rat}: the remainder is not a multiple of pi", _L(g.line))
            if g.angle is not None:
                d["angle"] = g.angle.to_json()
            if g.cond:
                d["cond"] = [{"param": b.param, "index": self.qref(b.idx, g.line), "access": b.access,
                              "value": b.value} for b in g.cond]
            return d

        def bit_below_n(c, kn, line):                   # bit index c + kn * n is at most n - 1
            if kn == 0:
                nmin.append(c + 1)
            elif kn != 1 or c > -1:
                raise Refuse("bit_index_out_of_range", "bit index is not below n: a has no such bit",
                             _L(line))

        for seg in self.col.segments:
            if seg["kind"] == "fixed":
                gates = []
                for g in seg["gates"]:
                    refs = [chk_fixed_index(q, g.line) for q in g.qubits]
                    for x in range(len(g.qubits)):
                        for y in range(x + 1, len(g.qubits)):
                            pair_fixed(g.qubits[x], g.qubits[y], g.line)
                    for b in g.cond:
                        if b.idx.kn == 0 and b.idx.c < 0:
                            raise Refuse("bit_index_negative", "negative bit index", _L(g.line))
                        if b.idx.kn and b.idx.c < 0:
                            nmin.append(-b.idx.c)
                        bit_below_n(b.idx.c, b.idx.kn, g.line)
                        (ends if b.idx.kn else consts).append(b.idx.c)
                    gates.append(gate_json(g, refs))
                segs_out.append({"kind": "fixed", "gates": gates})
                p3()
                continue
            rng: Rng = seg["range"]
            up = rng.step > 0
            st = abs(rng.step)
            # extreme counter values: lo_c is the constant end, the other end is n + hi_c
            if up:
                lo_c, hi_c = rng.start.c, rng.stop.c - 1
                nmin.append(lo_c - hi_c)          # nonempty: n + hi_c >= lo_c
            else:
                lo_c, hi_c = rng.stop.c + 1, rng.start.c
                nmin.append(lo_c - hi_c)
            win, fixed_q, fixed_bits, bitdep, body = set(), [], [], False, []    # the window counts qubits and bits
            for g in seg["body"]:
                refs = []
                for q in g.qubits:
                    r = self.qref(q, g.line)
                    if r["ref"] == "i":
                        if lo_c + q.c < 0:
                            raise Refuse("index_out_of_range", "index i+c is negative at the end of the loop",
                                         _L(g.line))
                        if hi_c + q.c > e - 1:
                            raise Refuse("index_out_of_range", "index i+c exceeds the number of qubits",
                                         _L(g.line))
                        win.add(q.c)
                    else:
                        chk_fixed_index(q, g.line)
                        key = (r["ref"], r["off"])
                        if key not in fixed_q:
                            fixed_q.append(key)
                    refs.append(r)
                qs = g.qubits
                for x in range(len(qs)):
                    for y in range(x + 1, len(qs)):
                        a, b = qs[x], qs[y]
                        if a.ki == 0 and b.ki == 0:
                            pair_fixed(a, b, g.line)
                            continue
                        if a.ki == 0:
                            a, b = b, a
                        if b.ki == 1:
                            if a.c == b.c:
                                raise Refuse("qubit_collision", "gate applied twice to one qubit", _L(g.line))
                        elif b.kn == 0:       # i + a.c == b.c
                            istar = b.c - a.c
                            hit = istar >= lo_c and (st == 1 or not up or (istar - lo_c) % st == 0)
                            if hit:
                                raise Refuse("qubit_collision",
                                             f"the window hits the fixed qubit {b.c} at i={istar}", _L(g.line))
                        else:                 # i + a.c == n + b.c
                            if b.c - a.c <= hi_c:
                                raise Refuse("qubit_collision",
                                             f"the window hits the fixed qubit n{b.c:+d} at the right end", _L(g.line))
                for b in g.cond:
                    if b.idx.ki:
                        bitdep = True
                        if lo_c + b.idx.c < 0:
                            raise Refuse("bit_index_negative", "bit index i+c is negative at the end of the loop",
                                         _L(g.line))
                        win.add(b.idx.c)
                    elif b.idx.kn == 0 and b.idx.c < 0:
                        raise Refuse("bit_index_negative", "negative bit index", _L(g.line))
                    elif b.idx.kn == 0 and not up:
                        raise Refuse("bit_const_in_down_ladder", "condition on a bit a[c] in a loop with step -s: after "
                                     "mirroring it is a bit at the end of the register (Definition 3(b))", _L(g.line))
                    elif b.idx.kn == 0:
                        consts.append(b.idx.c)
                        if b.idx.c not in fixed_bits:
                            fixed_bits.append(b.idx.c)
                    elif b.idx.kn:
                        raise Refuse("bit_on_end_fixed", "condition on a fixed bit a[n+c] at the end of the register: "
                                     "Definition 3(b) does not allow it", _L(g.line))
                    bit_below_n(b.idx.c + b.idx.ki * hi_c, b.idx.kn + b.idx.ki, g.line)   # at i = n + hi_c
                body.append(gate_json(g, refs))
            if not win:
                raise Refuse("body_without_window",
                             "the loop body does not depend on the counter: a repetition, not a ladder", _L(seg["line"]))
            # no window position may meet a fixed position at any i (even in different gates of the body)
            for href, hoff in fixed_q:
                for o in sorted(win):
                    if href == "const":
                        istar = hoff - o
                        hit = istar >= lo_c and (st == 1 or not up or (istar - lo_c) % st == 0)
                    else:
                        istar, hit = None, hoff - o <= hi_c
                    if hit:
                        where = f"at i={istar}" if istar is not None else "at the right end"
                        fixed = str(hoff) if href == "const" else f"n{hoff:+d}"
                        raise Refuse("fixed_overlaps_window",
                                     f"window position i{o:+d} meets the fixed position {fixed} {where}", _L(seg["line"]))
            for href, hoff in fixed_q:               # a start fixed position lies left of the first window, even when no window meets it
                if href == "const" and hoff >= lo_c + min(win):
                    raise Refuse("fixed_inside_ladder", f"fixed qubit {hoff} inside the ladder region (in a gap between "
                                 "windows): Definition 3(c) needs start fixed positions left of the first window",
                                 _L(seg["line"]))
            for c in fixed_bits:                    # a bit a[c] is a position of S, also left of the first window
                if c >= lo_c + min(win):
                    raise Refuse("fixed_inside_ladder", f"bit a[{c}] inside the ladder region: positions of S, including "
                                 "bits a[c], must lie left of the first window (Definition 3(c))", _L(seg["line"]))
            p3()
            w = max(win) - min(win) + 1
            fixed_json = [{"ref": r, "off": o} for r, o in fixed_q]
            seg_json = {
                "kind": "loop", "line": seg["line"],
                "range": {"start": _q(rng.start), "stop": _q(rng.stop), "step": rng.step},
                "order": "up" if up else "down",
                "i_first": {"ref": "const" if up else "n", "off": rng.start.c},
                "i_bound": {"ref": "n" if up else "const", "off": hi_c if up else lo_c},
                "s": st, "w": w, "h": len(fixed_q) + len(fixed_bits),
                "window": {"base": min(win), "offsets": sorted(win)},
                "fixed_qubits": fixed_json,
                "fixed_bits": fixed_bits,
                "bit_dependent": bitdep,
                "body": body,
            }
            segs_out.append(seg_json)
            loops_summary.append({"s": st, "w": w, "h": len(fixed_q) + len(fixed_bits), "fixed_qubits": fixed_json,
                                  "bit_dependent": bitdep})
        def nonneg(c0, kn, node):                       # c0 + kn * n >= 0 for all n >= n_min
            if kn > 0:
                nmin.append(-(c0 // kn))
            elif c0 < 0 or kn < 0:
                raise Refuse("index_out_of_range", "index or slice outside the list of qubits", node)
        for v, nd in self.bounds:
            nonneg(v.c, v.kn, nd)
        for v, nd, rng in self.loop_bounds:            # linear in i: the two ends of the counter range suffice
            lo, hi = (rng.start.c, rng.stop.c - 1) if rng.step > 0 else (rng.stop.c + 1, rng.start.c)
            nonneg(v.c + v.ki * lo, v.kn, nd)           # i = lo
            nonneg(v.c + v.ki * hi, v.kn + v.ki, nd)    # i = n + hi
        if not loops_summary and not any(q["ref"] == "n" for sg in segs_out for g in sg["gates"] for q in g["qubits"]):
            raise Refuse("no_ladder", "no loop and no qubit n + c: the code has a fixed size", self.fn)
        ops = {g["op"] for s in segs_out for g in (s.get("gates") or s.get("body"))}
        return {
            "schema": SCHEMA,
            "status": "ladder",
            "size_param": self.size_param,
            "bit_params": self.bit_params,
            "num_qubits": {"ref": "n", "off": e},
            "n_min": max(nmin),
            "gate_set": "core" if ops <= CORE_GATES else "extended",
            "ops": sorted(ops),
            "loops": loops_summary,
            "segments": segs_out,
        }


def _check_recursion(mod: Module, fdef: ast.FunctionDef):
    """call graph of the module functions reachable from the generator; a cycle is refused"""
    state = {}

    def visit(m, f, path):
        key = (m.path, f.name)
        if state.get(key) == "open":
            raise Refuse("recursion", "recursion: " + " -> ".join(path + [f.name]), f)
        if state.get(key) == "done":
            return
        state[key] = "open"
        for sub in ast.walk(f):
            if isinstance(sub, ast.Call) and isinstance(sub.func, ast.Name):
                res = m.resolve_func(sub.func.id)
                if res:
                    visit(res[0], res[1], path + [f.name])
        state[key] = "done"
    visit(mod, fdef, [])


def _q(v: Lin):
    return {"ref": "n" if v.kn else "const", "off": v.c}


class _L:
    def __init__(self, line):
        self.lineno = line


# ---------------------------------------------------------------- interface

def analyze(path, func: str) -> dict:
    path = Path(path)
    src = {"file": str(path), "function": func}
    try:
        mod = Module.load(path)
        if func not in mod.funcs:
            raise Refuse("unknown_call", f"function `{func}` not found", None)
        if mod.shadowed & (BUILTINS_USED | {"QuantumCircuit"}):
            raise Refuse("shadowed_builtin", f"the module redefines {sorted(mod.shadowed)}", None)
        _check_recursion(mod, mod.funcs[func])
        out = Analyzer(mod, mod.funcs[func]).run()
        return {"schema": SCHEMA, "source": src, **{k: v for k, v in out.items() if k != "schema"}}
    except Refuse as r:
        reason = {"code": r.code, "message": r.msg, "line": r.line}
        if r.file:
            reason["file"] = r.file
        if r.via:
            reason["via"] = r.via
        return {"schema": SCHEMA, "source": src, "status": "refused", "reason": reason}
    except RecursionError:
        return {"schema": SCHEMA, "source": src, "status": "refused",
                "reason": {"code": "recursion", "message": "inlining too deep", "line": None}}


# ---------------------------------------------------------------- expansion of a description (checks against Qiskit)

def ref_value(r: dict, n: int, i: int | None = None) -> int:
    base = {"const": 0, "n": n, "i": i}[r["ref"]]
    return base + r["off"]


def _bit(params, c, n, i):
    v = params[c["param"]]
    k = ref_value(c["index"], n, i)
    b = v[k] if c["access"] == "index" else (v >> k) & 1
    return int(bool(b)) == c["value"]


def expand(desc: dict, n: int, params: dict | None = None):
    """gate list (op, qubits, angle or None) for given n and parameters"""
    params = params or {}
    out = []

    def emit(g, i):
        if all(_bit(params, c, n, i) for c in g.get("cond", [])):
            ang = g["angle"]["value"] if "angle" in g else None
            out.append((g["op"], tuple(ref_value(q, n, i) for q in g["qubits"]), ang))

    for seg in desc["segments"]:
        if seg["kind"] == "fixed":
            for g in seg["gates"]:
                emit(g, None)
        else:
            r = seg["range"]
            for i in range(ref_value(r["start"], n), ref_value(r["stop"], n), r["step"]):
                for g in seg["body"]:
                    emit(g, i)
    return out


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: python -m ladder_ec.translator FILE.py FUNCTION", file=sys.stderr)
        sys.exit(2)
    print(json.dumps(analyze(sys.argv[1], sys.argv[2]), ensure_ascii=False, indent=1))
