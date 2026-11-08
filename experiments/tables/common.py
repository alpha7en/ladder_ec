"""Shared loaders and counting rules of the Section 7 tables. Only reads files under data/; runs nothing."""
import collections
import json
import os
import re
import statistics

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BENCH = os.path.join(ROOT, "data", "benchmark")
RUNS = os.path.join(ROOT, "data", "runs")
C = collections.Counter
median = statistics.median


def jl(path):
    return [json.loads(l) for l in open(path)]


def by_id(path, key=None):
    return {r["id"]: (r[key] if key else r) for r in jl(path)}


def run_path(name):
    """data/runs/<name>; the environment variable RUNS=<dir> reads the run files from another folder"""
    return os.path.join(os.environ.get("RUNS", RUNS), name)


class Benchmark:
    """pairs, ground truth and the method runs; pairs whose mutation is a no-op are left out (Sec. 7.1)"""

    def __init__(self):
        self.truth = by_id(os.path.join(BENCH, "truth.jsonl"))
        self.pairs_all = by_id(os.path.join(BENCH, "pairs.jsonl"))
        self.pairs = {i: p for i, p in self.pairs_all.items() if not self.truth[i]["noop"]}
        self.second = by_id(os.path.join(BENCH, "pairs_second.jsonl"))

    def run(self, name):
        return by_id(run_path(name), "ours")

    def truth_of(self, i, label):
        """ground truth of pair i for a relation key (the truth file writes E5 for E5strict)"""
        return self.truth[i]["truth"][{"E5strict": "E5"}.get(label, label)]

    def equal(self):
        return [i for i, p in self.pairs.items() if self.truth_of(i, p["label"])]

    def unequal(self):
        return [i for i, p in self.pairs.items() if self.truth_of(i, p["label"]) is False]


def answer(o, label):
    v = o["verdicts"]
    return v.get(label, "") if isinstance(v, dict) else str(v)


def kind(s):
    """answer string -> one of eq_exact, eq_eps, noans, neq, timeout, refused"""
    for prefix, k in (("equal, exact", "eq_exact"), ("equal, error", "eq_eps"), ("no answer", "noans"),
                      ("NOT", "neq"), ("timeout", "timeout"), ("refused", "refused")):
        if s.startswith(prefix):
            return k
    return "other:" + s


def eps_exponent(s):
    return int(re.search(r"2\^(-?\d+)", s).group(1))


def group(p, truth):
    """detection group of an unequal pair (Table 5): by the least differing size n_d (ground truth) and the kind
    of difference the benchmark generator recorded for the mutation (pair field expect.when)"""
    n = truth["first_difference"][p["label"]]["n"]
    if n > 20:
        return 5
    if n > 10:
        return 4
    if p["expect"]["when"] == "superposition":
        return 2
    if p["expect"]["when"] in ("some a", "all a", "asymmetric a", "a = 1...1"):
        return 3
    return 1


# Kind of each refusal code (Table S2, column "Kind"): out of class, out of fragment, or both.
KIND = {"branch_on_size": "class", "branch_on_loop_index": "class", "branch_on_param": "class",
        "conditional_expression": "both", "nonaffine_expression": "both", "nonaffine_index": "both",
        "nested_loop": "class", "while_loop": "both", "recursion": "both", "growing_arity": "class",
        "library_gate": "fragment", "nonconst_angle": "class", "two_registers": "class", "fold_index": "class",
        "stride_index": "class", "index_out_of_range": "class", "qubit_collision": "class",
        "bit_index_negative": "class", "bit_index_out_of_range": "class", "bit_on_end_fixed": "class",
        "bit_const_in_down_ladder": "class", "fixed_overlaps_window": "class", "fixed_inside_ladder": "class",
        "alias_reassigned": "fragment", "alias_under_condition": "fragment", "loop_var_after_loop": "fragment",
        "bit_param_misuse": "both", "range_not_ladder": "both", "body_without_window": "class",
        "no_size_dependence": "both", "no_ladder": "fragment", "unknown_call": "fragment", "unknown_name": "fragment",
        "unknown_gate": "fragment", "unsupported_statement": "fragment", "unsupported_expression": "fragment",
        "extra_params": "fragment", "shadowed_builtin": "fragment",
        "early_return": "fragment"}       # a return before the end of the function: a syntax restriction


def pct(k, m):
    return "%d%%" % round(100 * k / m) if m else "-"
