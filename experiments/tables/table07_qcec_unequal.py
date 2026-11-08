"""Table 7: unequal pairs (1717): answer of the method against QCEC on the sizes n <= n_d, the least differing size
of the ground truth (Sec. 7.4). "equivalent"* and "not checked"*: the solution raises an exception at n_d, so the
enumeration skips that size. Usage: python -I experiments/tables/table07_qcec_unequal.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Benchmark, C, answer, kind, median, by_id, run_path      # noqa: E402

COLS = [("not_equal", '"not equivalent"'), ("equal", '"equivalent"*'), ("error", "not checked*"),
        ("pair_timeout", "timeout 60 s")]
ROWS = [("neq", '"NOT equal"'), ("refused", "refusal"), ("timeout", "timeout 300 s")]


def q90(x):
    x = sorted(x)
    return x[int(0.9 * len(x))]


def main():
    B = Benchmark()
    M = B.run("method_main.jsonl")
    wit = by_id(run_path("qcec_witness.jsonl"))
    bad = B.unequal()
    assert set(wit) == set(bad)
    L = lambda i: B.pairs[i]["label"]
    tab = C((kind(answer(M[i], L(i))), wit[i]["verdict"]) for i in bad)
    print("%-16s" % "Method" + "".join("%18s" % c for _, c in COLS))
    for k, name in ROWS:
        print("%-16s" % name + "".join("%18d" % tab[(k, c)] for c, _ in COLS))
    rest = {key: v for key, v in tab.items() if key[0] not in dict(ROWS) or key[1] not in dict(COLS)}
    if rest:
        print("other:", rest)
    both = [i for i in bad if kind(answer(M[i], L(i))) == "neq" and wit[i]["verdict"] == "not_equal"]
    qt, mt = [wit[i]["wall_s"] for i in both], [M[i]["time_s"] for i in both]
    print("both find the error (%d): QCEC median %.3f s, 90%% %.3f s; method median %.3f s, 90%% %.3f s; "
          "QCEC faster on %d" % (len(both), median(qt), q90(qt), median(mt), q90(mt),
                                 sum(wit[i]["wall_s"] < M[i]["time_s"] for i in both)))
    print("source of the counterexample on these pairs:", dict(C(M[i]["counterexamples"][L(i)]["found_by"] for i in both)))


if __name__ == "__main__":
    main()
