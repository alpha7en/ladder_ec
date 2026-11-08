"""Table 5: share of the unequal pairs on which an error is found, by detection group (Sec. 7.3): the instance
checks (state, random-input, sampled; experiments/instance_checks.py) and the method, with the stage that found
the counterexample. Also the effect of the quick search (run without it, data/runs/method_no_quick_search.jsonl).
Usage: python -I experiments/tables/table05_detection.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Benchmark, C, answer, kind, group, pct, median, by_id, run_path      # noqa: E402

NAMES = {1: "1. at small n", 2: "2. only on a superposition", 3: "3. mutation of a condition on bits",
         4: "4. only at 10 < n <= 20", 5: "5. only at n > 20"}
CHECKS = ("judge_state", "judge_unitary", "random_test")       # columns State, Rand., Samp. (instance_checks.py)


def main():
    B = Benchmark()
    M, NQ = B.run("method_main.jsonl"), B.run("method_no_quick_search.jsonl")
    base = by_id(run_path("instance_checks.jsonl"))
    bad = B.unequal()
    L = lambda i: B.pairs[i]["label"]
    G = {}
    for i in bad:
        G.setdefault(group(B.pairs[i], B.truth[i]), []).append(i)
    print("%-36s %5s %6s %6s %6s %6s %9s %9s" % ("Group", "Pairs", "State", "Rand.", "Samp.", "Found", "Automata",
                                                 "Refusal,timeout"))
    for g in sorted(G):
        ids = G[g]
        found = [i for i in ids if kind(answer(M[i], L(i))) == "neq"]
        fb = C(M[i]["counterexamples"][L(i)]["found_by"] for i in found)
        print("%-36s %5d %6s %6s %6s %6s %9d %9s   %s" % (
            NAMES[g], len(ids), *(pct(sum(base[i][k]["caught"] for i in ids), len(ids)) for k in CHECKS),
            pct(len(found), len(ids)), fb["automata"], pct(len(ids) - len(found), len(ids)), dict(fb)))
    found = [i for i in bad if kind(answer(M[i], L(i))) == "neq"]
    print("counterexamples %d, found by %s" % (len(found), dict(C(M[i]["counterexamples"][L(i)]["found_by"]
                                                                   for i in found))))
    print("without the quick search: found %d; other answers on the pairs found with it: %s" % (
        sum(kind(answer(NQ[i], L(i))) == "neq" for i in bad),
        dict(C(kind(answer(NQ[i], L(i))) for i in found if kind(answer(NQ[i], L(i))) != "neq"))))
    qs = [i for i in found if M[i]["counterexamples"][L(i)]["found_by"] == "quick_search"]
    print("  %d pairs found by the quick search: median time %.3f s with it, %.3f s without" % (
        len(qs), median(M[i]["time_s"] for i in qs), median(NQ[i]["time_s"] for i in qs)))
    eq = [i for i in B.equal() if kind(answer(M[i], L(i))) in ("eq_exact", "eq_eps")]
    print("  %d equal pairs answered \"equal\": median time %.3f s with it, %.3f s without" % (
        len(eq), median(M[i]["time_s"] for i in eq), median(NQ[i]["time_s"] for i in eq)))


if __name__ == "__main__":
    main()
