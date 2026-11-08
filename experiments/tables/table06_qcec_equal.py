"""Table 6: equal pairs answered "equal" (444): QCEC enumeration of all sizes n <= n_cut (Corollary 1) and of all
n up to the smaller bound of Table 1 computed from the minimized automata, against the method (Sec. 7.4).
Cells "x / y": enumeration up to n_cut / up to the smaller bound. An enumeration not finished within 60 s counts
as 60 s. Usage: python -I experiments/tables/table06_qcec_equal.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Benchmark, answer, kind, median, by_id, run_path      # noqa: E402


def main():
    B = Benchmark()
    M = B.run("method_main.jsonl")
    nc, mn = by_id(run_path("qcec_ncut.jsonl")), by_id(run_path("qcec_ncut_min.jsonl"))
    L = lambda i: B.pairs[i]["label"]
    eq = {i for i in B.equal() if kind(answer(M[i], L(i))) in ("eq_exact", "eq_eps")}
    assert set(nc) == eq, (len(nc), len(eq))
    diff = [(i, M[i]["n_cut"][L(i)], nc[i]["cutoff"]) for i in nc if M[i]["n_cut"][L(i)] != nc[i]["cutoff"]]
    if diff:
        print("QCEC enumeration capped at n = 10^5 (n_cut of the method, cap):", diff)
    row = lambda i: ("with bits a" if nc[i]["bits"] else "no bits, n_cut <= 50" if nc[i]["cutoff"] <= 50
                     else "no bits, 51..1000" if nc[i]["cutoff"] <= 1000 else "no bits, > 1000")
    t = lambda r: r["wall_s"] if r["verdict"] == "equal" else 60.0
    print("%-22s %-22s %5s %13s %13s %15s" % ("Tasks", "n_cut", "Pairs", "QCEC finished", "Method median",
                                             "QCEC median"))
    for k in ("no bits, n_cut <= 50", "no bits, 51..1000", "no bits, > 1000", "with bits a"):
        ids = [i for i in nc if row(i) == k]
        best = {i: mn.get(i, nc[i]) for i in ids}
        print("%-22s %-22s %5d %6d / %-6d %13.2f %7.1f / %-6.1f" % (
            k, "%d..%.3g" % (min(M[i]["n_cut"][L(i)] for i in ids), max(M[i]["n_cut"][L(i)] for i in ids)), len(ids),
            sum(nc[i]["verdict"] == "equal" for i in ids), sum(best[i]["verdict"] == "equal" for i in ids),
            median(M[i]["time_s"] for i in ids), median(t(nc[i]) for i in ids), median(t(best[i]) for i in ids)))
    print("smaller bound below n_cut: %d of %d pairs; finished: %d up to n_cut, %d with the smaller bound" % (
        len(mn), len(nc), sum(r["verdict"] == "equal" for r in nc.values()),
        sum(mn.get(i, nc[i])["verdict"] == "equal" for i in nc)))
    d2 = [i for i in mn if mn[i]["verdict"] == "equal"]
    print("finished with the smaller bound (%d pairs): QCEC median %.2f s, method %.3f s; QCEC faster on %d" % (
        len(d2), median(mn[i]["wall_s"] for i in d2), median(M[i]["time_s"] for i in d2),
        sum(mn[i]["wall_s"] < M[i]["time_s"] for i in d2)))


if __name__ == "__main__":
    main()
