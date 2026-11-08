"""QCEC enumeration of sizes on the benchmark (Sec. 7.4, Tables 6 and 7, Fig. 6), with experiments/qcec/driver.py
in mode A: all n from 1 to the bound, for tasks with bits all a (at most 4096 per n); 60 s per instance and per
pair; 3 GB per worker process; one pair at a time (on the machine of the method).
Modes (first argument):
  ncut      equal pairs answered "equal": the bound n_cut the method reported (data/runs/method_main.jsonl),
            capped at 10^5;
  ncut_min  the same pairs with the bound of Corollary 1 and Table 1 computed from the dimensions d_G, d_R of the
            minimized family automata (modulo one prime) instead of the dimension bounds of the method:
            d = d_G + d_R (E1, E5, E3), 2 d_G d_R (E2), d_G^2 + d_R^2 (E5phase, E3phase); n_cut = n_0 + d.
            This is the smallest bound Corollary 1 gives while equality is unknown. The dimensions come from an
            instrumented run of the method; the input file with these bounds is shipped as
            data/runs/qcec_ncut_min_input.jsonl and used as it is (only pairs whose bound is below n_cut);
  witness   unequal pairs: the bound is the least differing size n_d of the ground truth, an enumeration that
            knows where to stop (a lower bound on the time of any enumeration).
Usage: python -I experiments/qcec/enumerate_sizes.py ncut|ncut_min|witness
       -> data/runs/qcec_<mode>.jsonl (per pair), data/runs/qcec_<mode>_instances.jsonl (per instance, not shipped)"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import driver as D   # noqa: E402

BENCH, RUNS = os.path.join(ROOT, "data", "benchmark"), os.path.join(ROOT, "data", "runs")
D.PAIR_LIMIT = 60.0
D.CONFIG["limit_pair_s"] = D.PAIR_LIMIT
D.MEM_LIMIT = 3 << 30
CUT_MAX = 10 ** 5


def write_input(mode, inp):
    pairs = [json.loads(l) for l in open(os.path.join(BENCH, "pairs.jsonl"))]
    truth = {r["id"]: r for r in map(json.loads, open(os.path.join(BENCH, "truth.jsonl")))}
    ours = {r["id"]: r["ours"] for r in map(json.loads, open(os.path.join(RUNS, "method_main.jsonl")))}
    with open(inp, "w") as f:
        for p in pairs:
            t, o, L = truth.get(p["id"]), ours.get(p["id"]), p["label"]
            if not t or t["noop"]:
                continue
            q = {k: p[k] for k in ("id", "task", "label", "ref", "sol", "origin")}
            if mode == "ncut":
                if not t["truth"][L] or not isinstance(o["verdicts"], dict) or \
                        not o["verdicts"].get(L, "").startswith("equal") or L not in o.get("n_cut", {}):
                    continue
                cut = o["n_cut"][L]
                q.update(cutoff={L: min(cut, CUT_MAX)}, n_cut_true=cut, n_min=1)
            else:
                if t["truth"][L] is not False or not t["first_difference"].get(L):
                    continue
                q.update(cutoff={L: t["first_difference"][L]["n"]}, witness_n=t["first_difference"][L]["n"], n_min=1)
            f.write(json.dumps(q, ensure_ascii=False) + "\n")


def main():
    mode = sys.argv[1]
    assert mode in ("ncut", "ncut_min", "witness"), mode
    D.CONFIG["tool"] += ("; run from experiments/qcec/enumerate_sizes.py %s (pair limit 60 s, worker memory 3 GB, "
                         "one pair at a time)" % mode)
    inp = os.path.join(RUNS, "qcec_%s_input.jsonl" % mode)
    if mode != "ncut_min":
        write_input(mode, inp)
    sys.argv = [sys.argv[0], "--mode", "A", "--pairs", inp, "--resume",
                "--out", os.path.join(RUNS, "qcec_%s.jsonl" % mode),
                "--inst", os.path.join(RUNS, "qcec_%s_instances.jsonl" % mode), "--jobs", "1"] + sys.argv[2:]
    D.main()


if __name__ == "__main__":
    main()
