"""Answers of the method on the AutoQ-Para pairs (Table 9): decide_core, three runs, the time is the median.
Usage: python -I experiments/autoq_para/run_ours.py -> data/runs/autoq_para/ours.jsonl"""
import json, os, statistics, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path[:0] = [HERE, ROOT]
from ladder_ec.check import decide_core   # noqa: E402
from pairs import PAIRS                   # noqa: E402

G = os.path.join(HERE, "gens.py")
decide_core((G, "ghz_ladder"), (G, "ghz_ladder"), ["E1"])          # warm-up (imports)
with open(os.path.join(ROOT, "data", "runs", "autoq_para", "ours.jsonl"), "w") as f:
    for name, r, s, rel, exp, aq, what in PAIRS:
        ts, res = [], None
        for _ in range(3):
            t = time.perf_counter()
            res = decide_core((G, r), (G, s), [rel])
            ts.append(time.perf_counter() - t)
        row = dict(pair=name, ref=r, sol=s, relation=rel, expected=exp, autoq=aq, verdict=res.get("verdicts"),
                   n_cut=res.get("n_cut"), n0=res.get("n0"), counterexamples=res.get("counterexamples"),
                   times_s=[round(x, 3) for x in ts], median_s=round(statistics.median(ts), 3),
                   ref_status=res.get("ref_status"))
        f.write(json.dumps(row, ensure_ascii=False, default=str) + "\n")
        print(name, row["verdict"], row["median_s"], flush=True)
