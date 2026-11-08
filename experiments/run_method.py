"""Run the method on the whole benchmark in one process (Sec. 7.1, Tables 3-7, Fig. 6).

Variants:
  main             the relation of each task (data/benchmark/pairs.jsonl)       -> data/runs/method_main.jsonl
  second           the second relation of each task (pairs_second.jsonl)        -> data/runs/method_second.jsonl
  no_quick_search  the relation of each task, without the quick search          -> data/runs/method_no_quick_search.jsonl
Each pair gets 30 s; pairs over the limit are written as "timeout 30 s" and rerun with 300 s by --retry.
An interrupted run resumes: pairs already in the output are skipped. --seeds FILE takes the seed of the probabilistic plan of
each pair from an earlier output, so that the run repeats it (times aside).
Usage: python -I experiments/run_method.py main|second|no_quick_search [--retry] [--seeds FILE] [--out FILE]
"""
import argparse, json, os, signal, sys, time, traceback

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ladder_ec.check import decide_core      # noqa: E402

LIMIT, RETRY_LIMIT = 30, 300
PAIRS = {"main": "pairs.jsonl", "second": "pairs_second.jsonl", "no_quick_search": "pairs.jsonl"}


def timeout(*_):
    raise TimeoutError


def run(p, limit, quick, seed=None):
    t0 = time.time()
    signal.alarm(limit)
    try:
        ours = decide_core((os.path.join(ROOT, p["ref"]["file"]), p["ref"]["function"]),
                           (os.path.join(ROOT, p["sol"]["file"]), p["sol"]["function"]), [p["label"]],
                           seed=seed, quick_search=quick)
    except TimeoutError:
        ours = {"verdicts": "timeout %d s" % limit}
    except AssertionError:                       # an overflow guard (Sec. 6.2) or another internal check
        f = traceback.extract_tb(sys.exc_info()[2])[-1]
        ours = {"verdicts": "crash: AssertionError at %s:%d" % (os.path.basename(f.filename), f.lineno),
                "message": str(sys.exc_info()[1])[:300]}
    finally:
        signal.alarm(0)
    ours["time_s"] = round(time.time() - t0, 3)
    print(p["id"], ours.get("verdicts"), ours["time_s"], flush=True)
    return {"id": p["id"], "ours": ours}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("variant", choices=sorted(PAIRS))
    ap.add_argument("--retry", action="store_true", help="rerun the pairs that timed out, with 300 s")
    ap.add_argument("--seeds", help="earlier output whose seeds to reuse")
    ap.add_argument("--out", help="output file (default data/runs/method_<variant>.jsonl)")
    a = ap.parse_args()
    quick = a.variant != "no_quick_search"
    out = a.out or os.path.join(ROOT, "data", "runs", "method_%s.jsonl" % a.variant)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    seeds = {}
    if a.seeds:
        seeds = {r["id"]: r["ours"].get("seed") for r in map(json.loads, open(a.seeds))}
    signal.signal(signal.SIGALRM, timeout)
    pairs = [json.loads(l) for l in open(os.path.join(ROOT, "data", "benchmark", PAIRS[a.variant]))]
    if a.retry:
        rows = [json.loads(l) for l in open(out)]
        by_id = {p["id"]: p for p in pairs}
        rows = [run(by_id[r["id"]], RETRY_LIMIT, quick, seeds.get(r["id"]))
                if str(r["ours"]["verdicts"]).startswith("timeout") else r for r in rows]
        with open(out + ".tmp", "w") as f:
            f.writelines(json.dumps(r, ensure_ascii=False, default=str) + "\n" for r in rows)
        os.replace(out + ".tmp", out)
        return
    run(pairs[0], LIMIT, quick)                 # warm-up pair (imports, caches); not written
    done = {json.loads(l)["id"] for l in open(out)} if os.path.exists(out) else set()
    with open(out, "a") as f:
        for p in pairs:
            if p["id"] not in done:
                f.write(json.dumps(run(p, LIMIT, quick, seeds.get(p["id"])), ensure_ascii=False, default=str) + "\n")
                f.flush()


if __name__ == "__main__":
    main()
