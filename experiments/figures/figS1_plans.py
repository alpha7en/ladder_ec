"""Fig. S1: number of primes and time of the exact plan without the cap T <= T_max (Thm. S.12) and of the
probabilistic plan (Thm. S.13) against n_cut (Sec. 7.4, Online Resource 1).
Population: pairs without a no-op mutation on which the probabilistic plan answered "equal" with error <= 2^-40
(224 pairs of data/runs/exact_plan.jsonl, written by experiments/exact_plan.py). Times are those of the automata
part only. exact_time = plan_s + exact_s; for an exact plan not finished within 30 min, exact_done = 0 and
exact_time = plan_s + 1800 (a lower bound). ncut_j: n_cut spread by 10^u, |u| <= 0.04, for drawing only.
Writes experiments/figures/out/figS1_plans.csv and figS1_plans.pdf (when matplotlib is installed); prints the
numbers of Sec. 7.4. Usage: python -I experiments/figures/figS1_plans.py"""
import csv, json, os, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RUNS, BENCH, OUT, jitter, pyplot      # noqa: E402

LIMIT_S = 1800


def main():
    truth = {r["id"]: r for r in map(json.loads, open(os.path.join(BENCH, "truth.jsonl")))}
    rows = [r for r in map(json.loads, open(os.path.join(os.environ.get("RUNS", RUNS), "exact_plan.jsonl")))
            if r["id"] in truth and not truth[r["id"]]["noop"] and r["random"]["error"] is not None
            and r["random"]["error"] <= 2 ** -40]
    os.makedirs(OUT, exist_ok=True)
    out = []
    for r in sorted(rows, key=lambda r: r["id"]):
        nc, done = r["random"]["n_cut"], r["exact"] is not None
        out.append([r["id"], nc, "%.6g" % (nc * 10 ** jitter(r["id"])), r["random"]["primes"], r["T"], r["random_s"],
                    "%.3f" % (r["plan_s"] + (r["exact_s"] if r["exact_s"] is not None else LIMIT_S)), int(done)])
    with open(os.path.join(OUT, "figS1_plans.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "ncut", "ncut_j", "rand_primes", "T", "rand_s", "exact_time", "exact_done"])
        w.writerows(out)
    done = [r for r in rows if r["exact"] is not None]
    T = sorted(r["T"] for r in done)                         # statistics of the exact plan: finished pairs
    print("pairs %d; exact plan answers \"equal\" with error 0 on %d (difference found on %d)" % (
        len(rows), sum(r["exact"]["verdict"] == "equal" and r["exact"]["error"] == 0 for r in done),
        sum(r["exact"]["verdict"] != "equal" for r in done)))
    print("not finished within %d s: %s" % (LIMIT_S, [(r["id"], r["T"]) for r in rows if r["exact"] is None]))
    print("exact plan, primes (finished pairs): median %g, 90%% %d, max %d" % (
        statistics.median(T), T[int(0.9 * len(T))], max(T)))
    tf = sorted(r["plan_s"] + r["exact_s"] for r in done)
    print("exact plan, time per finished pair: median %.1f s, 90%% within %.0f s, max %.0f s (%.0f min)" % (
        statistics.median(tf), tf[int(0.9 * len(tf))], max(tf), max(tf) / 60))
    rp = [r["random"]["primes"] for r in rows]
    print("probabilistic plan: %d..%d primes (median %g), median time %.2f s" % (
        min(rp), max(rp), statistics.median(rp), statistics.median(r["random_s"] for r in rows)))
    plt = pyplot()
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(8, 3))
        x = [float(o[2]) for o in out]
        ax[0].scatter(x, [o[4] for o in out], s=8, facecolors="none", edgecolors="C1", label="exact plan, T")
        ax[0].scatter(x, [o[3] for o in out], s=8, facecolors="none", edgecolors="C0", label="probabilistic plan")
        ax[0].set(xscale="log", yscale="log", xlabel="n_cut", ylabel="primes", title="a")
        ax[1].scatter(x, [float(o[6]) for o in out], s=8, facecolors="none", edgecolors="C1", label="exact plan")
        ax[1].scatter(x, [o[5] for o in out], s=8, facecolors="none", edgecolors="C0", label="probabilistic plan")
        ax[1].set(xscale="log", yscale="log", xlabel="n_cut", ylabel="time, s", title="b")
        ax[0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, "figS1_plans.pdf"))


if __name__ == "__main__":
    main()
