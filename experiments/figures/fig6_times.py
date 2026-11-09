"""Fig. 7: time of the method and of the QCEC enumeration of all n <= n_cut against n_cut on the 444 equal pairs
answered "equal" (Sec. 7.4; same pairs as Table 6). An enumeration not finished within 60 s is written as 60 s
(qcec_done = 0). ncut_j: n_cut spread by 10^u, |u| <= 0.04, for drawing only.
Writes experiments/figures/out/fig6_times.csv (the data of the paper figure) and fig6_times.pdf (a plain plot,
when matplotlib is installed); prints the numbers of the caption.
Usage: python -I experiments/figures/fig6_times.py"""
import csv, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RUNS, BENCH, OUT, jitter, pyplot      # noqa: E402


def main():
    runs = os.environ.get("RUNS", RUNS)
    ours = {r["id"]: r["ours"] for r in map(json.loads, open(os.path.join(runs, "method_main.jsonl")))}
    pairs = {p["id"]: p for p in map(json.loads, open(os.path.join(BENCH, "pairs.jsonl")))}
    q = {r["id"]: r for r in map(json.loads, open(os.path.join(runs, "qcec_ncut.jsonl")))}
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for i in sorted(q):
        done = q[i]["verdict"] == "equal"
        nc = ours[i]["n_cut"][pairs[i]["label"]]
        rows.append([i, int(q[i]["bits"]), nc, "%.6g" % (nc * 10 ** jitter(i)), ours[i]["time_s"],
                     q[i]["wall_s"] if done else 60.0, int(done)])
    with open(os.path.join(OUT, "fig6_times.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "bits", "ncut", "ncut_j", "method_s", "qcec_s", "qcec_done"])
        w.writerows(rows)
    nb, wb = [r for r in rows if not r[1]], [r for r in rows if r[1]]
    print("pairs %d: without bits %d, with bits %d; enumeration not finished: %d and %d" % (
        len(rows), len(nb), len(wb), sum(not r[6] for r in nb), sum(not r[6] for r in wb)))
    print("n_cut > 10^5: %d pairs, largest %d; longest method time %.3f s" % (
        sum(r[2] > 10 ** 5 for r in rows), max(r[2] for r in rows), max(r[4] for r in rows)))
    plt = pyplot()
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(8, 3), sharey=True, gridspec_kw={"width_ratios": [3, 2]})
        for a, part, title in ((ax[0], nb, "a  no bits a"), (ax[1], wb, "b  with bits a")):
            a.scatter([float(r[3]) for r in part], [r[4] for r in part], s=8, facecolors="none", edgecolors="C0",
                      label="method")
            fin = [r for r in part if r[6]]
            a.scatter([float(r[3]) for r in fin], [r[5] for r in fin], s=8, marker="^", facecolors="none",
                      edgecolors="C1", label="QCEC, finished")
            a.scatter([r[2] for r in part if not r[6]], [60] * sum(not r[6] for r in part), s=12, marker="|",
                      color="C1", label="QCEC, not finished in 60 s")
            a.axhline(60, ls="--", lw=0.8, color="gray")
            a.set(xscale="log", yscale="log", xlabel="n_cut", title=title)
        ax[0].set_ylabel("time, s")
        ax[0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, "fig6_times.pdf"))


if __name__ == "__main__":
    main()
