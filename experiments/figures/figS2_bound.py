"""Fig. S2: the dimension bound of Theorem 1 against the least dimension of the segment automaton on the 1087
distinct ladders of the benchmark (Sec. 7.2). Data: data/runs/ladder_dims.csv (experiments/checks/ladder_dims.py);
bound = bound of Thm. 1, min_dim = least dimension (rank modulo three primes), wh = w + h.
Writes experiments/figures/out/figS2_pairs.csv (distinct pairs (min_dim, bound), number of ladders, marker size
class: 1 one ladder, 2 from 2 to 9, 3 ten or more; part a) and figS2_box.csv (ratio bound / min_dim by w + h:
count, min, quartiles (statistics.quantiles, inclusive), max; part b), and figS2_bound.pdf when matplotlib is
installed. Usage: python -I experiments/figures/figS2_bound.py"""
import collections, csv, os, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RUNS, OUT, pyplot      # noqa: E402


def main():
    rows = list(csv.DictReader(open(os.path.join(os.environ.get("RUNS", RUNS), "ladder_dims.csv"))))
    assert len(rows) == 1087, len(rows)
    os.makedirs(OUT, exist_ok=True)
    groups = collections.defaultdict(list)
    for r in rows:
        groups[int(r["w"]) + int(r["h"])].append(int(r["bound"]) / int(r["min_dim"]))
    pairs = collections.Counter((int(r["min_dim"]), int(r["bound"])) for r in rows)
    with open(os.path.join(OUT, "figS2_pairs.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["min_dim", "bound", "count", "size"])
        for (m, b), c in sorted(pairs.items()):
            w.writerow([m, b, c, 1 if c == 1 else (2 if c < 10 else 3)])
    with open(os.path.join(OUT, "figS2_box.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["wh", "count", "min", "q1", "median", "q3", "max"])
        for k in sorted(groups):
            g = sorted(groups[k])
            q = statistics.quantiles(g, n=4, method="inclusive") if len(g) > 1 else [g[0]] * 3
            w.writerow([k, len(g)] + ["%.4g" % v for v in (g[0], q[0], q[1], q[2], g[-1])])
    exact = sum(int(r["bound"]) == int(r["min_dim"]) for r in rows)
    print("ladders %d; bound equal to the least dimension on %d; bound >= least dimension on all: %s" % (
        len(rows), exact, all(int(r["bound"]) >= int(r["min_dim"]) for r in rows)))
    print("median ratio by w + h:", {k: round(statistics.median(v), 3) for k, v in sorted(groups.items())})
    plt = pyplot()
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(8, 3))
        ax[0].scatter([m for m, b in pairs], [b for m, b in pairs], s=[6 * (1 + (c > 1) + (c >= 10)) for c in pairs.values()],
                      facecolors="none", edgecolors="C0")
        lim = max(b for m, b in pairs)
        ax[0].plot([1, lim], [1, lim], "--", lw=0.8, color="gray")
        ax[0].set(xscale="log", yscale="log", xlabel="least dimension", ylabel="bound of Theorem 1", title="a")
        ks = sorted(groups)
        ax[1].boxplot([groups[k] for k in ks], positions=ks)
        ax[1].set(yscale="log", xlabel="w + h", ylabel="bound / least dimension", title="b")
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, "figS2_bound.pdf"))


if __name__ == "__main__":
    main()
