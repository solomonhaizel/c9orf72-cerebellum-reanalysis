#!/usr/bin/env python3
# Statistical test for Figure 1b: total C9ORF72, control vs C9-ALS (GSE67196).
# Two-sided Mann-Whitney U in SciPy, with the method pinned to "exact" so a future
# SciPy release cannot change the result through automatic method selection.
# Reads data/total_c9_normalized.tsv, writes results/stats.tsv; the figure
# (plot_c9_panel.R) reads the p-value from there and does not recompute the test.
import csv, math, os, sys
import numpy as np
from scipy import stats

INFILE  = "data/total_c9_normalized.tsv"
OUTDIR  = "results"
GROUPS  = {"control", "c9ALS"}

def load(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    if not rows:
        sys.exit("ERROR: no rows in input.")
    acc = [r["sample"] for r in rows]
    if len(acc) != len(set(acc)):
        sys.exit("ERROR: sample accessions are not unique.")
    seen = set()
    data = {"control": [], "c9ALS": []}
    for r in rows:
        g = r["group"]
        if g not in GROUPS:
            sys.exit(f"ERROR: unexpected group '{g}' (allowed: {sorted(GROUPS)}).")
        try:
            x = float(r["norm_count"])
        except (ValueError, KeyError):
            sys.exit(f"ERROR: non-numeric norm_count for {r.get('sample')}.")
        if not math.isfinite(x) or x < 0:
            sys.exit(f"ERROR: norm_count for {r.get('sample')} is not finite/nonnegative ({x}).")
        data[g].append(x)
    return np.array(data["control"]), np.array(data["c9ALS"])

def main():
    ctrl, c9 = load(INFILE)
    # validation: exactly 8 per group, neither empty
    if len(ctrl) == 0 or len(c9) == 0:
        sys.exit("ERROR: a group is empty.")
    if len(ctrl) != 8 or len(c9) != 8:
        sys.exit(f"ERROR: expected 8 control and 8 C9-ALS; got {len(ctrl)} and {len(c9)}.")

    result = stats.mannwhitneyu(
        c9,
        ctrl,
        alternative="two-sided",
        method="exact",
    )
    p = result.pvalue
    fold = c9.mean() / ctrl.mean()

    os.makedirs(OUTDIR, exist_ok=True)   # ensure output directory exists
    with open(os.path.join(OUTDIR, "stats.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["comparison", "n_control", "n_c9als", "fold_change", "p_value", "test"])
        w.writerow(["C9-ALS vs control", len(ctrl), len(c9), round(fold, 4),
                    f"{p:.6g}", "two-sided Mann-Whitney U (exact)"])

    print(f"control n={len(ctrl)} mean={ctrl.mean():.1f}")
    print(f"c9ALS   n={len(c9)} mean={c9.mean():.1f}")
    print(f"fold={fold:.3f}  p={p:.4f}  (SciPy, exact) -> {OUTDIR}/stats.tsv")

if __name__ == "__main__":
    main()
