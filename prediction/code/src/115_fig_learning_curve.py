"""Fig. 5: learning curve of RedoxKO (macro AUPRC with genera held out against the number of labeled training genera),
with the log-linear fit and its slope interval, extrapolated only to twice the measured range.
Writes docs/figures_paper/fig15_learning_curve.{png,pdf}. (The label-noise simulation moved to the supplement.)"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/figures_paper"
INK, MUTED, GRID, BLUE = "#1f1f1f", "#6b6b6b", "#e6e6e6", "#2a78d6"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})

d = pd.read_csv(ROOT / "results/learning_curve_rf.csv")
d = d[(d.arm == "RF") & (d.split_level == "genus")]
# full size: the paper's folds and the five further random allocations (review round 3), not the paper's folds alone
al = pd.read_csv(ROOT / "results/allocations_rf_lx.csv").query("level == 'genus' and allocation != 'gtdb'")
full = pd.DataFrame({"arm": "RF", "split_level": "genus", "frac_genera": 1.0, "draw": range(len(al)),
                     "macro_auprc": al.macro_RF.to_numpy(), "mean_train_genera": 148.0})
d = pd.concat([d[d.frac_genera < 1], full], ignore_index=True)
g = d.groupby("frac_genera").agg(n=("mean_train_genera", "mean"), m=("macro_auprc", "mean"),
                                 s=("macro_auprc", "std")).reset_index()
x, y = np.log10(d.mean_train_genera), d.macro_auprc
b, a = np.polyfit(x, y, 1)
se = np.sqrt(((y - (a + b * x)) ** 2).sum() / (len(x) - 2) / ((x - x.mean()) ** 2).sum())
nmax = g.n.max()
fig, ax = plt.subplots(figsize=(4.2, 2.9))
xs_m = np.logspace(np.log10(g.n.min() * 0.9), np.log10(nmax), 50)
xs = xs_m
ax.fill_between(xs, a + b * np.log10(xs) - 1.96 * se * np.abs(np.log10(xs) - x.mean()),
                a + b * np.log10(xs) + 1.96 * se * np.abs(np.log10(xs) - x.mean()), color=BLUE, alpha=0.12, lw=0)
ax.plot(xs_m, a + b * np.log10(xs_m), color=BLUE, lw=1.6)
ax.errorbar(g.n, g.m, yerr=g.s.fillna(0), fmt="o", ms=5, color=BLUE, capsize=2, lw=1, zorder=3)
y2 = float("nan")
ax.set_xscale("log")
ax.set_xticks([37, 74, 111, 148], ["37", "74", "111", "148"])
ax.minorticks_off()
ax.set_xlabel("Labeled training genera (log scale)")
ax.set_ylabel("Macro AUPRC, genus-level CV (higher is better)")
ax.grid(color=GRID, lw=0.5)
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"fig15_learning_curve.{ext}", dpi=300, bbox_inches="tight")
print("ok", round(a, 3), round(b, 3), "slope CI", round(b - 1.96 * se, 3), round(b + 1.96 * se, 3))
print(g.round(3).to_string())
