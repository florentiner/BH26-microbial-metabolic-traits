"""Supplementary Fig. S2: effect of randomly hiding a fraction of training or test positives (closed-world labels).
Macro AUPRC (mean over categories of the mean over folds) per repeat; points are means over repeats and bars the
range across repeats. Input: results/trees_xai/da/headroom_label_noise.csv. Writes docs/figures_paper/figS2_label_noise.{png,pdf}."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/figures_paper"
INK, MUTED, GRID, BLUE = "#1f1f1f", "#6b6b6b", "#e6e6e6", "#2a78d6"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})

d = pd.read_csv(ROOT / "results/trees_xai/da/headroom_label_noise.csv")
cols = ["auprc_clean_vs_true", "auprc_noisytrain_vs_true", "auprc_clean_vs_noisytest"]
per = d.groupby(["rate", "seed", "trait"])[cols].mean().groupby(["rate", "seed"]).mean()  # macro per repeat
g = per.groupby("rate").agg(["mean", "min", "max"])
base = per["auprc_clean_vs_true"].mean()
fig, ax = plt.subplots(figsize=(4.2, 2.9))
r = np.r_[0, g.index.to_numpy()] * 100
# one model, one colour: the two conditions differ by line style and marker fill
for col, mk, ls, face, label in (("auprc_noisytrain_vs_true", "o", "-", BLUE, "Training labels"),
                                 ("auprc_clean_vs_noisytest", "s", "--", "white", "Test labels")):
    m = np.r_[base, g[(col, "mean")]]
    lo, hi = np.r_[base, g[(col, "min")]], np.r_[base, g[(col, "max")]]
    ax.errorbar(r, m, yerr=[m - lo, hi - m], fmt=mk, ls=ls, color=BLUE, mfc=face, mec=BLUE, ms=4.5, lw=1.6,
                capsize=2, label=label)
ax.legend(title="Positives hidden in", loc="lower left", fontsize=6.6, title_fontsize=6.6, frameon=False)
ax.set_xticks(r, [f"{int(x)}" for x in r])
ax.set_xlabel("Positives hidden at random (%)")
ax.set_ylabel("Macro AUPRC (higher is better)")
ax.grid(axis="y", color=GRID, lw=0.5)
for e in ("png", "pdf"):
    fig.savefig(OUT / f"figS2_label_noise.{e}", dpi=300, bbox_inches="tight")
print(g.round(3).to_string())
