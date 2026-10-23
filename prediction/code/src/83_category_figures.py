"""Per-category figures for both tasks: how much data each category has, and how well the best model does on it.

One figure per task, two panels: left = labelled positives per category (log scale), right = the best model's
per-category AUPRC on held-out genera with its bootstrap interval, the taxonomic null as a tick, and the positive
rate (chance level) as a grey dash. Categories are ordered by AUPRC so the strong and weak groups read at a glance.

usage: python src/83_category_figures.py
Output: docs/figures_paper/fig8_task1_categories.{png,pdf}, fig9_task2_categories.{png,pdf}
"""
from pathlib import Path

import matplotlib
import matplotlib.ticker
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = ROOT / "docs/figures_paper"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
                     "figure.dpi": 300, "savefig.bbox": "tight"})
BLUE, GREY, RED = "#2f6f9f", "#9aa0a6", "#c1666b"


def pretty(t):
    return (t.replace("donor_", "donor: ").replace("acceptor_", "acceptor: ").replace("_", " ")
             .replace("CO2 methanogen acetogen", "CO2 (methanogen/acetogen)"))


def draw(per_trait, arm, title, name, n_total):
    d = per_trait[(per_trait.arm == arm) & (per_trait.split_level == "genus")].copy()
    d = d.sort_values("auprc")
    y = np.arange(len(d))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.4, 0.32 * len(d) + 1.4), sharey=True,
                                 gridspec_kw={"width_ratios": [1, 1.6], "wspace": 0.08})
    a1.barh(y, d.n_pos, color=BLUE, height=0.62)
    for i, (n, p) in enumerate(zip(d.n_pos, d.prevalence)):
        a1.annotate(f"{int(n)}  ({p:.0%})", (n, i), xytext=(3, 0), textcoords="offset points", va="center",
                    fontsize=6.4, color="#444")
    a1.set_xscale("log")
    a1.set_xlim(right=d.n_pos.max() * 4)
    a1.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a1.set_xlabel(f"labelled positives (of {n_total}; % = positive rate)")
    a1.set_yticks(y)
    a1.set_yticklabels([pretty(t) for t in d.trait], fontsize=7)
    a1.set_title("data per category", fontsize=8, loc="left")
    for i, r in enumerate(d.itertuples()):
        a2.plot([r.auprc_ci_low, r.auprc_ci_high], [i, i], color=BLUE, lw=1.6, solid_capstyle="round")
        a2.plot(r.auprc, i, "o", ms=4.5, color=BLUE, zorder=3)
        a2.plot(r.null_taxNN_auprc, i, "|", ms=8, mew=1.4, color=RED, zorder=3)
        a2.plot(r.prevalence, i, "_", ms=8, mew=1.4, color=GREY, zorder=2)
    a2.set_xlim(0, 1)
    a2.set_xlabel(f"AUPRC on held-out genera")
    a2.set_title(f"{title}: best model per category", fontsize=8, loc="left")
    handles = [plt.Line2D([], [], color=BLUE, marker="o", ls="", label=f"{arm}"),
               plt.Line2D([], [], color=RED, marker="|", ls="", mew=1.4, label="taxonomic null"),
               plt.Line2D([], [], color=GREY, marker="_", ls="", mew=1.4, label="positive rate (chance)")]
    a2.legend(handles=handles, loc="lower right", frameon=False, fontsize=6.6)
    fig.savefig(OUT / f"{name}.png")
    fig.savefig(OUT / f"{name}.pdf")
    plt.close(fig)
    print("wrote", name, f"({len(d)} categories)")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    da = pd.read_csv(RES / "metrics_per_trait.csv")
    draw(da, "B2XLE_large", "Task 1 — electron donors and acceptors", "fig8_task1_categories", 358)
    pw = pd.concat([pd.read_csv(RES / "pathways/metrics_per_trait.csv"),
                    pd.read_csv(RES / "trees_xai/pw/metrics_per_trait.csv")]).drop_duplicates(["arm", "split_level", "trait"])
    best = "PB3" if (pw.arm == "PB3").any() else "PB2"
    draw(pw, best, "Task 2 — energy-metabolism pathways", "fig9_task2_categories", 3378)


if __name__ == "__main__":
    main()
