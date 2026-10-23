"""Composition of the two training sets: genomes per phylum and labelled positives per category, side by side.

Left column: task 1 (358 genomes, 14 donor/acceptor categories). Right column: task 2 (3 378 genomes, 20 pathway
traits). Top row: genomes per phylum (top 12, others pooled). Bottom row: positives per category with the fraction
of positives inherited from FAPROTAX shaded for task 2 (taxonomy-inferred labels, a disclosure the paper needs).

usage: python src/85_dataset_figure.py
Output: docs/figures_paper/fig10_training_data.{png,pdf}
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
OUT = ROOT / "docs/figures_paper"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
                     "figure.dpi": 300, "savefig.bbox": "tight"})
BLUE, LIGHT, GREY = "#2f6f9f", "#a9c4dc", "#9aa0a6"


def phyla_panel(ax, counts, title):
    s = counts.sort_values(ascending=False)
    if len(s) > 12:
        s = pd.concat([s.iloc[:12], pd.Series({f"{len(s) - 12} other phyla": s.iloc[12:].sum()})])
    y = np.arange(len(s))[::-1]
    ax.barh(y, s.values, color=BLUE, height=0.65)
    for yi, v in zip(y, s.values):
        ax.annotate(f"{int(v)}", (v, yi), xytext=(3, 0), textcoords="offset points", va="center", fontsize=6.4)
    ax.set_yticks(y)
    ax.set_yticklabels(s.index, fontsize=6.8)
    ax.set_xlim(right=s.values.max() * 1.18)
    ax.set_xlabel("genomes")
    ax.set_title(title, fontsize=8, loc="left")


def category_panel(ax, df, name_col, title, shade_col=None):
    d = df.sort_values("positives", ascending=True)
    y = np.arange(len(d))
    ax.barh(y, d.positives, color=LIGHT if shade_col else BLUE, height=0.65, label="labelled positives")
    if shade_col:
        ax.barh(y, d.positives * d[shade_col], color=BLUE, height=0.65, label="of which inherited from FAPROTAX")
        ax.legend(frameon=False, fontsize=6.4, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1)
    for yi, (v, p) in enumerate(zip(d.positives, d.positive_rate)):
        ax.annotate(f"{int(v)} ({p:.0%})", (v, yi), xytext=(3, 0), textcoords="offset points", va="center", fontsize=6.2)
    ax.set_yticks(y)
    ax.set_yticklabels([str(x).replace("donor_", "D: ").replace("acceptor_", "A: ").replace("_", " ") for x in d[name_col]],
                       fontsize=6.6)
    ax.set_xscale("log")
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlim(right=d.positives.max() * 6)
    ax.set_xlabel("labelled positives (log; % = positive rate)")
    ax.set_title(title, fontsize=8, loc="left")


def main():
    from pathway_common import load_pathways
    from redox_common import load_labelled
    lab, _, _ = load_labelled()
    lab2, _, _ = load_pathways()
    t1 = pd.read_csv(ROOT / "results/dataset_task1_categories.csv")
    t2 = pd.read_csv(ROOT / "results/dataset_task2_traits.csv")
    fig, ax = plt.subplots(2, 2, figsize=(7.6, 8.2), gridspec_kw={"height_ratios": [1, 1.35], "hspace": 0.42, "wspace": 0.75})
    phyla_panel(ax[0, 0], lab.phylum.value_counts(), f"Task 1 training set — {len(lab)} genomes, {lab.genus.nunique()} genera, {lab.phylum.nunique()} phyla")
    phyla_panel(ax[0, 1], lab2.phylum.value_counts(), f"Task 2 training set — {len(lab2)} genomes, {lab2.genus.nunique()} genera, {lab2.phylum.nunique()} phyla")
    category_panel(ax[1, 0], t1, "category", "Task 1 — 14 donor/acceptor categories")
    category_panel(ax[1, 1], t2, "trait", "Task 2 — 20 pathway traits", shade_col="fraction_of_positives_from_FAPROTAX")
    fig.savefig(OUT / "fig10_training_data.png")
    fig.savefig(OUT / "fig10_training_data.pdf")
    print("wrote fig10_training_data")


if __name__ == "__main__":
    main()
