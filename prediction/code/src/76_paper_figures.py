"""Publication figures for the BH26 manuscript (docs/paper_tables_figures.md).

Every panel is drawn from a result file; nothing is recomputed. Output: docs/figures_paper/*.png (300 dpi) and
the same figures as .pdf for the manuscript.

usage: python src/76_paper_figures.py
"""
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
V = ROOT / "validation_new_organisms"
OUT = ROOT / "docs/figures_paper"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25,
                     "grid.linewidth": 0.5, "figure.dpi": 300, "savefig.bbox": "tight"})
# one categorical hue per model family, assigned in fixed order and never cycled
C = {"ML": "#2f6f9f", "DL": "#c1666b", "ML+DL": "#6b8f3a", "rule": "#8a7f8d", "null": "#9aa0a6", "LLM": "#d68c45"}


def save(fig, name):
    fig.savefig(OUT / f"{name}.png")
    fig.savefig(OUT / f"{name}.pdf")
    plt.close(fig)
    print("wrote", name)


def kind_of(arm):
    from importlib import import_module
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    return import_module("67_metrics_tables").ARMS.get(arm, ("ML", ""))[0]


def fig1_arms():
    """Macro AUPRC of every arm at both split levels, with cluster-bootstrap intervals."""
    da = pd.read_csv(RES / "metrics_macro.csv")
    tx = pd.read_csv(RES / "trees_xai/da/metrics_macro.csv")
    d = pd.concat([da, tx]).drop_duplicates(["arm", "split_level"], keep="first")
    g = d[d.split_level == "genus"].set_index("arm")
    f = d[d.split_level == "family"].set_index("arm")
    keep = [a for a in g.index if a in f.index]
    g, f = g.loc[keep].sort_values("macro_auprc"), f.loc[keep]
    fig, ax = plt.subplots(figsize=(5.2, 7.0))
    y = np.arange(len(g))
    for i, arm in enumerate(g.index):
        k = kind_of(arm)
        ax.plot([g.macro_auprc_ci_low[arm], g.macro_auprc_ci_high[arm]], [i, i], color=C.get(k, "#444"), lw=1.6,
                solid_capstyle="round", alpha=0.85)
        ax.plot(g.macro_auprc[arm], i, "o", ms=4.5, color=C.get(k, "#444"), zorder=3)
        ax.plot(f.macro_auprc[arm], i, "|", ms=7, mew=1.3, color=C.get(k, "#444"), alpha=0.75, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(g.index, fontsize=6.5)
    ax.set_xlabel("Macro AUPRC (dot = held-out genera with 95% CI, tick = held-out families)")
    ax.set_ylim(-1, len(g))
    handles = [plt.Line2D([], [], color=v, marker="o", ls="", label=k) for k, v in C.items() if k in
               {kind_of(a) for a in g.index}]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=7)
    save(fig, "fig1_arms_macro_auprc")


def fig2_per_trait():
    """Per-trait AUPRC against the number of positives, showing where the macro average is lost."""
    d = pd.read_csv(RES / "metrics_per_trait.csv")
    d = d[(d.arm == "B2XL") & (d.split_level == "genus")]
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.scatter(d.n_pos, d.auprc, s=26, color=C["ML"], zorder=3)
    for _, r in d.iterrows():
        ax.annotate(r.trait.replace("donor_", "D:").replace("acceptor_", "A:"), (r.n_pos, r.auprc),
                    textcoords="offset points", xytext=(4, 3), fontsize=5.6, color="#444")
    ax.plot(d.n_pos, d.prevalence, "_", ms=8, color="#9aa0a6", label="positive rate (chance)")
    ax.set_xscale("log")
    ax.set_xlabel("positives among 358 labelled genomes (log scale)")
    ax.set_ylabel("AUPRC, held-out genera")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    save(fig, "fig2_per_trait_auprc")


def fig3_learning_curve():
    """The learning curve and its extrapolation, with the fit's uncertainty."""
    d = pd.read_csv(RES / "learning_curve.csv")
    g = d.groupby("frac_genera").agg(n=("mean_train_genera", "mean"), m=("macro_auprc", "mean"),
                                     s=("macro_auprc", "std")).reset_index()
    x, y = np.log10(d.mean_train_genera), d.macro_auprc
    b, a = np.polyfit(x, y, 1)
    n = len(x)
    resid = y - (a + b * x)
    se = np.sqrt((resid ** 2).sum() / (n - 2) / ((x - x.mean()) ** 2).sum())
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    xs = np.logspace(np.log10(30), np.log10(1600), 100)
    ax.fill_between(xs, a + (b - 1.96 * se) * np.log10(xs), a + (b + 1.96 * se) * np.log10(xs),
                    color=C["ML"], alpha=0.13, lw=0)
    ax.plot(xs, a + b * np.log10(xs), color=C["ML"], lw=1.6)
    ax.errorbar(g.n, g.m, yerr=g.s.fillna(0), fmt="o", ms=5, color=C["ML"], capsize=2, lw=1, zorder=3)
    ax.axvspan(30, g.n.max(), color="#000", alpha=0.04, lw=0)
    lo, hi = ax.get_ylim()
    ax.annotate("measured", (55, lo + 0.04 * (hi - lo)), fontsize=7, color="#555")
    ax.annotate("extrapolated", (330, lo + 0.04 * (hi - lo)), fontsize=7, color="#555")
    ax.set_xscale("log")
    ax.set_xlabel("training genera (log scale)")
    ax.set_ylabel("macro AUPRC, held-out genera")
    ax.set_title(f"AUPRC = {a:.3f} + {b:.3f}·log10(genera),  slope 95% CI ±{1.96 * se:.3f}", fontsize=7.5)
    save(fig, "fig3_learning_curve")


def fig4_label_noise():
    """Hiding positives costs little in training and a great deal in measurement."""
    d = pd.read_csv(RES / "trees_xai/da/headroom_label_noise.csv")
    g = d.groupby("rate")[["auprc_clean_vs_true", "auprc_noisytrain_vs_true", "auprc_clean_vs_noisytest"]].mean()
    fig, ax = plt.subplots(figsize=(4.6, 3.1))
    r = np.r_[0, g.index.to_numpy()]
    base = g.auprc_clean_vs_true.iloc[0]
    ax.plot(r * 100, np.r_[base, g.auprc_noisytrain_vs_true], "o-", ms=4, color=C["ML"],
            label="trained on damaged labels, scored against truth")
    ax.plot(r * 100, np.r_[base, g.auprc_clean_vs_noisytest], "s-", ms=4, color=C["DL"],
            label="clean model, scored against damaged labels")
    ax.set_xlabel("% of positives hidden")
    ax.set_ylabel("macro AUPRC")
    ax.legend(frameon=False, fontsize=6.6, loc="lower left")
    save(fig, "fig4_label_noise")


def fig5_validation():
    """Prospective validation: accuracy per category with counts, both models."""
    fs = {a: V / f"scored178_{a}_per_category.csv" for a in ("B2", "B2XLE_large")}
    if not all(f.exists() for f in fs.values()):
        return
    a = pd.read_csv(fs["B2"]).set_index("category")
    b = pd.read_csv(fs["B2XLE_large"]).set_index("category")
    keep = a[a.n >= 5].sort_values("n").index
    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    y = np.arange(len(keep))
    ax.barh(y - 0.2, a.loc[keep, "accuracy"], height=0.38, color=C["ML"], label="B2")
    ax.barh(y + 0.2, b.loc[keep, "accuracy"], height=0.38, color=C["ML+DL"], label="B2XL + Bacformer-large")
    for i, c in enumerate(keep):
        x = max(a.loc[c, "accuracy"], b.loc[c, "accuracy"])
        ax.annotate(f"n={int(a.loc[c, 'n'])}", (x + 0.015, i), va="center", fontsize=6, color="#555")
    ax.set_yticks(y)
    ax.set_yticklabels([c.replace("donor_", "D:").replace("acceptor_", "A:") for c in keep], fontsize=6.8)
    ax.set_xlabel("fraction of published growth tests called correctly")
    ax.set_xlim(0, 1)
    ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    save(fig, "fig5_prospective_validation")


def fig6_llm():
    """LLM versus trained models on both tasks, gene evidence only and with the organism name."""
    da = pd.read_csv(RES / "baseline2_llm_summary.csv")
    rows = []
    for arm, lab in [("H_genome", "Claude, evidence only"), ("H_named", "Claude, + organism name"),
                     ("H_gpt_genome", "GPT, evidence only"), ("H_gpt_named", "GPT, + organism name")]:
        s = da[(da.llm_arm == arm) & (da.model == arm)]
        if len(s):
            rows.append(("donors/acceptors", lab, float(s.macro_auprc.iloc[0]), "LLM"))
    for arm, val in [("B2", None), ("B2XLE_large", None), ("null_taxNN", None)]:
        s = da[(da.llm_arm == "H_genome") & (da.model == arm)]
        if len(s):
            rows.append(("donors/acceptors", arm, float(s.macro_auprc.iloc[0]),
                         "null" if arm.startswith("null") else "ML+DL" if "E_large" in arm else "ML"))
    pw = RES / "pathway_llm_summary.csv"
    if pw.exists():
        p = pd.read_csv(pw).set_index("arm").macro_auprc
        name = {"P_H_genome": "Claude, evidence only", "P_H_named": "Claude, + organism name"}
        for arm, v in p.items():
            rows.append(("pathways", name.get(arm, arm), float(v),
                         "LLM" if arm.startswith("P_H") else "null" if arm.startswith("null") else
                         "rule" if arm == "PA" else "ML"))
    d = pd.DataFrame(rows, columns=["task", "arm", "auprc", "kind"])
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.4), sharex=False)
    for ax, task in zip(axes, ["donors/acceptors", "pathways"]):
        s = d[d.task == task].sort_values("auprc")
        ax.barh(np.arange(len(s)), s.auprc, color=[C.get(k, "#444") for k in s.kind], height=0.62)
        ax.set_yticks(np.arange(len(s)))
        ax.set_yticklabels(s.arm, fontsize=6.6)
        ax.set_title(task, fontsize=8)
        ax.set_xlabel("macro AUPRC on the 96-genome subset")
    save(fig, "fig6_llm_comparison")


def fig7_funnel():
    """How the prospective validation set was built."""
    steps = [("PubMed species descriptions\nreporting redox chemistry", 553),
             ("binomials parsed", 570), ("genus/species/genome absent\nfrom all project data", 194),
             ("annotated RefSeq genome", 178), ("≥1 published growth test", 167), ("scored", 162)]
    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    y = np.arange(len(steps))[::-1]
    ax.barh(y, [n for _, n in steps], color=C["ML"], height=0.6)
    for i, (lab, n) in zip(y, steps):
        ax.annotate(f"{n}", (n + 8, i), va="center", fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels([s for s, _ in steps], fontsize=6.4)
    ax.set_xlabel("organisms / papers")
    ax.set_xlim(0, 640)
    save(fig, "fig7_validation_funnel")


if __name__ == "__main__":
    for fn in (fig1_arms, fig2_per_trait, fig3_learning_curve, fig4_label_noise, fig5_validation, fig6_llm,
               fig7_funnel):
        try:
            fn()
        except Exception as e:                       # a missing result file must not stop the rest
            print(f"skipped {fn.__name__}: {type(e).__name__}: {e}")
