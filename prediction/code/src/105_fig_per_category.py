"""Fig. 3: RedoxKO and METABOLIC per category on two tests, as balanced accuracy ((sensitivity + specificity) / 2),
which scores a tool that calls every genome positive at 0.5. (a) 255 growth tests of the 100 primary new species;
(b) 562 external genomes of unseen genera (4,962 genome-category pairs). 95% percentile bootstrap intervals
(2,000 resamples within each category; papers for new species, genomes for external), seed 20261017.
Writes docs/figures_paper/fig12_new_species_per_category.{png,pdf}."""
import sys
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from metabolic_calls import calls_from_rows  # noqa: E402

OUT = ROOT / "docs/figures_paper"
INK, MUTED, GRID, BLUE, RULE = "#1f1f1f", "#6b6b6b", "#e6e6e6", "#2a78d6", "#4a4a4a"
SEED, B = 20261017, 2000
NAMES = {"acceptor_sulfate": "Sulfate reduction", "acceptor_CO2_methanogen_acetogen": "Methanogen/acetogen",
         "acceptor_oxygen": "O₂ respiration", "acceptor_sulfur": "Sulfur reduction",
         "donor_sulfur_oxidation": "Sulfur oxidation", "acceptor_nitrate": "Nitrate reduction",
         "donor_hydrogen_oxidation": "H₂ oxidation", "donor_iron_oxidation": "Fe(II) oxidation",
         "acceptor_iron": "Fe(III) reduction"}
import json as _json
RELEASED = set(_json.loads((ROOT / "results/models/RF/metadata.json").read_text())["reliable_traits"])  # final model RedoxKO
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": INK})


def bacc(y, c):
    return 0.5 * ((c[y == 1] == 1).mean() + (c[y == 0] == 0).mean())


def stats(d, cluster):
    rows = []
    rng = np.random.default_rng(SEED)
    for t, g in d.groupby("category"):
        if t not in NAMES or g.y.nunique() < 2:
            continue
        r = {"category": t, "n": len(g)}
        keys = g[cluster].unique()
        groups = {k: x for k, x in g.groupby(cluster)}
        bs = {"R": [], "M": []}
        for _ in range(B):
            s = pd.concat([groups[k] for k in rng.choice(keys, len(keys))])
            if s.y.nunique() < 2:
                continue
            bs["R"].append(bacc(s.y.to_numpy(), s.c_R.to_numpy()))
            bs["M"].append(bacc(s.y.to_numpy(), s.c_M.to_numpy()))
        for k, col in (("R", "c_R"), ("M", "c_M")):
            r[k] = bacc(g.y.to_numpy(), g[col].to_numpy())
            r[k + "_lo"], r[k + "_hi"] = np.percentile(bs[k], [2.5, 97.5])
        rows.append(r)
    return pd.DataFrame(rows).set_index("category")


x1 = pd.read_csv(ROOT / "results/additional/rf_paper/x1_calls_rf.csv")
new = pd.DataFrame({"category": x1.category, "y": x1.expected.astype(int), "c_R": x1.call_RF.astype(int),
                    "c_M": x1.call_METABOLIC.astype(int), "cl": x1.paper})
thr = {t: b["threshold"] for t, b in joblib.load(ROOT / "results/models/RF/model.joblib").items()}
p = pd.read_csv(ROOT / "results/external_validation_missing_rows_predictions.csv").query("arm == 'Random forest'").merge(
    pd.read_csv(ROOT / "results/external_validation_predictions.csv").query("arm == 'B2'")[["genome_id", "trait", "y"]],
    on=["genome_id", "trait"])
Me = calls_from_rows(pd.read_parquet(ROOT / "data/additional/metabolic_hmm/external.parquet"))
p = p[p.genome_id.isin(Me.index)]
ext = pd.DataFrame({"category": p.trait.to_numpy(), "y": p.y.astype(int).to_numpy(),
                    "c_R": (p.score.to_numpy() >= p.trait.map(thr).to_numpy()).astype(int),
                    "c_M": [int(Me.at[g, c]) for g, c in zip(p.genome_id, p.trait)], "cl": p.genome_id.to_numpy()})
sn, se = stats(new, "cl"), stats(ext, "cl")
order = se.sort_values("R").index.tolist()  # shared row order, by RedoxKO on the external genomes
sn.round(3).to_csv(OUT / "fig12_values_new_species.csv")
se.round(3).to_csv(OUT / "fig12_values_external.csv")

fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), sharey=True, gridspec_kw={"wspace": 0.2})
for ax, s, title in ((axes[0], sn, "a  100 newly described species (255 tests)"),
                     (axes[1], se, "b  562 external genomes of unseen genera")):
    y = np.arange(len(order))
    for k, color, mk, dy, face, lab in (("R", BLUE, "o", 0.15, BLUE, "RedoxKO"),
                                         ("M", RULE, "D", -0.15, "white", "METABOLIC")):
        for yi, t in zip(y, order):
            if t not in s.index:
                continue
            r = s.loc[t]
            a = 1.0 if r.n >= 10 else 0.35
            ax.errorbar(r[k], yi + dy, xerr=[[r[k] - r[k + "_lo"]], [r[k + "_hi"] - r[k]]], fmt=mk, color=color,
                        mfc=face, ms=4.2, capsize=1.8, lw=0.9, alpha=a, label=lab if yi == 0 else None)
    # number of tests per category of the new species, at the right edge of panel a
    for yi, t in zip(y, order):
        if ax is axes[0] and t in s.index:
            ax.text(1.015, yi, f"{int(s.loc[t].n)}", transform=ax.get_yaxis_transform(), fontsize=6, color=MUTED,
                    va="center", ha="left", clip_on=False)
    if ax is axes[0]:
        ax.text(1.015, len(order) - 0.45, "n", transform=ax.get_yaxis_transform(), fontsize=6, color=MUTED, va="center",
                ha="left", style="italic")
    ax.axvline(0.5, color=MUTED, lw=0.6, ls=(0, (3, 2)))
    ax.set_xlim(0.0, 1.02)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0], ["0", "0.25", "0.5", "0.75", "1"])
    ax.set_ylim(-0.6, len(order) - 0.2)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_title(title, loc="left", fontsize=8, color=INK)
    ax.set_xlabel("Balanced accuracy (higher is better)")
axes[0].set_yticks(np.arange(len(order)), [f"{NAMES[t]}{' †' if t in RELEASED else ''}" for t in order], fontsize=7)
for ax in axes:
    ax.text(0.51, len(order) - 0.45, "chance", fontsize=6, color=MUTED, va="center")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, ["RedoxKO", "METABOLIC"], loc="lower center", bbox_to_anchor=(0.55, 0.96), ncol=2,
           fontsize=7, frameon=False)
for ext_ in ("png", "pdf"):
    fig.savefig(OUT / f"fig12_new_species_per_category.{ext_}", dpi=300, bbox_inches="tight")
print(sn.round(2).to_string())
print(se.round(2).to_string())
