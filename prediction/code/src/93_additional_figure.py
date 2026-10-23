"""Fig. 4: fraction of 300 high-quality habitat MAGs predicted positive by the final model RedoxKO, per category (X4).
Categories with released yes/no calls (blue) first, then score-only categories (grey). For each expectation recorded
before the analysis, a solid outline marks the habitat expected to be higher and a dashed outline the habitat it was
compared with (results/additional/rf_habitat/x4_tests.csv).
MAG-isolate agreement by completeness is reported in Table 5. Writes docs/figures_paper/fig11_mags_habitats.{png,pdf}."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results/additional"
OUT = ROOT / "docs/figures_paper"
INK, MUTED = "#1f1f1f", "#6b6b6b"
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.edgecolor": MUTED, "xtick.color": INK,
                     "ytick.color": INK})

fr = pd.read_csv(R / "rf_habitat/x4_predicted_fraction_by_habitat.csv", index_col=0)
n = fr.loc["n_genomes"]
fr = fr.drop(index="n_genomes").astype(float)
released = {"acceptor_oxygen": "O₂ respiration", "acceptor_sulfate": "Sulfate reduction",
            "acceptor_CO2_methanogen_acetogen": "Methanogen/acetogen", "acceptor_sulfur": "Sulfur reduction",
            "donor_sulfur_oxidation": "Sulfur oxidation", "donor_ammonia_oxidation": "Ammonia oxidation"}
other = {"donor_hydrogen_oxidation": "H₂ oxidation", "donor_iron_oxidation": "Fe(II) oxidation",
         "acceptor_nitrate": "Nitrate reduction",
         "donor_CO_oxidation": "CO oxidation", "acceptor_iron": "Fe(III) reduction",
         "donor_nitrite_oxidation": "Nitrite oxidation"}
import json
assert set(released) == set(json.load(open(ROOT / "results/models/RF/metadata.json"))["reliable_traits"])
cols = ["hydrothermal", "acid_mine", "marine_surface", "digester", "human_gut"]
colnames = ["Hydrothermal\nvents", "Acid mine\ndrainage", "Seawater", "Anaerobic\ndigesters", "Human\ngut"]
rows = list(released) + list(other)
M = fr.loc[rows, cols].to_numpy()
tests = pd.read_csv(R / "rf_habitat/x4_tests.csv")
# categories tested on raw scores (Mann-Whitney): show median scores in their rows
SCORE_TESTED = sorted(set(tests[tests.kind.str.startswith("Mann")].trait))
W = ROOT / "data/additional/habitat"
sel = pd.read_csv(W / "selection.tsv", sep="\t").set_index("mag")
Pr = pd.concat([pd.read_csv(p) for p in sorted((W / "pred_rf").glob("batch_*.csv"))], ignore_index=True)
gcol = "genome" if "genome" in Pr else Pr.columns[0]
S = Pr.pivot_table(index=gcol, columns="trait", values="score")
S["habitat"] = sel.habitat.reindex(S.index)
for k in SCORE_TESTED:
    med = S.groupby("habitat")[k].median()
    M[rows.index(k)] = [med[c] for c in cols]
expected = {(t, h) for t, h in zip(tests.trait, tests.higher_in)}
compared = {(t, h) for t, h in zip(tests.trait, tests.lower_in) if t != "donor_iron_oxidation"}

blue = LinearSegmentedColormap.from_list("b", ["#f4f8fd", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
grey = LinearSegmentedColormap.from_list("g", ["#f7f7f7", "#d9d9d9", "#a6a6a6", "#6e6e6e", "#3d3d3d"])
fig, ax = plt.subplots(figsize=(4.8, 4.0))
nr = len(released)
for i, k in enumerate(rows):
    cm = blue if i < nr else grey
    for j, h in enumerate(cols):
        v = M[i, j]
        y = i + (0.5 if i >= nr else 0)
        ax.add_patch(Rectangle((j - 0.5, y - 0.5), 1, 1, fc=cm(v), ec="white", lw=1.2))
        ax.text(j, y, f"{v * 100:.0f}" + ("*" if k in SCORE_TESTED else ""), ha="center", va="center",
                fontsize=6.8, color="white" if v > 0.55 else INK)
        if (k, h) in expected:
            ax.add_patch(Rectangle((j - 0.46, y - 0.46), 0.92, 0.92, fc="none", ec="#1f1f1f", lw=1.6, zorder=5))
        elif (k, h) in compared:
            ax.add_patch(Rectangle((j - 0.42, y - 0.42), 0.84, 0.84, fc="none", ec="#1f1f1f", lw=1.0, ls=(0, (2, 1.5)),
                                   zorder=5))
ypos = [i + (0.5 if i >= nr else 0) for i in range(len(rows))]
ax.set_yticks(ypos, [released.get(k, other.get(k)) for k in rows], fontsize=7)
for lab in ax.get_yticklabels()[nr:]:
    lab.set_color(MUTED)
assert len(set(n[k] for k in cols)) == 1   # equal MAG count per habitat, given in the caption
ax.set_xticks(range(len(cols)), colnames, fontsize=7)
ax.xaxis.tick_top()
ax.set_xlim(-0.5, len(cols) - 0.5)
ax.set_ylim(len(rows), -0.5)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
ax.text(-0.02, -0.75 + 0, "", transform=ax.transData)
ax.text(len(cols) - 0.45, nr - 0.5 + 0.1, "", fontsize=6)
ax.text(1.02, (nr / 2 - 0.5), "yes/no\ncalls", transform=ax.get_yaxis_transform(), fontsize=6.5, color="#1c5cab",
        va="center")
ax.text(1.02, nr + 0.5 + (len(rows) - nr) / 2 - 0.5, "score\nonly", transform=ax.get_yaxis_transform(),
        fontsize=6.5, color=MUTED, va="center")

for ext in ("png", "pdf"):
    fig.savefig(OUT / f"fig11_mags_habitats.{ext}", dpi=300, bbox_inches="tight")
print("ok")
