"""Figure 2: what the KO models rely on (donor and acceptor categories).
(a) Marker-gene enrichment among stable top features for the released RedoxKO (exact linear SHAP / TreeSHAP,
refitted per genus-blocked fold; results/additional/xai_redoxko) and for random forests on the same features
(results/trees_xai/da); (b) RedoxKO's stable features for four categories, mean |SHAP| relative to the category's top
feature; (c) random forest versus short rules (MCC). Supplementary figure: pathway enrichment (random forest).
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
X = ROOT / "results/trees_xai"
OUT = ROOT / "docs/figures_paper"
INK, MUTED, GRID = "#1f1f1f", "#6b6b6b", "#e6e6e6"
T1, T2, GREY = "#2a78d6", "#1baf7a", "#b5b5b5"
plt.rcParams.update({"font.size": 7, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False})

DA = {"donor_sulfur_oxidation": "Sulfur oxidation", "donor_ammonia_oxidation": "Ammonia oxidation",
      "donor_nitrite_oxidation": "Nitrite oxidation", "donor_hydrogen_oxidation": "H₂ oxidation",
      "donor_CO_oxidation": "CO oxidation", "donor_iron_oxidation": "Fe(II) oxidation", "donor_other": "Other donors",
      "acceptor_oxygen": "O₂ respiration", "acceptor_nitrate": "Nitrate reduction", "acceptor_sulfate": "Sulfate reduction",
      "acceptor_sulfur": "Sulfur reduction", "acceptor_CO2_methanogen_acetogen": "Methanogen/acetogen",
      "acceptor_iron": "Fe(III) reduction", "acceptor_other": "Other acceptors"}
PW = {"denitrification": "Denitrification", "dnra": "DNRA", "nitrate_reduction": "Nitrate reduction",
      "nitrogen_fixation": "Nitrogen fixation", "methanogenesis": "Methanogenesis", "acetogenesis": "Acetogenesis",
      "sulfate_reduction": "Sulfate reduction", "sulfur_reduction": "Sulfur reduction",
      "sulfur_oxidation": "Sulfur oxidation", "iron_reduction": "Fe(III) reduction", "iron_oxidation": "Fe(II) oxidation",
      "manganese_reduction": "Mn(IV) reduction", "arsenate_reduction": "Arsenate reduction",
      "fumarate_reduction": "Fumarate reduction", "hydrogen_oxidation": "H₂ oxidation", "fermentation": "Fermentation",
      "aerobic_chemoheterotrophy": "Aerobic chemoheterotrophy", "photoautotrophy": "Photoautotrophy",
      "photoheterotrophy": "Photoheterotrophy", "methylotrophy": "Methylotrophy"}
SYM = {"K17222": "soxA", "K17223": "soxX", "K17224": "soxB", "K17218": "sqr", "K27925": "tetH", "K00424": "cydX",
       "K04561": "norB", "K00373": "narJ", "K02575": "narK", "K02568": "napB", "K00374": "narI", "K02305": "norC",
       "K00362": "nirB", "K00376": "nosZ", "K11181": "dsrB", "K00124": "fdoH", "K27187": "dsrM", "K27188": "dsrK",
       "K00177": "korC", "K00127": "fdoI", "K00176": "korD", "K08352": "phsA", "K27190": "dsrO", "K00958": "sat",
       "K02257": "ctaB", "K03389": "hdrB2", "K00197": "cdhE", "K14127": "mvhD", "K03388": "hdrA2", "K16950": "asrA", "K14138": "acsB"}


SYM.update({"K02567": "napA", "K00371": "narH", "K00459": "ncd2", "K15635": "apgM", "K00937": "ppk1", "K00194": "cdhD",
            "K11261": "fwdE", "K00338": "nuoI", "K00239": "sdhA", "K05601": "hcp", "K17227": "soxZ", "K10945": "amoB",
            "K10946": "amoC", "K10535": "hao"})
RK = ROOT / "results/additional/xai_redoxko"
REDOX = "#2a78d6"


def enrich_panel(ax, rows_spec, legend=True):
    """rows_spec: list of (label, [(q, k, n, nk, color, filled, marker, dy)])"""
    y = np.arange(len(rows_spec))[::-1]
    for yi, (lab, pts) in zip(y, rows_spec):
        if all(p[3] == 0 for p in pts):
            ax.scatter(0.15, yi, s=14, marker="x", color=GREY, lw=0.8, zorder=3, clip_on=False)
            ax.text(0.5, yi, "no marker gene among features", va="center", fontsize=5.5, color=MUTED)
            continue
        for q, k, n, nk, col, filled, mk, dy, lbl in pts:
            v = -np.log10(max(q, 1e-12))
            sig = q < 0.05
            ax.plot([0, v], [yi + dy, yi + dy], color=col if sig else GREY, lw=0.9, zorder=2)
            ax.scatter(v, yi + dy, s=18, marker=mk, color=(col if sig else "white") if filled else "white",
                       edgecolor=col if sig else GREY, lw=0.9, zorder=3, clip_on=False)
            if lbl:
                ax.text(v + 0.25, yi + dy, f"{int(k)}/{int(n)}", va="center", fontsize=5.3, color=INK if sig else MUTED)
    ax.axvline(-np.log10(0.05), color=MUTED, lw=0.6, ls=(0, (3, 2)))
    ax.set_yticks(y, [r[0] for r in rows_spec], fontsize=6.3)
    ax.set_ylim(-0.7, len(rows_spec) - 0.3)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_xlabel("Marker-gene enrichment among\nstable top KOs, −log10 q (BH)")
    ax.text(-np.log10(0.05) + 0.15, len(rows_spec) - 0.45, "q = 0.05", fontsize=5.3, color=MUTED, va="bottom")


fig = plt.figure(figsize=(7.2, 6.2), dpi=300)
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.05], height_ratios=[1.6, 1.0], wspace=0.75, hspace=0.5)

# a  enrichment, Task 1, RedoxKO and random forest
a = fig.add_subplot(gs[:, 0])
er = pd.read_csv(RK / "enrichment.csv").set_index("trait")
ef = pd.read_csv(X / "da/xai_enrichment_rf.csv").set_index("trait")
spec = []
for k, lab in DA.items():
    r, f = er.loc[k], ef.loc[k]
    spec.append((lab, [(f.q_bh if f.n_known_markers_in_features else 1, f.n_stable_top_known, f.n_stable_top,
                        f.n_known_markers_in_features, REDOX, True, "o", 0.0, True)]))
enrich_panel(a, spec)
a.set_xlim(0, 11)
a.set_title("a  Marker genes among top features", loc="left", fontsize=8, color=INK, x=-0.55, pad=8)

# b  RedoxKO stable features, four categories
b = fig.add_subplot(gs[0, 1])
st = pd.read_csv(X / "da/xai_stable_rf.csv")   # final model RedoxKO (random forest)
show = [("acceptor_nitrate", "Nitrate reduction", 5), ("acceptor_sulfate", "Sulfate reduction", 5),
        ("acceptor_CO2_methanogen_acetogen", "Methanogen/acetogen", 5), ("acceptor_oxygen", "O₂ respiration", 5)]
yy, ys, seps = 0.0, [], []
for k, title, n in show:
    s = st[st.trait == k].sort_values("mean_abs_shap", ascending=False).head(n)
    top = s.mean_abs_shap.max()
    yy += 1.3
    seps.append((yy - 1.0, title))
    for _, r in s.iterrows():
        v = r.mean_abs_shap / top
        b.barh(-yy, v, height=0.72, color=REDOX if r.known_marker else "white", edgecolor=REDOX, lw=0.7)
        pos = r.shap_direction > 0
        b.text(v + 0.03, -yy, "+" if pos else "−", va="center", fontsize=6.5, color=INK if pos else "#c0392b")
        ys.append((-yy, SYM.get(r.KO, r.KO)))
        yy += 1
b.set_yticks([p for p, _ in ys], [l for _, l in ys], fontsize=5.8, style="italic")
for y0, title in seps:
    b.text(-0.02, -y0 + 0.25, title, fontsize=6.3, weight="bold", color=INK, transform=b.get_yaxis_transform(),
           ha="right", va="center")
b.set_ylim(-yy + 0.3, 0.2)
b.set_xlim(0, 1.18)
b.set_xticks([0, 0.5, 1.0])
b.set_xlabel("Mean |SHAP|, relative to top feature")
b.grid(axis="x", color=GRID, lw=0.5)
b.tick_params(axis="y", length=0)
from matplotlib.patches import Patch
b.legend(handles=[Patch(fc=REDOX, ec=REDOX, label="marker gene"),
                  Patch(fc="white", ec=REDOX, lw=0.7, label="not a marker gene")],
         loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=2, fontsize=5.5, frameon=False)
b.set_title("b  Most stable features of RedoxKO", loc="left", fontsize=8, color=INK, x=-0.62, pad=8)
b.text(1.0, 1.01, "+/− effect of gene presence", transform=b.transAxes, ha="right", fontsize=5.3, color=MUTED)

# c  model vs short rules
d = fig.add_subplot(gs[1, 1])
rv = pd.read_csv(X / "da/rules_vs_best.csv")
order = ["donor_ammonia_oxidation", "acceptor_CO2_methanogen_acetogen", "donor_sulfur_oxidation", "acceptor_nitrate",
         "acceptor_sulfate", "acceptor_oxygen", "donor_hydrogen_oxidation", "acceptor_sulfur"]
for i, k in enumerate(order[::-1]):
    s = rv[rv.trait == k]
    model = s.best_mcc.iloc[0]
    tree = s[s.rule_set.str.startswith("tree")].rule_mcc
    mark = s[s.rule_set.str.startswith("marker")].rule_mcc
    pts = [model] + list(tree) + list(mark)
    # 95% genus-bootstrap interval of the gap RedoxKO minus rule, drawn around the rule as model - interval
    for sel, dy, col in (("marker", -0.17, INK), ("tree", 0.17, GREY)):
        r = s[s.rule_set.str.startswith(sel)]
        if len(r):
            d.plot([model - r.delta_ci_high.iloc[0], model - r.delta_ci_low.iloc[0]], [i + dy, i + dy], color=col,
                   lw=0.7, zorder=2)
    if len(mark):
        d.scatter(mark.iloc[0], i - 0.17, marker="D", s=14, color="white", edgecolor=INK, lw=0.8, zorder=3)
    if len(tree):
        d.scatter(tree.iloc[0], i + 0.17, marker="s", s=14, color=GREY, zorder=3)
    d.scatter(model, i, s=18, color=REDOX, zorder=4)
d.set_yticks(range(len(order)), [DA[k] for k in order[::-1]], fontsize=6.0)
d.set_xlim(0, 1)
d.set_xlabel("MCC, genus-level CV (higher is better)")
d.grid(axis="x", color=GRID, lw=0.5)
d.scatter([], [], s=18, color=REDOX, label="RedoxKO")
d.scatter([], [], marker="s", s=14, color=GREY, label="tree, depth ≤ 3")
d.scatter([], [], marker="D", s=14, color="white", edgecolor=INK, lw=0.8, label="marker-gene rule")
d.legend(loc="lower center", fontsize=5.5, frameon=False, scatterpoints=1, bbox_to_anchor=(0.45, 1.0), ncol=3,
         handletextpad=0.2, columnspacing=0.8)
d.set_title("c  Learned model vs marker rules", loc="left", fontsize=8, color=INK, x=-0.62, pad=16)
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"fig14_xai.{ext}", bbox_inches="tight", pad_inches=0.03)

# Supplementary: pathway traits, random forest
fs, ax = plt.subplots(figsize=(3.6, 4.2), dpi=300)
ep = pd.read_csv(X / "pw/xai_enrichment_rf.csv").set_index("trait")
spec = [(lab, [(ep.loc[k].q_bh if ep.loc[k].n_known_markers_in_features else 1, ep.loc[k].n_stable_top_known,
                ep.loc[k].n_stable_top, ep.loc[k].n_known_markers_in_features, T2, True, "o", 0, True)])
        for k, lab in PW.items()]
enrich_panel(ax, spec)
ax.set_xlim(0, 12.5)
ax.set_title("Pathway traits (random forest)", loc="left", fontsize=8, color=INK)
for ext in ("png", "pdf"):
    fs.savefig(OUT / f"figS_xai_pathways.{ext}", bbox_inches="tight", pad_inches=0.03)
print("ok")
