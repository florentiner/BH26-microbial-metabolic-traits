"""Figure 1: workflow schematic (data -> features -> models -> evaluations), with an inset showing the core
evaluation idea (whole lineages held out). One headline number per evaluation; details in the tables."""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

OUT = Path(__file__).resolve().parents[1] / "docs/figures_paper"
INK, MUTED, LINE = "#1f1f1f", "#6b6b6b", "#b8b8b8"
T1, T2, BF = "#2a78d6", "#1baf7a", "#eb6834"
T1_BG, T2_BG, BF_BG, NEUTRAL_BG = "#eaf2fc", "#e6f6ef", "#fdeee6", "#f3f3f3"
plt.rcParams.update({"font.size": 7, "font.family": "DejaVu Sans"})

fig = plt.figure(figsize=(7.2, 5.25), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(-15.6, 62.6)
ax.axis("off")


def box(x, y, w, h, title, body, edge, bg, title_color=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.8", fc=bg, ec=edge, lw=0.8))
    ax.text(x + 1.0, y + h - 1.3, title, ha="left", va="top", fontsize=6.8, weight="bold", color=title_color)
    ax.text(x + 1.0, y + h - 4.3, body, ha="left", va="top", fontsize=6, color=INK, linespacing=1.35)


def arrow(x0, y0, x1, y1, color=LINE, ls="-"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7, lw=0.9, color=color,
                                 shrinkA=0, shrinkB=0, linestyle=ls))


# a  panel title and a numbered step banner (chevrons) over the four columns
PANEL_A = "a  Development of RedoxKO and EnergyKO-M"
ax.text(1, 61.2, PANEL_A, ha="left", va="center", fontsize=8, weight="bold", color=INK)
from matplotlib.patches import Polygon as _Poly
STEPS = [(1, 22.3, "1", "Data"), (23, 43.3, "2", "Features"), (44, 64.3, "3", "Models"), (65, 99, "4", "Evaluation")]
for k, (x0, x1, n, name) in enumerate(STEPS):
    y0, y1, p = 55.3, 58.5, 1.4
    pts = [(x0, y0), (x1 - p, y0), (x1, (y0 + y1) / 2), (x1 - p, y1), (x0, y1)]
    if k > 0:
        pts.append((x0 + p, (y0 + y1) / 2))
    ax.add_patch(_Poly(pts, closed=True, fc="#e9e9e9", ec="white", lw=1.2))
    ax.text(x0 + (2.6 if k else 1.4), (y0 + y1) / 2, n, ha="left", va="center", fontsize=6.4, color=MUTED, weight="bold")
    ax.text(x0 + (5.0 if k else 3.8), (y0 + y1) / 2, name, ha="left", va="center", fontsize=7, color=INK, weight="bold")

# a  data
box(1, 23, 20, 29, "Donor and acceptor\ngenomes", "\n358 genomes,\nlabeled for\n14 categories", T1, T1_BG, T1)
box(1, 6, 20, 12, "Pathway genomes", "3,378 genomes, labeled\nfor 20 traits", T2, T2_BG, T2)

# b  features (sequence on top so that no arrows cross)
box(24, 40, 18, 12, "Protein sequence", "ESM-C vectors →\nBacformer-large\ngenome vector", BF, BF_BG, BF)
box(24, 23, 18, 12, "KO panel", "852 energy-\nmetabolism KOs", T1, T1_BG, T1)
box(24, 6, 18, 12, "KO profile", "10,517 KOs + KEGG\nmodule completeness", T2, T2_BG, T2)

# c  models
box(45, 40, 18, 12, "RedoxKO-ens-BF", "two KO learners +\nBacformer-large", BF, BF_BG, BF)
box(45, 23, 18, 12, "RedoxKO", "final KO model", T1, T1_BG, T1)
box(45, 6, 18, 12, "EnergyKO-M", "final pathway\nmodel", T2, T2_BG, T2)

arrow(21, 46, 24, 46, T1)   # the same 358 genomes give the protein sequences for Bacformer
arrow(21, 29, 24, 29, T1)
arrow(21, 12, 24, 12, T2)
arrow(42, 46, 45, 46, BF)
arrow(42, 30, 45, 44, T1)
arrow(42, 29, 45, 29, T1)
arrow(42, 12, 45, 12, T2)

# 4  evaluations: each value in the colour of its model (blue RedoxKO, orange RedoxKO-ens-BF, green EnergyKO-M)
for x, c, name in [(66, T1, "RedoxKO"), (75.5, BF, "RedoxKO-ens-BF"), (89, T2, "EnergyKO-M")]:
    ax.text(x, 53.2, "●", color=c, fontsize=7, va="center")
    ax.text(x + 1.6, 53.2, name, color=c, fontsize=6, va="center", weight="bold")


def ebox(y, h, title, rows, edge=LINE):
    """rows: (label, value RedoxKO or pathway, value RedoxKO-ens-BF, colour of first value)"""
    ax.add_patch(FancyBboxPatch((66, y - h), 33, h, boxstyle="round,pad=0,rounding_size=0.6", fc="white", ec=edge,
                                lw=0.8))
    ax.text(68, y - 1.1, title, ha="left", va="top", fontsize=6, color=INK, weight="bold")
    for k, (lab, v1, v2, c1) in enumerate(rows):
        yy = y - 3.3 - 2.0 * k
        if v1 is None:                       # a grey note line
            ax.text(68, yy, lab, ha="left", va="top", fontsize=5.6, color=MUTED)
            continue
        ax.text(68, yy, lab, ha="left", va="top", fontsize=6, color=INK)
        ax.text(78, yy, v1, ha="left", va="top", fontsize=6, color=c1, weight="bold")
        if v2:
            ax.text(88.5, yy, v2, ha="left", va="top", fontsize=6, color=BF, weight="bold")


y = 51.0
for title, rows, h in [
        ("Genus- and family-level CV", [("AUPRC", "0.612 / 0.562", "0.625 / 0.585", T1)], 5.4),
        ("562 genomes of unseen genera", [("AUPRC", "0.555", "0.499", T1)], 5.4),
        ("100 newly described species", [("AUROC", "0.735", "0.746", T1), ("MCC", "0.420", "0.370", T1),
                                         ("other tools: MCC 0.519 METABOLIC, 0.194 FAPROTAX", None, None, None)], 9.0),
        ("321 MAGs, ≥ 70% complete", [("agreement", "≥ 96%", "≥ 96%", T1)], 5.4)]:
    ebox(y, h, title, rows)
    y -= h + 1.8
ebox(17.6, 7.4, "Pathway tests", [("genus-level CV", "AUPRC 0.501", "", T2), ("99 new species", "AUROC 0.673", "", T2)],
     edge=T2)
ax.plot([64.6, 64.6], [y + 1.8, 50.8], color=MUTED, lw=1.0)
arrow(63, 29, 64.6, 29, T1)
arrow(63, 44, 64.6, 46, BF)
arrow(63, 12, 66, 13.9, T2)

# e  using the released models on a new genome: inputs (grey pills) -> processes (white boxes) -> data (tinted)
#    -> outputs (orange documents), with the two alternative input routes and the two model branches drawn as forks
from matplotlib.patches import Polygon
ax.plot([1, 99], [2.6, 2.6], color=LINE, lw=0.6)
ax.text(1, 0.6, "b  Prediction for a new organism", ha="left", va="center", fontsize=8, weight="bold", color=INK)
TOP, BOT, H = -1.6, -9.0, 6.0   # top edges of the two rows and box height


def pbox(x, y, w, h, title, body, edge, fc="white", style="round,pad=0,rounding_size=0.8", tcol=None):
    ax.add_patch(FancyBboxPatch((x, y - h), w, h, boxstyle=style, fc=fc, ec=edge, lw=0.8))
    inset = 2.2 if "rounding_size=2.4" in style else 0.8   # keep text clear of the rounded ends of input pills
    ax.text(x + inset, y - 0.9, title, ha="left", va="top", fontsize=6.3, weight="bold", color=tcol or INK)
    if body:
        ax.text(x + inset, y - 3.0, body, ha="left", va="top", fontsize=5.5, color=INK, linespacing=1.25)


def doc(x, y, w, h, title, body, col, bg):
    f = 1.6
    ax.add_patch(Polygon([(x, y), (x + w - f, y), (x + w, y - f), (x + w, y - h), (x, y - h)], closed=True,
                         fc=bg, ec=col, lw=0.8))
    ax.plot([x + w - f, x + w - f, x + w], [y, y - f, y - f], color=col, lw=0.8)
    ax.text(x + 0.8, y - 0.9, title, ha="left", va="top", fontsize=6.3, weight="bold", color=col)
    ax.text(x + 0.8, y - 3.0, body, ha="left", va="top", fontsize=5.5, color=INK, linespacing=1.25)


# inputs (pills)
pbox(1, TOP, 12.5, H, "Accession", "NCBI\nassembly", MUTED, NEUTRAL_BG, "round,pad=0,rounding_size=2.4")
pbox(1, BOT, 12.5, H, "FASTA", "genome\nor MAG", MUTED, NEUTRAL_BG, "round,pad=0,rounding_size=2.4")
# protein routes
pbox(16, TOP, 13, H, "NCBI proteins", "from the NCBI\nannotation", MUTED)
pbox(16, BOT, 13, H, "Pyrodigal", "gene calling", MUTED)
arrow(13.5, TOP - H / 2, 16, TOP - H / 2, MUTED)
arrow(13.5, BOT - H / 2, 16, BOT - H / 2, MUTED)
# KO assignment (spans both rows) and KO profile (data)
yk, hk = TOP, TOP - BOT + H
pbox(32, yk, 15.5, hk, "KO assignment", "pyhmmer search of\nKOfam profiles with\nKofamScan thresholds", MUTED)
arrow(29, TOP - H / 2, 32, yk - hk / 2 + 1.5, MUTED)
arrow(29, BOT - H / 2, 32, yk - hk / 2 - 1.5, MUTED)
pbox(50.5, yk, 11.5, hk, "KO profile", "which KOs the\ngenome carries", MUTED, NEUTRAL_BG)
arrow(47.5, yk - hk / 2, 50.5, yk - hk / 2, MUTED)
# model branches
pbox(65, TOP, 13.5, H, "RedoxKO", "852-KO panel", T1, T1_BG, tcol=T1)
pbox(65, BOT, 13.5, H, "EnergyKO-M", "10,517 KOs and\nKEGG modules", T2, T2_BG, tcol=T2)
arrow(62, yk - hk / 2 + 1.5, 65, TOP - H / 2, MUTED)
arrow(62, yk - hk / 2 - 1.5, 65, BOT - H / 2, MUTED)   # arrows take the colour of the box they leave
# outputs (documents)
doc(81.5, TOP, 17.5, H, "Donor/acceptor table", "score, call,\nrelatives in training", T1, T1_BG)
doc(81.5, BOT, 17.5, H, "Pathway table", "score, call,\nmodule completeness", T2, T2_BG)
arrow(78.5, TOP - H / 2, 81.5, TOP - H / 2, T1)
arrow(78.5, BOT - H / 2, 81.5, BOT - H / 2, T2)

for ext in ("png", "pdf"):
    fig.savefig(OUT / f"fig13_workflow.{ext}", bbox_inches="tight", pad_inches=0.03)
print("ok")
