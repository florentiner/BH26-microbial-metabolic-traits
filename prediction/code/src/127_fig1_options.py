"""Ten alternative designs for the top part (a-d) of Figure 1, for the author to choose from. Panel e is unchanged.
Numbers are the paper's headline values (Tables 1, 3, 5, 6 and the Results text).
Output: docs/figures_paper/fig1_options/option_XX.png and contact_sheet.png
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle

OUT = Path(__file__).resolve().parents[1] / "docs/figures_paper/fig1_options"
OUT.mkdir(parents=True, exist_ok=True)
INK, MUTED, LINE = "#1f1f1f", "#6b6b6b", "#b8b8b8"
T1, T2, BF = "#2a78d6", "#1baf7a", "#eb6834"
T1_BG, T2_BG, BF_BG, NBG = "#eaf2fc", "#e6f6ef", "#fdeee6", "#f3f3f3"
TEST = "#7A5195"
plt.rcParams.update({"font.size": 7, "font.family": "DejaVu Sans"})


def canvas(h=3.4):
    fig = plt.figure(figsize=(7.2, h), dpi=200)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100 * h / 7.2)
    ax.axis("off")
    return fig, ax


def box(ax, x, y, w, h, title, body="", edge=MUTED, bg="white", tcol=INK, fs=6.6, r=0.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=bg, ec=edge, lw=0.8))
    ax.text(x + 1.0, y + h - 1.2, title, ha="left", va="top", fontsize=fs, weight="bold", color=tcol)
    if body:
        ax.text(x + 1.0, y + h - 3.9, body, ha="left", va="top", fontsize=fs - 0.8, color=INK, linespacing=1.3)


def arrow(ax, x0, y0, x1, y1, color=LINE, lw=0.9, ls="-"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7, lw=lw, color=color,
                                 shrinkA=0, shrinkB=0, linestyle=ls))


def label(ax, x, y, s, fs=8):
    ax.text(x, y, s, ha="left", va="center", fontsize=fs, weight="bold", color=INK)


def save(fig, n, title):
    fig.savefig(OUT / f"option_{n:02d}.png", bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    TITLES[n] = title


TITLES = {}
EVAL_DA = [("Genus-level CV", "AUPRC 0.576"), ("Family-level CV", "AUPRC 0.560"), ("562 external genomes", "AUPRC 0.511"),
           ("100 new species", "AUROC 0.73-0.75"), ("372 MAG-isolate pairs", "≥ 94% agree"),
           ("7 tools and baselines", "beat taxonomy tools")]
EVAL_PW = [("Genus-level CV", "AUPRC 0.501"), ("99 new species", "AUROC 0.67")]


# ---------------------------------------------------------------------------------------------- 1 grid (current)
def opt1():
    fig, ax = canvas(3.4)
    for x, s in [(1, "a  Data"), (24, "b  Features"), (45, "c  Models"), (66, "d  Evaluation")]:
        label(ax, x, 45.5, s)
        box(ax, 1, 17, 20, 10, "Donor and acceptor labels", "358 genomes, 14 categories", T1, T1_BG, T1)
    box(ax, 1, 3, 20, 10, "Pathway labels", "3,378 genomes, 20 traits", T2, T2_BG, T2)
    box(ax, 24, 31, 18, 10, "Protein sequence", "ESM-C -> Bacformer-large", BF, BF_BG, BF)
    box(ax, 24, 17, 18, 10, "KO panel", "852 KOs", T1, T1_BG, T1)
    box(ax, 24, 3, 18, 10, "KO profile", "10,517 KOs + modules", T2, T2_BG, T2)
    box(ax, 45, 31, 18, 10, "RedoxKO-ens-BF", "KO learners + Bacformer", BF, BF_BG, BF)
    box(ax, 45, 17, 18, 10, "RedoxKO", "final KO model", T1, T1_BG, T1)
    box(ax, 45, 3, 18, 10, "EnergyKO-M", "final pathway model", T2, T2_BG, T2)
    for y, c in [(36, BF), (22, T1), (8, T2)]:
        arrow(ax, 21, y, 24, y, c); arrow(ax, 42, y, 45, y, c)
    y = 41
    for t, m in EVAL_DA:
        box(ax, 66, y - 4.6, 33, 4.6, t, "", T1, "white", INK, 6)
        ax.text(98, y - 2.4, m, ha="right", va="center", fontsize=5.6, color=INK)
        y -= 5.4
    for t, m in EVAL_PW[:1]:
        box(ax, 66, 3, 33, 4.6, "Pathways: CV / 99 new species", "", T2, "white", INK, 6)
        ax.text(98, 4.4, "AUPRC 0.501 / AUROC 0.67", ha="right", va="center", fontsize=5.6)
    save(fig, 1, "Grid (current layout, tidied)")


# ---------------------------------------------------------------------------------------------- 2 swimlanes
def opt2():
    fig, ax = canvas(3.3)
    for y0, name, col, bg, steps, evals in [
            (22, "Donors and acceptors", T1, T1_BG,
             [("358 genomes", "14 categories"), ("852-KO panel", "presence/absence"), ("RedoxKO", "final model")],
             ["CV genus 0.576", "external 0.511", "new species 0.73", "MAGs ≥ 94%"]),
            (3, "Pathways", T2, T2_BG,
             [("3,378 genomes", "20 traits"), ("10,517 KOs", "+ KEGG modules"), ("EnergyKO-M", "final model")],
             ["CV genus 0.501", "new species 0.67"])]:
        ax.add_patch(Rectangle((0.5, y0 - 1), 99, 15, fc=bg, ec="none", alpha=0.6))
        ax.text(1.5, y0 + 12.4, name, fontsize=7.5, weight="bold", color=col, va="center")
        x = 2
        for i, (t, b) in enumerate(steps):
            box(ax, x, y0, 15, 9.5, t, b, col, "white", col)
            if i < 2:
                arrow(ax, x + 15, y0 + 4.7, x + 18, y0 + 4.7, col)
            x += 18
        arrow(ax, x - 3, y0 + 4.7, x + 1, y0 + 4.7, col)
        for j, e in enumerate(evals):
            ax.add_patch(FancyBboxPatch((x + 2 + (j % 2) * 18.5, y0 + 5.4 - (j // 2) * 5.2), 17.5, 4.0,
                                        boxstyle="round,pad=0,rounding_size=1.5", fc="white", ec=col, lw=0.7))
            ax.text(x + 3 + (j % 2) * 18.5, y0 + 7.4 - (j // 2) * 5.2, e, fontsize=5.8, va="center")
    # Bacformer add-on
    box(ax, 38, 38, 22, 6.5, "+ Bacformer-large vector", "RedoxKO-ens-BF", BF, BF_BG, BF, 6.2)
    arrow(ax, 44, 38, 44, 31.6, BF, ls="--")
    save(fig, 2, "Swimlanes: one lane per task, evaluations as chips")


# ---------------------------------------------------------------------------------------------- 3 strictness ladder
def opt3():
    fig, ax = canvas(3.2)
    label(ax, 1, 42, "a  Data, features and models")
    for i, (t, b, c, bg) in enumerate([("Genome", "RefSeq or MAG", MUTED, NBG), ("KO profile", "KOfam, pyhmmer", T1, T1_BG),
                                       ("RedoxKO / EnergyKO-M", "trained per trait", T1, T1_BG)]):
        box(ax, 1 + i * 0, 30 - i * 11, 34, 8.5, t, b, c, bg, c if c != MUTED else INK)
        if i < 2:
            arrow(ax, 18, 30 - i * 11, 18, 27.5 - i * 11, c)
    ax.text(1, 3.5, "labels: 358 genomes (donors, acceptors)\n3,378 genomes (pathways)", fontsize=5.8, color=MUTED)
    label(ax, 40, 42, "b  Accuracy falls as tests get stricter (RedoxKO)")
    sub = fig.add_axes([0.45, 0.13, 0.53, 0.68])
    names = ["random\nsplit", "genus\nCV", "family\nCV", "external\ngenera"]
    m = [0.647, 0.576, 0.560, 0.511]
    tax = [np.nan, 0.472, 0.380, 0.380]
    xs = np.arange(4)
    sub.plot(xs, m, "-o", color=T1, lw=1.2, ms=4, label="RedoxKO")
    sub.plot(xs, tax, "--s", color=MUTED, lw=1, ms=3.5, label="taxonomic baseline")
    for x, v in zip(xs, m):
        sub.text(x, v + 0.02, f"{v:.3f}", ha="center", fontsize=5.8, color=T1)
    sub.set_xticks(xs, names, fontsize=6)
    sub.set_ylabel("macro AUPRC", fontsize=6.5)
    sub.set_ylim(0.3, 0.72)
    sub.tick_params(labelsize=6)
    for s in ("top", "right"):
        sub.spines[s].set_visible(False)
    sub.legend(fontsize=5.8, frameon=False, loc="upper right")
    sub.text(3.15, 0.33, "new species: AUROC 0.73\nMAGs: ≥ 94% agreement", fontsize=5.6, color=INK, ha="right")
    save(fig, 3, "Strictness ladder: compact pipeline + accuracy vs test strictness")


# ---------------------------------------------------------------------------------------------- 4 data-set map
def opt4():
    fig, ax = canvas(3.2)
    label(ax, 1, 42, "a  Training and test sets")
    sets = [("Pathway training", 3378, T2, "train"), ("Donor/acceptor training", 358, T1, "train"),
            ("External genomes", 562, TEST, "test"), ("New species (100 + 43)", 143, TEST, "test"),
            ("New species, pathways", 99, TEST, "test"), ("MAG-isolate pairs", 372, TEST, "test")]
    y = 37
    for name, n, c, role in sets:
        L = 6 * np.log10(n)
        ax.add_patch(Rectangle((24, y - 2), L, 3.2, fc=c, ec="none", alpha=0.85))
        ax.text(23, y - 0.4, name, ha="right", va="center", fontsize=6)
        ax.text(24.6 + L, y - 0.4, f"{n:,}", va="center", fontsize=5.8)
        y -= 5.4
    ax.text(24, 4.0, "bar length: genomes (log scale)\nblue, green: training   purple: never in training",
            fontsize=5.6, color=MUTED, va="center", linespacing=1.3)
    label(ax, 52, 42, "b  Models")
    box(ax, 52, 27, 22, 10, "RedoxKO", "852 KOs, final", T1, T1_BG, T1)
    box(ax, 52, 14, 22, 10, "EnergyKO-M", "10,517 KOs + modules", T2, T2_BG, T2)
    box(ax, 52, 2, 22, 9, "RedoxKO-ens-BF", "+ Bacformer-large", BF, BF_BG, BF)
    label(ax, 77, 42, "c  Headline results")
    rows = ["CV genus  0.576", "external  0.511", "new species  0.73", "MAGs  ≥ 94%", "pathways  0.501"]
    for i, r in enumerate(rows):
        ax.text(78, 36 - i * 6, r, fontsize=6.3, va="center", color=T2 if "pathways" in r else INK)
    save(fig, 4, "Data-set map (Koblitz-style data overview) + models + results")


# ---------------------------------------------------------------------------------------------- 5 shared trunk
def opt5():
    fig, ax = canvas(3.2)
    box(ax, 1, 18, 15, 9, "Genome", "RefSeq or MAG", MUTED, NBG)
    box(ax, 21, 18, 15, 9, "KO profile", "KOfam / pyhmmer", MUTED, NBG)
    arrow(ax, 16, 22.5, 21, 22.5, MUTED)
    # branches
    box(ax, 42, 31, 17, 9, "852-KO panel", "energy KOs", T1, T1_BG, T1)
    box(ax, 42, 5, 17, 9, "Full profile", "+ KEGG modules", T2, T2_BG, T2)
    box(ax, 42, 18, 17, 9, "Protein vectors", "ESM-C, Bacformer", BF, BF_BG, BF)
    arrow(ax, 36, 24, 42, 35, T1); arrow(ax, 36, 21, 42, 9.5, T2)
    ax.plot([8.5, 8.5, 39], [27, 30, 30], color=BF, lw=0.9, ls="--")
    arrow(ax, 39, 30, 42, 25, BF, ls="--")
    box(ax, 63, 31, 15, 9, "RedoxKO", "donors, acceptors", T1, T1_BG, T1)
    box(ax, 63, 18, 15, 9, "RedoxKO-ens-BF", "fusion", BF, BF_BG, BF)
    box(ax, 63, 5, 15, 9, "EnergyKO-M", "pathways", T2, T2_BG, T2)
    arrow(ax, 59, 35.5, 63, 35.5, T1); arrow(ax, 59, 22.5, 63, 22.5, BF); arrow(ax, 59, 9.5, 63, 9.5, T2)
    arrow(ax, 59, 33, 63, 25, T1)
    for i, (t, m) in enumerate([("CV genus", "0.576"), ("external", "0.511"), ("new species", "0.73"), ("MAGs", "≥94%")]):
        ax.text(81, 39 - i * 3.6, f"{t}: {m}", fontsize=5.9, color=T1, va="center")
    ax.text(81, 9.5, "CV genus: 0.501\nnew species: 0.67", fontsize=5.9, color=T2, va="center")
    label(ax, 1, 43, "One annotation, three branches")
    save(fig, 5, "Shared trunk with branches (Flamholz-style pipeline)")


# ---------------------------------------------------------------------------------------------- 6 icon pipeline
def opt6():
    fig, ax = canvas(2.9)
    xs = [10, 33, 56, 79]
    titles = ["Genome", "KO profile", "Trained models", "Strict tests"]
    bodies = ["RefSeq or MAG", "852 or 10,517 KOs\nper genome", "RedoxKO, EnergyKO-M\n(+ Bacformer fusion)",
              "unseen genera, new\nspecies, MAGs, tools"]
    # genome icon: circle chromosome
    ax.add_patch(Circle((xs[0], 25), 7, fc="none", ec=MUTED, lw=2))
    for a in np.linspace(0, 2 * np.pi, 9)[:-1]:
        ax.add_patch(Rectangle((xs[0] + 7 * np.cos(a) - 0.8, 25 + 7 * np.sin(a) - 0.8), 1.6, 1.6, fc=T1, ec="none"))
    # KO matrix icon
    rng = np.random.default_rng(3)
    for i in range(6):
        for j in range(8):
            ax.add_patch(Rectangle((xs[1] - 8 + j * 2, 19 + i * 2), 1.8, 1.8, fc=T1 if rng.random() < 0.4 else "#e5e5e5",
                                   ec="none"))
    # model icon: small tree
    for (x0, y0, x1, y1) in [(0, 6, -5, 1), (0, 6, 5, 1), (-5, 1, -7, -4), (-5, 1, -3, -4), (5, 1, 3, -4), (5, 1, 7, -4)]:
        ax.plot([xs[2] + x0, xs[2] + x1], [25 + y0, 25 + y1], color=T2, lw=1.4)
    # tests icon: checklist
    for k in range(4):
        ax.add_patch(Rectangle((xs[3] - 6, 30 - k * 3.4), 2.2, 2.2, fc="white", ec=INK, lw=0.7))
        ax.plot([xs[3] - 5.6, xs[3] - 5, xs[3] - 3.9], [31 - k * 3.4, 30.4 - k * 3.4, 32 - k * 3.4], color=TEST, lw=1)
        ax.plot([xs[3] - 2.8, xs[3] + 6], [31.1 - k * 3.4] * 2, color=LINE, lw=1.4)
    for x, t, b in zip(xs, titles, bodies):
        ax.text(x, 13, t, ha="center", fontsize=7.2, weight="bold")
        ax.text(x, 9.5, b, ha="center", va="top", fontsize=5.8, linespacing=1.3)
    for a, b in zip(xs[:-1], xs[1:]):
        arrow(ax, a + 9, 25, b - 9, 25, MUTED)
    ax.text(50, 38, "CV genus 0.576  |  external 0.511  |  new species AUROC 0.73  |  pathways 0.501",
            ha="center", fontsize=6.3, color=INK)
    save(fig, 6, "Icon pipeline (minimal, teaser-style) with one results line")


# ---------------------------------------------------------------------------------------------- 7 data matrix (X, y)
def opt7():
    fig, ax = canvas(3.3)
    label(ax, 1, 44, "a  Genome x KO matrix and labels")
    rng = np.random.default_rng(7)
    rows, cols = 14, 22
    genus = np.repeat(np.arange(7), 2)
    for i in range(rows):
        held = genus[i] == 5
        for j in range(cols):
            on = rng.random() < 0.35
            ax.add_patch(Rectangle((2 + j * 1.5, 38 - i * 2.3), 1.35, 2.1,
                                   fc=(TEST if held else T1) if on else "#ececec", ec="none"))
        for j in range(4):
            ax.add_patch(Rectangle((38 + j * 1.6, 38 - i * 2.3), 1.45, 2.1,
                                   fc=(T2 if rng.random() < 0.3 else "#ececec"), ec="none"))
    ax.text(18, 41.2, "KOs (852 or 10,517 features)", ha="center", fontsize=5.8)
    ax.text(41, 41.2, "labels", ha="center", fontsize=5.8)
    ax.text(18, 4.5, "genomes (rows); purple rows = test genus, whole genus left out of training", ha="center",
            fontsize=5.8, color=TEST)
    label(ax, 58, 44, "b  Models")
    box(ax, 58, 30, 18, 9, "RedoxKO", "logistic reg. or XGBoost", T1, T1_BG, T1)
    box(ax, 58, 18, 18, 9, "EnergyKO-M", "pathways", T2, T2_BG, T2)
    box(ax, 58, 6, 18, 9, "RedoxKO-ens-BF", "+ Bacformer", BF, BF_BG, BF)
    label(ax, 80, 44, "c  Tests")
    for i, (t, m) in enumerate(EVAL_DA[:5]):
        ax.text(80, 38 - i * 6.2, t, fontsize=6, weight="bold")
        ax.text(80, 35.8 - i * 6.2, m, fontsize=5.7)
    save(fig, 7, "ML data-matrix view (X = KOs, y = labels, test genus highlighted)")


# ---------------------------------------------------------------------------------------------- 8 results dot plot
def opt8():
    fig, ax = canvas(3.2)
    label(ax, 1, 42, "a  Pipeline")
    box(ax, 1, 31, 24, 7, "Genome", "RefSeq or MAG", MUTED, NBG)
    box(ax, 1, 20, 24, 7, "KO profile", "KOfam / pyhmmer", T1, T1_BG, T1)
    box(ax, 1, 6, 11.5, 9, "RedoxKO", "donors,\nacceptors", T1, T1_BG, T1, 6.2)
    box(ax, 13.5, 6, 11.5, 9, "EnergyKO-M", "pathways", T2, T2_BG, T2, 6.2)
    arrow(ax, 13, 31, 13, 27, MUTED)
    arrow(ax, 9, 20, 6.5, 15, T1)
    arrow(ax, 17, 20, 19.5, 15, T2)
    label(ax, 30, 42, "b  RedoxKO and baselines on each test")
    sub = fig.add_axes([0.40, 0.12, 0.57, 0.70])
    tests = ["genus CV", "family CV", "external genera"]
    red = [0.576, 0.560, 0.511]
    tax = [0.472, 0.380, 0.380]
    freq = [0.185, 0.178, np.nan]
    y = np.arange(3)[::-1]
    sub.hlines(y, 0.1, 0.7, color="#eeeeee", lw=6)
    sub.plot(red, y, "o", color=T1, ms=6, label="RedoxKO")
    sub.plot(tax, y, "s", color=MUTED, ms=5, label="taxonomic baseline")
    sub.plot(freq, y, "^", color=LINE, ms=5, label="frequency baseline")
    sub.set_yticks(y, tests, fontsize=6.3)
    sub.set_xlabel("macro AUPRC", fontsize=6.5)
    sub.set_xlim(0.1, 0.7)
    sub.tick_params(labelsize=6)
    for s in ("top", "right", "left"):
        sub.spines[s].set_visible(False)
    sub.legend(fontsize=5.8, frameon=False, loc="lower right", ncol=3, bbox_to_anchor=(1, -0.42))
    save(fig, 8, "Results-first: short pipeline + dot plot against baselines")


# ---------------------------------------------------------------------------------------------- 9 hub
def opt9():
    fig, ax = canvas(3.4)
    cx, cy = 50, 24
    ax.add_patch(Circle((cx, cy), 8, fc=T1_BG, ec=T1, lw=1))
    ax.text(cx, cy + 1.5, "KO profile", ha="center", fontsize=7.5, weight="bold", color=T1)
    ax.text(cx, cy - 2.5, "of each genome", ha="center", fontsize=5.8)
    left = [("New genome (RefSeq or MAG)", MUTED), ("358 labeled genomes", T1), ("3,378 labeled genomes", T2)]
    right = [("RedoxKO  (donors, acceptors)", T1), ("EnergyKO-M  (pathways)", T2), ("+ Bacformer: RedoxKO-ens-BF", BF)]
    bottom = ["genus/family CV", "562 external", "100 new species", "372 MAG pairs", "7 tools, baselines"]
    for i, (t, c) in enumerate(left):
        y = 38 - i * 12
        box(ax, 2, y - 2.5, 26, 5, t, "", c, "white", c, 6.3)
        arrow(ax, 28, y, cx - 8, cy + (1 - i) * 4, c)
    for i, (t, c) in enumerate(right):
        y = 38 - i * 12
        box(ax, 72, y - 2.5, 27, 5, t, "", c, "white", c, 6.3)
        arrow(ax, cx + 8, cy + (1 - i) * 4, 72, y, c)
    for i, t in enumerate(bottom):
        x = 8 + i * 18
        ax.add_patch(FancyBboxPatch((x, 0.5), 15.5, 4, boxstyle="round,pad=0,rounding_size=1.5", fc="white", ec=TEST,
                                    lw=0.7))
        ax.text(x + 7.75, 2.5, t, ha="center", va="center", fontsize=5.8)
    ax.text(50, 6.6, "tests on genomes never seen in training", ha="center", fontsize=6, color=TEST)
    arrow(ax, cx, cy - 8, cx, 8, TEST)
    save(fig, 9, "Hub: KO profile at the centre, inputs left, models right, tests below")


# ---------------------------------------------------------------------------------------------- 10 protocol timeline
def opt10():
    fig, ax = canvas(3.0)
    label(ax, 1, 38, "a  Pipeline")
    for i, (t, b, c, bg) in enumerate([("Genome", "RefSeq\nor MAG", MUTED, NBG), ("KO profile", "KOfam", T1, T1_BG),
                                       ("Model", "per trait", T1, T1_BG)]):
        box(ax, 1 + i * 11, 25, 9.5, 8.5, t, b, c, bg, c if c != MUTED else INK, 6.2)
        if i < 2:
            arrow(ax, 10.5 + i * 11, 29.5, 12 + i * 11, 29.5, c)
    ax.text(1, 21, "RedoxKO (donors, acceptors)\nEnergyKO-M (pathways)\nRedoxKO-ens-BF (+ Bacformer)", fontsize=5.8,
            va="top", linespacing=1.4)
    label(ax, 37, 38, "b  Registered protocol and tests")
    ax.plot([38, 98], [22, 22], color=MUTED, lw=1)
    events = [(38, "17 Sep", "protocol\nregistered", INK), (46, "", "genus/family\nCV 0.576", T1),
              (55, "18 Sep", "external\nconfirmation 0.511", T1), (65, "18 Sep", "pathway task\nadded, 0.501", T2),
              (74, "", "100 new species\nAUROC 0.73", TEST), (84, "", "MAG pairs,\ntool comparison", TEST),
              (95, "17-19 Oct", "further\nregistered tests", INK)]
    for i, (x, d, t, c) in enumerate(events):
        ax.add_patch(Circle((x, 22), 0.9, fc=c, ec="none"))
        ax.text(x, 25, d, ha="center", fontsize=5.6, color=MUTED)
        ax.text(x, 19 if i % 2 == 0 else 13, t, ha="center", va="top", fontsize=5.6, color=c, linespacing=1.25)
    save(fig, 10, "Protocol timeline: pre-registered steps with headline results")


for f in (opt1, opt2, opt3, opt4, opt5, opt6, opt7, opt8, opt9, opt10):
    f()

# contact sheet
from PIL import Image, ImageDraw, ImageFont
ims = [Image.open(OUT / f"option_{n:02d}.png").convert("RGB") for n in range(1, 11)]
W = 1100
tiles = []
for n, im in enumerate(ims, 1):
    im = im.resize((W, int(im.height * W / im.width)))
    canvas_im = Image.new("RGB", (W, im.height + 44), "white")
    canvas_im.paste(im, (0, 44))
    d = ImageDraw.Draw(canvas_im)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 26)
    except OSError:
        font = ImageFont.load_default()
    d.text((8, 8), f"Option {n}: {TITLES[n]}", fill="black", font=font)
    tiles.append(canvas_im)
cols = 2
rows_ = (len(tiles) + 1) // 2
hs = [max(t.height for t in tiles[r * 2:r * 2 + 2]) for r in range(rows_)]
sheet = Image.new("RGB", (W * cols + 30, sum(hs) + 15 * rows_ + 15), "#dddddd")
y = 15
for r in range(rows_):
    for c in range(cols):
        k = r * 2 + c
        if k < len(tiles):
            sheet.paste(tiles[k], (10 + c * (W + 10), y))
    y += hs[r] + 15
sheet.save(OUT / "contact_sheet.png")
print("ok", OUT)
