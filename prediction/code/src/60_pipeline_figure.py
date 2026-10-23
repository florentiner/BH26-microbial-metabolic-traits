"""Illustrated pipeline diagram of the released model (B2) and the best internal model (B2XL + Bacformer-large), for
docs/inorganic_prediction_report.docx (src/report/make_report.js). Definitions follow src/06_run_arms.py (arms B2,
B2XLE_large) and the Kaggle feature kernels (kaggle/esmc_bacformer_large).

usage: python src/60_pipeline_figure.py
Output: results/figures/pipelines_b2_b2xle_large.png
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle  # noqa: E402

from pathway_common import RES  # noqa: E402

GREY, BLUE, ORANGE, GREEN = "#F1F1F1", "#E3EEFA", "#FDEBD6", "#E2F2E6"
EDGE = {GREY: "#9E9E9E", BLUE: "#2A78D6", ORANGE: "#D98200", GREEN: "#2E7D32"}
INK, DARK_BLUE, DARK_GR = "#333333", "#1F5FAE", "#2E7D32"
plt.rcParams.update({"font.family": "DejaVu Sans"})
rng = np.random.default_rng(3)


# ------------------------------------------------------------------ icons (drawn inside x0..x0+w, y0..y0+h)
def icon_dna(ax, x0, y0, w, h):
    t = np.linspace(0, 1, 200)
    for ph, c in ((0, "#2A78D6"), (np.pi, "#D98200")):
        ax.plot(x0 + t * w, y0 + h / 2 + 0.38 * h * np.sin(2 * np.pi * 2 * t + ph), color=c, lw=1.1)
    for tt in np.linspace(0.04, 0.96, 14):
        a, b = np.sin(2 * np.pi * 2 * tt), np.sin(2 * np.pi * 2 * tt + np.pi)
        ax.plot([x0 + tt * w] * 2, [y0 + h / 2 + 0.38 * h * a, y0 + h / 2 + 0.38 * h * b], color="#9E9E9E", lw=0.5)


def icon_proteins(ax, x0, y0, w, h):
    cols = ["#6A9FD8", "#E7A04B", "#7BB77F", "#C084C8"]
    for r in range(3):
        n = [7, 5, 6][r]
        xs = x0 + 0.08 * w + np.arange(n) * 0.12 * w
        ys = y0 + (0.8 - r * 0.3) * h + 0.05 * h * np.sin(np.arange(n))
        ax.plot(xs, ys, color="#888888", lw=0.6, zorder=1)
        for i, (x, y) in enumerate(zip(xs, ys)):
            ax.add_patch(Circle((x, y), 0.045 * min(w, h * 2.2), color=cols[(i + r) % 4], zorder=2))


def icon_hmm(ax, x0, y0, w, h):
    cols = ["#2A78D6", "#D98200", "#2E7D32", "#C0392B"]
    for i in range(9):
        heights = rng.dirichlet(np.ones(4) * 0.6) * (0.35 + 0.6 * rng.random())
        base = y0 + 0.08 * h
        for k in np.argsort(heights):
            ax.add_patch(Rectangle((x0 + 0.06 * w + i * 0.098 * w, base), 0.08 * w, heights[k] * h * 0.85,
                                   color=cols[k], lw=0))
            base += heights[k] * h * 0.85


def icon_vector(ax, x0, y0, w, h, n=14):
    bits = rng.random(n) < 0.35
    cw = w * 0.92 / n
    for i, b in enumerate(bits):
        ax.add_patch(Rectangle((x0 + 0.04 * w + i * cw, y0 + 0.35 * h), cw * 0.9, 0.4 * h,
                               fc=DARK_BLUE if b else "white", ec="#777777", lw=0.4))
        ax.text(x0 + 0.04 * w + i * cw + cw * 0.45, y0 + 0.12 * h, "1" if b else "0", fontsize=3.6, ha="center",
                color="#555555")


def icon_trees_sigmoid(ax, x0, y0, w, h, aspect=1.0):
    t = np.linspace(-4, 4, 60)
    ax.plot(x0 + 0.05 * w + (t + 4) / 8 * 0.38 * w, y0 + 0.12 * h + 0.76 * h / (1 + np.exp(-2 * t)), color=DARK_BLUE,
            lw=1.2)
    for cx in (0.62, 0.86):
        nodes = [(cx, 0.85), (cx - 0.08, 0.5), (cx + 0.08, 0.5), (cx - 0.12, 0.15), (cx - 0.04, 0.15)]
        for a, b in ((0, 1), (0, 2), (1, 3), (1, 4)):
            ax.plot([x0 + nodes[a][0] * w, x0 + nodes[b][0] * w], [y0 + nodes[a][1] * h, y0 + nodes[b][1] * h],
                    color="#666666", lw=0.6)
        for nx, ny in nodes:
            ax.add_patch(Circle((x0 + nx * w, y0 + ny * h), 0.022 * w / aspect, color="#7BB77F" if ny < 0.3 else DARK_BLUE))


def icon_output(ax, x0, y0, w, h):
    vals = [0.92, 0.15, 0.71, 0.05, 0.38]
    for i, v in enumerate(vals):
        ax.add_patch(Rectangle((x0 + 0.1 * w + i * 0.17 * w, y0 + 0.1 * h), 0.12 * w, v * 0.8 * h,
                               color=DARK_GR if v >= 0.5 else "#B9D8BE", lw=0))
    ax.plot([x0 + 0.05 * w, x0 + 0.95 * w], [y0 + 0.1 * h + 0.5 * 0.8 * h] * 2, color="#C0392B", lw=0.7, ls="--")


def icon_embedding(ax, x0, y0, w, h):
    xs = x0 + 0.08 * w + np.arange(5) * 0.09 * w
    ax.plot(xs, [y0 + 0.5 * h] * 5, color="#888888", lw=0.6)
    for i, x in enumerate(xs):
        ax.add_patch(Circle((x, y0 + 0.5 * h), 0.03 * w,
                            color=["#6A9FD8", "#E7A04B", "#7BB77F", "#C084C8", "#6A9FD8"][i]))
    ax.annotate("", xy=(x0 + 0.66 * w, y0 + 0.5 * h), xytext=(x0 + 0.5 * w, y0 + 0.5 * h),
                arrowprops=dict(arrowstyle="-|>", lw=0.7, color="#555555", mutation_scale=5))
    cm = plt.get_cmap("Oranges")
    for i in range(8):
        ax.add_patch(Rectangle((x0 + 0.72 * w, y0 + 0.08 * h + i * 0.105 * h), 0.14 * w, 0.1 * h,
                               color=cm(0.25 + 0.7 * rng.random()), lw=0))


def icon_transformer(ax, x0, y0, w, h):
    for i in range(4):
        ax.add_patch(Rectangle((x0 + 0.1 * w, y0 + 0.05 * h + i * 0.17 * h), 0.8 * w, 0.12 * h, fc="#F6C98F",
                               ec="#B36B00", lw=0.5))
    xs = x0 + np.linspace(0.18, 0.82, 6) * w
    top = y0 + 0.78 * h
    for i in range(6):
        for j in range(i + 1, 6):
            if rng.random() < 0.45:
                t = np.linspace(0, np.pi, 20)
                ax.plot((xs[i] + xs[j]) / 2 + (xs[j] - xs[i]) / 2 * np.cos(t),
                        top + (0.05 + 0.15 * (j - i) / 5) * h * np.sin(t), color="#B36B00", lw=0.4)
    for x in xs:
        ax.add_patch(Circle((x, top), 0.022 * w, color="#B36B00"))


def icon_pool(ax, x0, y0, w, h):
    cm = plt.get_cmap("Oranges")
    for c in range(4):
        for i in range(6):
            ax.add_patch(Rectangle((x0 + 0.06 * w + c * 0.1 * w, y0 + 0.1 * h + i * 0.13 * h), 0.08 * w, 0.12 * h,
                                   color=cm(0.2 + 0.6 * rng.random()), lw=0))
    ax.text(x0 + 0.6 * w, y0 + 0.62 * h, "mean", fontsize=4.2, ha="center", color="#555555")
    ax.annotate("", xy=(x0 + 0.74 * w, y0 + 0.45 * h), xytext=(x0 + 0.47 * w, y0 + 0.45 * h),
                arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#555555", mutation_scale=5))
    for i in range(6):
        ax.add_patch(Rectangle((x0 + 0.78 * w, y0 + 0.1 * h + i * 0.13 * h), 0.1 * w, 0.12 * h,
                               color=cm(0.3 + 0.6 * rng.random()), lw=0))


# ------------------------------------------------------------------ boxes
def box(ax, x, y, w, h, title, detail, fill, icon=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.01", fc=fill,
                                ec=EDGE[fill], lw=0.9))
    ih, iw = h * 0.36, min(w * 0.8, 0.13)
    if icon:
        icon(ax, x + (w - iw) / 2, y + h - ih - h * 0.06, iw, ih)
    ax.text(x + w / 2, y + h * 0.44, title, ha="center", va="center", fontsize=6.2, weight="bold", color=INK)
    ax.text(x + w / 2, y + h * 0.2, detail, ha="center", va="center", fontsize=4.8, color="#444444",
            linespacing=1.18)
    return (x, y, w, h)


def arrow(ax, a, b):
    ax_, ay, aw, ah = a
    bx, by, bw, bh = b
    ax.annotate("", xy=(bx, by + bh / 2), xytext=(ax_ + aw, ay + ah / 2),
                arrowprops=dict(arrowstyle="-|>", lw=0.9, color="#555555", shrinkA=1, shrinkB=1, mutation_scale=8))


def main():
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    ax.set_xlim(-0.005, 1.005)
    ax.set_ylim(-0.01, 1.04)
    ax.axis("off")
    h, w, g = 0.25, 0.152, 0.018
    X = [i * (w + g) for i in range(6)]
    # ---- (a) B2
    ax.text(0.0, 1.035, "a   B2 — released model (CPU only, about 2 minutes per genome)", fontsize=7.6,
            weight="bold", va="top")
    y = 0.73
    A = [box(ax, X[0], y, w, h, "Genome", "RefSeq accession\nor nucleotide FASTA", GREY, icon_dna),
         box(ax, X[1], y, w, h, "Proteins", "NCBI annotation, or\ngene calls by Pyrodigal", GREY, icon_proteins),
         box(ax, X[2], y, w, h, "KofamScan", "KEGG KOfam profile\nHMMs, per-KO\nadaptive thresholds", BLUE, icon_hmm),
         box(ax, X[3], y, w, h, "852 KOs: 1 / 0", "energy-metabolism panel,\npresence / absence", BLUE, icon_vector),
         box(ax, X[4], y, w, h, "Model per trait", "L2 logistic regression\nor XGBoost (depth 2 or 4);\n"
             "best of 6 in inner CV", BLUE, icon_trees_sigmoid),
         box(ax, X[5], y, w, h, "P(trait) → call", "cut-off maximising MCC\n(inner CV); 14 traits", GREEN, icon_output)]
    for p, q in zip(A, A[1:]):
        arrow(ax, p, q)
    # ---- (b) B2XL + Bacformer-large
    ax.text(0.0, 0.675, "b   B2XL + Bacformer-large — best internal model (GPU needed for the language models)",
            fontsize=7.6, weight="bold", va="top")
    yt, yb = 0.35, 0.0
    ym = (yt + yb) / 2
    b1 = box(ax, X[0], ym, w, h, "Genome", "RefSeq accession\nor nucleotide FASTA", GREY, icon_dna)
    b2 = box(ax, X[1], ym, w, h, "Proteins", "in genome order\n(contig, position)", GREY, icon_proteins)
    k1 = box(ax, X[2], yt, w, h, "852 KOs: 1 / 0", "KofamScan; the same\nKO vector as in B2", BLUE, icon_vector)
    k2 = box(ax, X[3], yt, 2 * w + g, h, "KO model",
             "L2 logistic regression + XGBoost on the 852 KOs,\nscores averaged; C and depth chosen in inner CV", BLUE,
             icon_trees_sigmoid)
    e1 = box(ax, X[2], yb, w, h, "ESM-C 300M", "protein language model:\none 960-number vector\nper protein", ORANGE,
             icon_embedding)
    e2 = box(ax, X[3], yb, w, h, "Bacformer-large", "30-layer transformer over\nall proteins in genome\norder "
             "(6 000 per window)", ORANGE, icon_transformer)
    e3 = box(ax, X[4], yb, w, h, "Genome vector", "mean over proteins →\n960 numbers → L2\nlogistic regression", ORANGE,
             icon_pool)
    m = box(ax, X[5], ym, w, h, "½ + ½ → P(trait)", "average of the KO model\nand the Bacformer model;\nMCC cut-off",
            GREEN, icon_output)
    for p, q in ((b1, b2), (b2, k1), (b2, e1), (k1, k2), (e1, e2), (e2, e3), (k2, m), (e3, m)):
        arrow(ax, p, q)
    out = RES.parent / "results" / "figures"
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "pipelines_b2_b2xle_large.png", dpi=300, bbox_inches="tight", pad_inches=0.03)
    print("written", out / "pipelines_b2_b2xle_large.png")


if __name__ == "__main__":
    main()
