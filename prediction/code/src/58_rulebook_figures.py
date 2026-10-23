"""Figures for the rulebook manuscript (docs/rulebook_gene_rules.docx).

Fig. 1  held-out vs external precision of every positive donor/acceptor rule with an external test (>= 5 external
        genomes fired on), coloured by trust level, sized by external support.
Fig. 2  precision of the whole positive rule list per trait on held-out genera: main procedure vs strictly nested
        procedure (candidate genes chosen inside each training fold).

usage: python src/58_rulebook_figures.py
Output: results/rulebook/figures/fig1_transfer.png, fig2_strict.png (+ .csv with the plotted values)
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from pathway_common import RES  # noqa: E402

OUT = RES / "rulebook" / "figures"
# Okabe-Ito colours (colour-blind safe), one per trust level, fixed order
LEVEL_COL = {"mechanistic": "#0072B2", "physiological": "#009E73", "contains unexplained conditions": "#E69F00"}
LEVEL_MARK = {"mechanistic": "o", "physiological": "s", "contains unexplained conditions": "^"}
PRETTY = {"acceptor_CO2_methanogen_acetogen": "acceptor: CO₂ (methanogens, acetogens)", "uses_O2": "oxygen use"}
NICE = lambda t: PRETTY.get(t, t.replace("donor_", "donor: ").replace("acceptor_", "acceptor: ").replace("_", " "))  # noqa: E731

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": 0.8, "savefig.dpi": 300})


def fig1():
    r = pd.read_csv(RES / "rulebook" / "da_rules.csv")
    r = r[(r.cls == "has trait") & r.external_precision.notna() & r.heldout_precision.notna()
          & (r.external_support >= 5)].copy()
    r.to_csv(OUT / "fig1_transfer.csv", index=False, float_format="%.4f")
    fig, ax = plt.subplots(figsize=(3.5, 3.4))
    ax.plot([0, 1], [0, 1], color="#9A9A9A", lw=0.8, ls="--", zorder=0)
    for lvl, col in LEVEL_COL.items():
        s = r[r.interpretation_level == lvl]
        if len(s):
            ax.scatter(s.heldout_precision, s.external_precision, s=12 + s.external_support.clip(upper=150) * 0.9,
                       c=col, marker=LEVEL_MARK[lvl], edgecolors="white", linewidths=0.6, alpha=0.9,
                       label=f"{lvl} (n = {len(s)})")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Precision on held-out genera")
    ax.set_ylabel("Precision on external genomes (novel genera)")
    ax.legend(frameon=False, fontsize=7, loc="lower left", handletextpad=0.3)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(OUT / "fig1_transfer.png")
    plt.close(fig)
    return len(r)


def fig2():
    rows = []
    for task, name in (("da", "Donors & acceptors"), ("ox", "Oxygen use")):
        f = RES / "rulebook" / f"{task}_procedure.csv"
        if not f.exists():
            continue
        p = pd.read_csv(f)
        w = p.pivot(index="trait", columns="selection", values="heldout_precision")
        fired = p.pivot(index="trait", columns="selection", values="fired")
        for t in w.index:
            if fired.loc[t].min() > 0:            # a precision exists under both procedures
                rows.append(dict(task=name, trait=NICE(t), main=w.at[t, "main"], strict=w.at[t, "strict"]))
    d = pd.DataFrame(rows).sort_values(["task", "main"])
    d.to_csv(OUT / "fig2_strict.csv", index=False, float_format="%.4f")
    lo = max(0.0, (min(d["main"].min(), d["strict"].min()) // 0.1) * 0.1 - 0.05)
    fig, ax = plt.subplots(figsize=(5.0, 0.3 * len(d) + 1.3), layout="constrained")
    y = list(range(len(d)))
    for i, r in enumerate(d.itertuples()):
        ax.plot([r.main, r.strict], [i, i], color="#BDBDBD", lw=1.5, zorder=1)
    ax.scatter(d["main"], y, color="#0072B2", s=26, zorder=2, label="main procedure")
    ax.scatter(d["strict"], y, color="#E69F00", marker="D", s=22, zorder=3,
               label="strictly nested (candidate genes chosen inside training folds)")
    ax.set_yticks(y)
    ax.set_yticklabels(["oxygen use" if r.task == "Oxygen use" else r.trait for r in d.itertuples()])
    ax.set_xlim(lo, 1.0)
    ax.set_xlabel("Held-out precision")
    fig.legend(frameon=False, fontsize=7.5, loc="outside lower center", ncol=1)
    fig.savefig(OUT / "fig2_strict.png")
    plt.close(fig)
    return len(d)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    print("fig1 rules:", fig1(), "| fig2 traits:", fig2())
