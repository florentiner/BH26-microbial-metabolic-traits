"""Paper figures and the taxonomic-distance decay table (primary analysis).

Fig 1  macro AUPRC per arm, genus vs family hold-out, 95% cluster-bootstrap CI; null reference lines.
Fig 2  per-trait AUPRC minus taxonomic null for the released arm (genus hold-out), sorted, with 95% CI.
Fig 3  decay: macro AUPRC by the lowest taxonomic rank a test genome shares with its training fold (genus hold-out).
Fig 4  per-trait AUPRC heatmap, arms x traits (genus hold-out).
Output: results/figures/fig{1..4}_*.{png,pdf}, results/decay_by_rank.csv
"""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from redox_common import RES, load_labelled, safe_ap  # noqa: E402

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # validated categorical slots 1-4 (light)
TEXT, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
ARM_ORDER = ["A", "B1", "B1r", "B2", "C", "D", "E", "E_MAG", "ENS", "B2C", "B2E", "B2LGE", "null_taxNN",
             "null_prevalence"]
ARM_LABEL = {"A": "A rules", "B1": "B1 markers (provided)", "B1r": "B1r markers (recomputed)",
             "B2": "B2 energy KO panel", "C": "C pooled ESM-2", "D": "D ABMIL", "E": "E Bacformer",
             "E_MAG": "E Bacformer (MAG ckpt)", "ENS": "B2 + LithoGenie", "B2C": "B2 + ESM-2",
             "B2E": "B2 + Bacformer (released)", "B2LGE": "B2 + LithoGenie + Bacformer",
             "null_taxNN": "taxonomic null", "null_prevalence": "prevalence null"}
RANKS = ["family", "order", "class", "phylum", "domain"]

plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": TEXT, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "savefig.dpi": 300, "savefig.bbox": "tight", "pdf.fonttype": 42})


def save(fig, name):
    out = RES / "figures"
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"{name}.{ext}")
    plt.close(fig)


def fig_macro(macro):
    arms = [a for a in ARM_ORDER if a in set(macro.arm) and not a.startswith("null")]
    fig, ax = plt.subplots(figsize=(8.4, 3.4))
    x = np.arange(len(arms))
    for k, (lvl, off) in enumerate((("genus", -0.12), ("family", 0.12))):
        d = macro[macro.split_level == lvl].set_index("arm").loc[arms]
        ax.errorbar(x + off, d.macro_auprc, yerr=[d.macro_auprc - d.macro_auprc_ci_low,
                                                  d.macro_auprc_ci_high - d.macro_auprc],
                    fmt="o", ms=5, color=SERIES[k], elinewidth=1.5, capsize=0, label=f"held-out {lvl}")
        null = macro[(macro.split_level == lvl) & (macro.arm == "null_taxNN")].macro_auprc.iloc[0]
        ax.axhline(null, color=SERIES[k], lw=1, ls="--")
        ax.text(len(arms) - 0.5, null, f"taxonomic null ({lvl})", color=MUTED, fontsize=7, va="bottom", ha="right")
    ax.set_xticks(x, [ARM_LABEL[a] for a in arms], rotation=30, ha="right")
    ax.set_ylabel("macro AUPRC (14 traits)")
    ax.legend(frameon=False, loc="upper left")
    save(fig, "fig1_macro_auprc")


def fig_delta(per_trait, arm):
    d = per_trait[(per_trait.arm == arm) & (per_trait.split_level == "genus")].sort_values("delta_vs_null")
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    y = np.arange(len(d))
    ax.errorbar(d.delta_vs_null, y, xerr=[d.delta_vs_null - d.delta_vs_null_ci_low,
                                          d.delta_vs_null_ci_high - d.delta_vs_null],
                fmt="o", ms=5, color=SERIES[0], elinewidth=1.5, capsize=0)
    ax.axvline(0, color=MUTED, lw=1)
    ax.set_yticks(y, d.trait)
    ax.set_xlabel(f"AUPRC − taxonomic null  ({ARM_LABEL[arm]}, held-out genus)")
    for yi, (_, r) in zip(y, d.iterrows()):
        if r.delta_vs_null_q_bh < 0.05:
            ax.text(r.delta_vs_null_ci_high + 0.01, yi, "q<0.05", fontsize=6, color=MUTED, va="center")
    save(fig, "fig2_delta_vs_null")


def decay_table(lab, traits, Y, arms):
    """For genus hold-out: lowest rank shared by each test genome with any genome in its training fold."""
    fold = lab.fold_genus.to_numpy()
    shared = np.empty(len(lab), dtype=object)
    for k in np.unique(fold):
        tr = fold != k
        for i in np.where(fold == k)[0]:
            shared[i] = "none"
            for r in RANKS:
                v = lab.at[i, r]
                if pd.notna(v) and (lab.loc[tr, r] == v).any():
                    shared[i] = r
                    break
    gid = {g: i for i, g in enumerate(lab.genome_id)}
    rows = []
    for arm in arms:
        d = pd.read_parquet(RES / "oof" / (f"{arm}.parquet" if not arm.startswith("null") else "nulls.parquet"))
        d = d[(d.arm == arm) & (d.level == "genus")]
        P = np.full(Y.shape, np.nan)
        for j, t in enumerate(traits):
            dt = d[d.trait == t]
            P[dt.genome_id.map(gid).to_numpy(), j] = dt.prob
        for r in RANKS + ["none"]:
            m = shared == r
            if m.sum() < 10:
                continue
            aps = [safe_ap(Y[m, j], P[m, j]) for j in range(len(traits))]
            rows.append(dict(arm=arm, lowest_shared_rank=r, n_genomes=int(m.sum()),
                             n_traits_evaluable=int(np.sum(~np.isnan(aps))), macro_auprc=np.nanmean(aps),
                             macro_prevalence=np.nanmean([Y[m, j].mean() for j in range(len(traits))
                                                          if not np.isnan(aps[j])])))
    out = pd.DataFrame(rows)
    out.to_csv(RES / "decay_by_rank.csv", index=False, float_format="%.4f")
    return out


def fig_decay(dec, arms):
    fig, ax = plt.subplots(figsize=(4.8, 3.2))
    order = [r for r in RANKS + ["none"] if r in set(dec.lowest_shared_rank)]
    for k, arm in enumerate(arms):
        d = dec[dec.arm == arm].set_index("lowest_shared_rank").reindex(order)
        ax.plot(range(len(order)), d.macro_auprc, marker="o", ms=5, lw=2, color=SERIES[k], label=ARM_LABEL[arm])
    n = dec.drop_duplicates("lowest_shared_rank").set_index("lowest_shared_rank").reindex(order).n_genomes
    ax.set_xticks(range(len(order)), [f"{r}\n(n={int(c)})" for r, c in zip(order, n)])
    ax.set_xlabel("lowest taxonomic rank shared with training genomes (held-out genus)")
    ax.set_ylabel("macro AUPRC")
    ax.legend(frameon=False, fontsize=7)
    save(fig, "fig3_decay_by_rank")


def fig_heatmap(per_trait):
    d = per_trait[per_trait.split_level == "genus"]
    arms = [a for a in ARM_ORDER if a in set(d.arm)]
    M = d.pivot(index="trait", columns="arm", values="auprc")[arms]
    fig, ax = plt.subplots(figsize=(9.6, 4.6))
    im = ax.imshow(M.to_numpy(), cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(arms)), [ARM_LABEL[a] for a in arms], rotation=35, ha="right")
    ax.set_yticks(range(len(M)), M.index)
    ax.grid(False)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.iat[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6, color="white" if v > 0.6 else TEXT)
    fig.colorbar(im, ax=ax, fraction=0.03, label="AUPRC (held-out genus)")
    save(fig, "fig4_auprc_heatmap")


def main():
    lab, traits, Y = load_labelled()
    macro = pd.read_csv(RES / "metrics_macro.csv")
    per_trait = pd.read_csv(RES / "metrics_per_trait.csv")
    released_file = RES / "models" / "RELEASED.json"
    released = json.loads(released_file.read_text())["released_arm"] if released_file.exists() else "B2"
    fig_macro(macro)
    fig_delta(per_trait, released)
    arms_present = set(macro.arm)
    decay_arms = [a for a in (released, "B2", "E", "null_taxNN") if a in arms_present]
    decay_arms = list(dict.fromkeys(decay_arms))[:4]
    fig_decay(decay_table(lab, traits, Y, decay_arms), decay_arms)
    fig_heatmap(per_trait)
    print("figures written")


if __name__ == "__main__":
    main()
