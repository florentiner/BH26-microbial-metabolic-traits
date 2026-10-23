"""Metrics for the pathway arms: same estimator, bootstrap and inference as the donor/acceptor task.

Bootstrap unit = split group (genus or family); 1000 resamples shared across arms and traits at a level, so paired
contrasts use identical resamples. Outputs into results/pathways/.
usage: python src/33_pathway_metrics.py [--out-suffix ""]
"""
import argparse
import importlib
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef

from pathway_common import RES, load_pathways
from redox_common import LEVELS, point_metrics, safe_ap

m8 = importlib.import_module("08_metrics")
bh, boot_p, ci = m8.bh, m8.boot_p, m8.ci
N_BOOT, SEED = 1000, 20260918
NULLS = ("null_prevalence", "null_taxNN")
CONTRASTS = [("PB2", "PA"), ("PB2", "PB1"), ("PB3", "PB2"), ("PB1", "PA"), ("PB3", "PA"),
             ("PA", "null_taxNN"), ("PB1", "null_taxNN"), ("PB2", "null_taxNN"), ("PB3", "null_taxNN"),
             ("PE", "PB2"), ("PE", "null_taxNN")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--oof-dir", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--subset", default=None, choices=["no_faprotax", "informative_negatives"],
                    help="restrict the evaluation set (models are unchanged; docs/deviations.md #19). "
                         "no_faprotax: genomes whose Madin record is not FAPROTAX-derived (labels not "
                         "taxonomy-propagated). informative_negatives: genomes with >= 2 annotated pathway terms, "
                         "where an absent term is weaker evidence of a true negative")
    args = ap.parse_args()
    oof_dir = args.oof_dir or RES / "pathways/oof"
    out = RES / (args.out or ("pathways" if not args.subset else f"pathways/sensitivity_{args.subset}"))
    out.mkdir(parents=True, exist_ok=True)

    lab, traits, Y = load_pathways()
    if args.subset:
        if args.subset == "no_faprotax":
            keep = ~lab.faprotax_source.astype(bool).to_numpy()
        else:  # informative_negatives
            n_terms = lab.pathways.fillna("").str.split(",").apply(lambda xs: len([x for x in xs if x.strip()]))
            keep = (n_terms >= 2).to_numpy()
        print(f"{args.subset}: evaluating {keep.sum()} / {len(keep)} genomes")
        lab, Y = lab[keep].reset_index(drop=True), Y[keep]
        # a trait with too few positives left is dropped from this sensitivity analysis, not from the primary one
        traits = [t for j, t in enumerate(traits) if Y[:, j].sum() >= 10 and (Y[:, j] == 0).sum() >= 10]
        Y = np.column_stack([Y[:, i] for i, t in enumerate(load_pathways()[1]) if t in traits])
    oof = pd.concat([pd.read_parquet(f) for f in sorted(Path(oof_dir).glob("*.parquet"))], ignore_index=True)
    oof = oof[oof.genome_id.isin(set(lab.genome_id)) & oof.trait.isin(set(traits))]
    arms = [a for a in oof.arm.unique() if a not in NULLS] + [a for a in NULLS if a in set(oof.arm)]
    gid = {g: i for i, g in enumerate(lab.genome_id)}
    rng = np.random.default_rng(SEED)
    per_trait, macro, paired, per_fold = [], [], [], []

    for level in LEVELS:
        groups = lab[f"group_{level}"].to_numpy()
        ug = np.unique(groups)
        members = {g: np.where(groups == g)[0] for g in ug}
        boots = [np.concatenate([members[g] for g in rng.choice(ug, len(ug), replace=True)]) for _ in range(N_BOOT)]
        P, T, F, CFG = {}, {}, {}, {}
        for arm in arms:
            d = oof[(oof.arm == arm) & (oof.level == level)]
            if d.empty:
                continue
            probs, thr = np.full(Y.shape, np.nan), np.full(Y.shape, np.nan)
            fold, ychk, cfg = np.full(Y.shape, -1), np.full(Y.shape, np.nan), {}
            for j, t in enumerate(traits):
                dt = d[d.trait == t]
                if dt.empty:
                    continue
                idx = dt.genome_id.map(gid).to_numpy()
                probs[idx, j], thr[idx, j], fold[idx, j], ychk[idx, j] = dt.prob, dt.threshold, dt.fold, dt.y
                cfg[t] = dt.config.value_counts().index[0]
            if np.isnan(probs).any():
                print(f"skipping {arm}/{level}: {int(np.isnan(probs).sum())} labels unscored")
                continue
            assert np.array_equal(ychk, Y), f"{arm}/{level}: OOF labels differ from label matrix"
            P[arm], T[arm], F[arm], CFG[arm] = probs, thr, fold, cfg
        B = {a: np.stack([np.array([safe_ap(Y[b, j], P[a][b, j]) for b in boots]) for j in range(len(traits))], 1)
             for a in P}
        for arm in P:
            rows = []
            for j, t in enumerate(traits):
                y, p, th = Y[:, j], P[arm][:, j], T[arm][:, j]
                pm = point_metrics(y, p, th)
                has_thr = not np.all(np.isnan(th))
                bmcc = [matthews_corrcoef(Y[b, j], P[arm][b, j] >= T[arm][b, j]) for b in boots] if has_thr else []
                d_null = B[arm][:, j] - B["null_taxNN"][:, j]
                row = dict(trait=t, arm=arm, split_level=level, n_pos=int(y.sum()), n_neg=int((y == 0).sum()),
                           prevalence=y.mean(), auprc=pm["auprc"])
                row["auprc_ci_low"], row["auprc_ci_high"] = ci(B[arm][:, j])
                row["auprc_lift"] = pm["auprc"] / y.mean()
                row["auroc"], row["mcc"] = pm["auroc"], pm["mcc"]
                row["mcc_ci_low"], row["mcc_ci_high"] = ci(bmcc) if has_thr else (np.nan, np.nan)
                row["f1"], row["brier"] = pm["f1"], pm["brier"]
                row["null_taxNN_auprc"] = safe_ap(y, P["null_taxNN"][:, j])
                row["delta_vs_null"] = pm["auprc"] - row["null_taxNN_auprc"]
                row["delta_vs_null_ci_low"], row["delta_vs_null_ci_high"] = ci(d_null)
                row["delta_vs_null_p_boot"] = boot_p(d_null) if arm != "null_taxNN" else np.nan
                row["modal_config"] = CFG[arm].get(t, "")
                rows.append(row)
                for k in sorted(set(F[arm][:, j])):
                    fk = F[arm][:, j] == k
                    per_fold.append(dict(trait=t, arm=arm, split_level=level, fold=int(k),
                                         n_pos=int(Y[fk, j].sum()), auprc=safe_ap(Y[fk, j], P[arm][fk, j])))
            for r, qq in zip(rows, bh([r["delta_vs_null_p_boot"] for r in rows])):
                r["delta_vs_null_q_bh"] = qq
            per_trait += rows
            mac = np.nanmean(B[arm], axis=1)
            lo, hi = ci(mac)
            macro.append(dict(arm=arm, split_level=level, n_traits=len(traits),
                              macro_auprc=np.mean([r["auprc"] for r in rows]), macro_auprc_ci_low=lo,
                              macro_auprc_ci_high=hi, macro_auroc=np.mean([r["auroc"] for r in rows]),
                              macro_mcc=np.mean([r["mcc"] for r in rows]),
                              macro_brier=np.mean([r["brier"] for r in rows]),
                              macro_delta_vs_null=np.mean([r["delta_vs_null"] for r in rows]),
                              n_traits_better_than_null_q_lt_0_05=int(sum(
                                  r["delta_vs_null_q_bh"] < 0.05 and r["delta_vs_null"] > 0 for r in rows))))
        for a, b in CONTRASTS:
            if a not in P or b not in P:
                continue
            dm = np.nanmean(B[a], axis=1) - np.nanmean(B[b], axis=1)
            lo, hi = ci(dm)
            point = np.mean([safe_ap(Y[:, j], P[a][:, j]) - safe_ap(Y[:, j], P[b][:, j]) for j in range(len(traits))])
            paired.append(dict(split_level=level, contrast=f"{a} - {b}", trait="MACRO", delta_auprc=point, ci_low=lo,
                               ci_high=hi, p_boot=boot_p(dm), q_bh=np.nan, tied=bool(lo <= 0 <= hi)))
            rows = []
            for j, t in enumerate(traits):
                dj = B[a][:, j] - B[b][:, j]
                lo, hi = ci(dj)
                rows.append(dict(split_level=level, contrast=f"{a} - {b}", trait=t,
                                 delta_auprc=safe_ap(Y[:, j], P[a][:, j]) - safe_ap(Y[:, j], P[b][:, j]),
                                 ci_low=lo, ci_high=hi, p_boot=boot_p(dj), tied=bool(lo <= 0 <= hi)))
            for r, qq in zip(rows, bh([r["p_boot"] for r in rows])):
                r["q_bh"] = qq
            paired += rows
        print(f"{level}: done", flush=True)

    pd.DataFrame(per_trait).to_csv(out / "metrics_per_trait.csv", index=False, float_format="%.4f")
    pd.DataFrame(macro).to_csv(out / "metrics_macro.csv", index=False, float_format="%.4f")
    pd.DataFrame(per_fold).to_csv(out / "metrics_per_fold.csv", index=False, float_format="%.4f")
    pd.DataFrame(paired).to_csv(out / "paired_comparisons.csv", index=False, float_format="%.4f")
    mm = pd.DataFrame(macro)
    print(mm[["arm", "split_level", "macro_auprc", "macro_auprc_ci_low", "macro_auprc_ci_high",
              "macro_delta_vs_null", "n_traits_better_than_null_q_lt_0_05"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
