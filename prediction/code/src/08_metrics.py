"""Aggregate out-of-fold predictions into metric tables with cluster-bootstrap CIs.

usage: python src/08_metrics.py [--split frozen] [--label-mode main] [--subset quality_ok] [--out DIR]

Bootstrap unit = split group (genus or family). 1000 resamples are shared by every arm and trait at a level,
so paired differences use identical resamples. The bootstrap resamples the evaluation set only (training and
partition variance are not included; partition variance is examined in sensitivity S3).
AUPRC is computed on predictions pooled over the 5 outer folds (primary); mean per-fold AUPRC is in
metrics_per_fold.csv. Per-trait p-values: two-sided bootstrap p = 2*min(P(delta<=0), P(delta>=0)), BH-adjusted
within each (contrast, split level) family of 14 traits. Macro contrasts are the primary inferential results.

Outputs (in results/ for the primary analysis, otherwise results/sensitivity/<split>__<label-mode>[__<subset>]/):
  metrics_per_trait.csv, metrics_macro.csv, metrics_per_fold.csv, paired_comparisons.csv, selected_configs.csv,
  run_manifest.json
"""
import argparse
import json
import subprocess

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef

from redox_common import LEVELS, PROC, RES, ROOT, load_labelled, point_metrics, safe_ap, sha256

N_BOOT = 1000
SEED = 20260917
NULLS = ("null_prevalence", "null_taxNN")
CONTRASTS = [  # (a, b): a - b
    ("B1", "B1r"), ("B2", "B1r"), ("B1r", "A"), ("B2", "A"), ("C", "B2"), ("C", "B1r"), ("D", "C"), ("D", "B2"),
    ("E", "C"), ("E", "B2"), ("E_MAG", "E"),
    ("A", "null_taxNN"), ("B1", "null_taxNN"), ("B1r", "null_taxNN"), ("B2", "null_taxNN"), ("C", "null_taxNN"),
    ("D", "null_taxNN"), ("E", "null_taxNN"), ("E_MAG", "null_taxNN"),
    ("B2E", "B2"), ("B2E", "null_taxNN"), ("B2C", "B2"), ("ENS", "B2"), ("B2LGE", "B2E"),
    ("B2XL", "B2"), ("B2XL", "null_taxNN"), ("B2E", "B2XL"), ("B2XLE", "B2XL"), ("B2XLC", "B2XL"), ("B2B1r", "B2"),
    ("B2RP", "B2"),
]
THRESHOLD_SOURCE = {"A": "theta in {>0,0.5,1} maximising MCC on outer training fold",
                    "D": "inner grouped CV (5-seed ensemble), max MCC"}


def ci(a):
    a = np.asarray(a, dtype=float)
    a = a[~np.isnan(a)]
    return (np.percentile(a, 2.5), np.percentile(a, 97.5)) if len(a) else (np.nan, np.nan)


def boot_p(d):
    d = np.asarray(d, dtype=float)
    d = d[~np.isnan(d)]
    if not len(d):
        return np.nan
    # (k + 1) / (n + 1) correction: never reports p = 0 from a finite bootstrap (minimum 2/1001)
    k = min((d <= 0).sum(), (d >= 0).sum())
    return float(min(1.0, 2 * (k + 1) / (len(d) + 1)))


def bh(p):
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(pv)
    out[order] = np.minimum(ranked, 1)
    q[ok] = out
    return q


def n_seeds(config):
    model = json.loads(config).get("model") if config.startswith("{") else None
    return {"xgb": 3, "abmil": 5}.get(model, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", default="frozen")
    ap.add_argument("--label-mode", default="main")
    ap.add_argument("--subset", default=None, help="evaluation restricted to genomes with column==True in "
                                                   "data/processed/genome_qc.tsv, e.g. quality_ok")
    args = ap.parse_args()
    primary = args.split == "frozen" and args.label_mode == "main"
    oof_dir = RES / "oof" if primary else RES / "sensitivity" / f"{args.split}__{args.label_mode}" / "oof"
    out = RES if primary and not args.subset else \
        RES / "sensitivity" / (f"{args.split}__{args.label_mode}" + (f"__{args.subset}" if args.subset else ""))
    out.mkdir(parents=True, exist_ok=True)

    lab, traits, Y = load_labelled(args.split, args.label_mode)
    if args.subset:
        qc = pd.read_csv(PROC / "genome_qc.tsv", sep="\t").set_index("genome_id")[args.subset].astype(bool)
        keep = lab.genome_id.map(qc).fillna(False).to_numpy(dtype=bool)
        Y = Y.copy()
        Y[~keep] = np.nan
        print(f"subset {args.subset}: evaluating {keep.sum()} / {len(keep)} genomes")
    obs = ~np.isnan(Y)
    oof = pd.concat([pd.read_parquet(f) for f in sorted(oof_dir.glob("*.parquet"))], ignore_index=True)
    arms = [a for a in oof.arm.unique() if a not in NULLS] + list(NULLS)
    gid_index = {g: i for i, g in enumerate(lab.genome_id)}
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
            probs = np.full(Y.shape, np.nan)
            thr = np.full(Y.shape, np.nan)
            fold = np.full(Y.shape, -1)
            ychk = np.full(Y.shape, np.nan)
            cfg = {}
            for j, t in enumerate(traits):
                dt = d[d.trait == t]
                idx = dt.genome_id.map(gid_index).to_numpy()
                probs[idx, j], thr[idx, j], fold[idx, j], ychk[idx, j] = dt.prob, dt.threshold, dt.fold, dt.y
                cfg[t] = dt.config.value_counts().index[0]
            # every observed label must have a prediction and agree with the label matrix; arms with partial genome
            # coverage (B3, B3ENS: genomes absent from TogoMCP) are reported in research_ledger.csv instead
            if np.isnan(probs[obs]).any():
                print(f"skipping {arm}/{level}: partial coverage ({int(np.isnan(probs[obs]).sum())} labels unscored)")
                continue
            assert np.array_equal(ychk[obs], Y[obs]), f"{arm}/{level}: OOF labels differ from label matrix"
            P[arm], T[arm], F[arm], CFG[arm] = probs, thr, fold, cfg

        def ap_boot(arm, j):
            return np.array([safe_ap(Y[b, j][obs[b, j]], P[arm][b, j][obs[b, j]]) for b in boots])

        B = {arm: np.stack([ap_boot(arm, j) for j in range(len(traits))], axis=1) for arm in P}
        for arm in P:
            rows = []
            for j, t in enumerate(traits):
                m = obs[:, j]
                y, p, th = Y[m, j], P[arm][m, j], T[arm][m, j]
                pm = point_metrics(y, p, th)
                has_thr = not np.all(np.isnan(th))
                boot_mcc = [matthews_corrcoef(Y[b, j][obs[b, j]], P[arm][b, j][obs[b, j]] >= T[arm][b, j][obs[b, j]])
                            for b in boots] if has_thr else []
                d_null = B[arm][:, j] - B["null_taxNN"][:, j]
                null_tax = safe_ap(y, P["null_taxNN"][m, j])
                row = dict(trait=t, arm=arm, split_level=level, n_pos=int(y.sum()), n_neg=int((y == 0).sum()),
                           n_unobserved=int((~m).sum()), prevalence=y.mean(), auprc=pm["auprc"])
                row["auprc_ci_low"], row["auprc_ci_high"] = ci(B[arm][:, j])
                row["auprc_lift"] = pm["auprc"] / y.mean()
                row["auroc"] = pm["auroc"]
                row["mcc"] = pm["mcc"]
                row["mcc_ci_low"], row["mcc_ci_high"] = ci(boot_mcc) if has_thr else (np.nan, np.nan)
                row["f1"], row["brier"] = pm["f1"], pm["brier"]
                row["null_prevalence_auprc_analytic"] = y.mean()
                row["null_prevalence_auprc_pooled"] = safe_ap(y, P["null_prevalence"][m, j])
                row["null_taxNN_auprc"] = null_tax
                row["delta_vs_null"] = pm["auprc"] - null_tax
                row["delta_vs_null_ci_low"], row["delta_vs_null_ci_high"] = ci(d_null)
                row["delta_vs_null_p_boot"] = boot_p(d_null) if arm != "null_taxNN" else np.nan
                row["n_seeds"] = n_seeds(CFG[arm][t])
                row["modal_config"] = CFG[arm][t]
                row["threshold_source"] = ("none (ranking only)" if not has_thr else
                                           THRESHOLD_SOURCE.get(arm, "inner grouped CV, max MCC"))
                row["n_bootstrap_valid"] = int((~np.isnan(B[arm][:, j])).sum())
                rows.append(row)
                for k in sorted(set(F[arm][m, j])):
                    fk = (F[arm][:, j] == k) & obs[:, j]
                    per_fold.append(dict(trait=t, arm=arm, split_level=level, fold=int(k),
                                         n_pos=int(Y[fk, j].sum()), auprc=safe_ap(Y[fk, j], P[arm][fk, j])))
            q = bh([r["delta_vs_null_p_boot"] for r in rows])
            for r, qq in zip(rows, q):
                r["delta_vs_null_q_bh"] = qq
            per_trait += rows
            mac = np.nanmean(B[arm], axis=1)
            lo, hi = ci(mac)
            macro.append(dict(arm=arm, split_level=level, n_traits=len(traits),
                              macro_auprc=np.mean([r["auprc"] for r in rows]), macro_auprc_ci_low=lo,
                              macro_auprc_ci_high=hi, macro_auroc=np.mean([r["auroc"] for r in rows]),
                              macro_mcc=np.mean([r["mcc"] for r in rows]) if all(pd.notna(r["mcc"]) for r in rows)
                              else np.nan,
                              macro_brier=np.mean([r["brier"] for r in rows]),
                              macro_delta_vs_null=np.mean([r["delta_vs_null"] for r in rows]),
                              # secondary (post hoc, declared in docs/deviations.md): 12 defined categories,
                              # excluding the heterogeneous "other" donor/acceptor catch-alls
                              macro12_auprc=np.mean([r["auprc"] for r in rows if not r["trait"].endswith("_other")]),
                              macro12_delta_vs_null=np.mean([r["delta_vs_null"] for r in rows
                                                             if not r["trait"].endswith("_other")]),
                              n_traits_better_than_null_q_lt_0_05=int(sum(
                                  r["delta_vs_null_q_bh"] < 0.05 and r["delta_vs_null"] > 0 for r in rows))))
        for a, b in CONTRASTS:
            if a not in P or b not in P:
                continue
            dm = np.nanmean(B[a], axis=1) - np.nanmean(B[b], axis=1)
            lo, hi = ci(dm)
            point = np.mean([safe_ap(Y[obs[:, j], j], P[a][obs[:, j], j]) - safe_ap(Y[obs[:, j], j], P[b][obs[:, j], j])
                             for j in range(len(traits))])
            paired.append(dict(split_level=level, contrast=f"{a} - {b}", trait="MACRO", delta_auprc=point,
                               ci_low=lo, ci_high=hi, p_boot=boot_p(dm), q_bh=np.nan,
                               tied=bool(lo <= 0 <= hi)))
            rows = []
            for j, t in enumerate(traits):
                dj = B[a][:, j] - B[b][:, j]
                lo, hi = ci(dj)
                m = obs[:, j]
                rows.append(dict(split_level=level, contrast=f"{a} - {b}", trait=t,
                                 delta_auprc=safe_ap(Y[m, j], P[a][m, j]) - safe_ap(Y[m, j], P[b][m, j]),
                                 ci_low=lo, ci_high=hi, p_boot=boot_p(dj), tied=bool(lo <= 0 <= hi)))
            for r, qq in zip(rows, bh([r["p_boot"] for r in rows])):
                r["q_bh"] = qq
            paired += rows
        print(f"{level}: done", flush=True)

    pd.DataFrame(per_trait).to_csv(out / "metrics_per_trait.csv", index=False, float_format="%.4f")
    pd.DataFrame(macro).to_csv(out / "metrics_macro.csv", index=False, float_format="%.4f")
    pd.DataFrame(per_fold).to_csv(out / "metrics_per_fold.csv", index=False, float_format="%.4f")
    pd.DataFrame(paired).to_csv(out / "paired_comparisons.csv", index=False, float_format="%.4f")
    cfg = (oof[~oof.arm.isin(NULLS)].drop_duplicates(["arm", "level", "trait", "fold"])
           .groupby(["arm", "level", "config"]).size().rename("n_outer_folds_selected").reset_index())
    cfg.to_csv(out / "selected_configs.csv", index=False)
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain", "src"], cwd=ROOT, capture_output=True,
                                    text=True).stdout.strip())
    except OSError:
        commit, dirty = None, None
    manifest = {"analysis": {"split": args.split, "label_mode": args.label_mode, "subset": args.subset},
                "git_commit": commit, "src_dirty": dirty,
                "code_sha256": {f.name: sha256(f) for f in sorted((ROOT / "src").glob("*.py"))},
                "splits_sha256": sha256(PROC / "splits.tsv"),
                "labels_sha256": sha256(PROC / "labels.tsv"),
                "input_table_sha256": sha256(ROOT / "ko_presence_586genomes.tsv"),
                "ko_matrix_sha256": sha256(PROC / "ko_matrix_panel.tsv"),
                "kofam_ko_list_sha256": sha256(ROOT / "data/raw/kofam/ko_list.gz"),
                "preregistration_sha256": sha256(ROOT / "docs/preregistration.md"),
                "deviations_sha256": sha256(ROOT / "docs/deviations.md"),
                "requirements_lock_sha256": sha256(ROOT / "requirements.lock"),
                "n_bootstrap": N_BOOT, "bootstrap_seed": SEED,
                "oof_files": {f.name: sha256(f) for f in sorted(oof_dir.glob("*.parquet"))}}
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(pd.DataFrame(macro).to_string())


if __name__ == "__main__":
    main()
