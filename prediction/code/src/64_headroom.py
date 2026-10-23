"""Headroom diagnostics (docs/preregistration_trees_xai.md, Addendum P): what limits the metrics - capacity,
lineage transfer, or the labels? Diagnostic only; no released model comes from this file.

H1 capacity      : in-sample fit on all genomes (train = test). If this is ~1.0 the KO feature space can express the
                   labels, so model capacity is not the binding constraint.
H3 label noise   : hide a fraction r of the positives in the TRAINING folds (simulating traits that were never
                   reported), retrain, and score against the true labels -> the cost of false negatives in training.
                   Separately score the clean model against corrupted TEST labels -> how far unreported positives
                   push the MEASURED metric below the true skill.

usage: python src/64_headroom.py --task da [--rates 0.1 0.2 0.3] [--seeds 3] [--jobs 5]
Output: results/trees_xai/<task>/headroom_capacity.csv, headroom_label_noise.csv
"""
import argparse
import importlib

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.ensemble import RandomForestClassifier

tx = importlib.import_module("47_trees_xai")
from redox_common import safe_ap  # noqa: E402


def rf(seed=0):
    return RandomForestClassifier(n_estimators=300, min_samples_leaf=1, class_weight="balanced_subsample",
                                  n_jobs=2, random_state=seed)


def capacity(X, Y, traits):
    """In-sample AUPRC (train = test): an upper bound on what this feature space can represent at all."""
    rows = []
    for j, t in enumerate(traits):
        y = Y[:, j]
        if not 0 < y.sum() < len(y):
            continue
        p = rf().fit(X, y).predict_proba(X)[:, 1]
        rows.append(dict(trait=t, n_pos=int(y.sum()), in_sample_auprc=safe_ap(y, p)))
    return pd.DataFrame(rows)


def one_noise(X, Y, traits, fold, j, r, seed):
    """Fold-wise: train on labels with a fraction r of training positives hidden; score against true test labels."""
    rng = np.random.default_rng(1000 * seed + j)
    y = Y[:, j]
    out = []
    for k in sorted(np.unique(fold)):
        tr, te = np.where(fold != k)[0], np.where(fold == k)[0]
        ytr = y[tr].copy()
        pos = np.where(ytr == 1)[0]
        if len(pos) < 2 or y[te].sum() == 0:
            continue
        hide = rng.choice(pos, size=int(round(r * len(pos))), replace=False)
        ytr_noisy = ytr.copy()
        ytr_noisy[hide] = 0
        if ytr_noisy.sum() == 0:
            continue
        p_noisy = rf(seed).fit(X[tr], ytr_noisy).predict_proba(X[te])[:, 1]
        p_clean = rf(seed).fit(X[tr], ytr).predict_proba(X[te])[:, 1]
        # corrupted test labels: the same hiding rate applied to the test fold's positives
        yte_noisy = y[te].copy()
        tpos = np.where(yte_noisy == 1)[0]
        yte_noisy[rng.choice(tpos, size=int(round(r * len(tpos))), replace=False)] = 0
        out.append(dict(trait=traits[j], rate=r, seed=seed, fold=int(k), n_pos_test=int(y[te].sum()),
                        auprc_clean_vs_true=safe_ap(y[te], p_clean),
                        auprc_noisytrain_vs_true=safe_ap(y[te], p_noisy),
                        auprc_clean_vs_noisytest=safe_ap(yte_noisy, p_clean) if yte_noisy.sum() else np.nan))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=["da", "pw"], default="da")
    ap.add_argument("--rates", type=float, nargs="+", default=[0.1, 0.2, 0.3])
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--jobs", type=int, default=5)
    args = ap.parse_args()
    lab, traits, Y, X = tx.load_task(args.task)
    Xv = X.to_numpy()
    out = tx.OUT / args.task
    print(f"{args.task}: {Xv.shape[0]} genomes, {len(traits)} traits, {Xv.shape[1]} features", flush=True)

    cap = capacity(Xv, Y, traits)
    cap.to_csv(out / "headroom_capacity.csv", index=False, float_format="%.4f")
    print(f"H1 in-sample AUPRC: mean {cap.in_sample_auprc.mean():.3f}, min {cap.in_sample_auprc.min():.3f}", flush=True)

    fold = lab["fold_genus"].to_numpy()
    res = Parallel(n_jobs=args.jobs)(delayed(one_noise)(Xv, Y, traits, fold, j, r, s)
                                     for r in args.rates for s in range(args.seeds) for j in range(len(traits)))
    d = pd.DataFrame([x for sub in res for x in sub])
    d.to_csv(out / "headroom_label_noise.csv", index=False, float_format="%.4f")
    g = d.groupby("rate")[["auprc_clean_vs_true", "auprc_noisytrain_vs_true", "auprc_clean_vs_noisytest"]].mean()
    print("H3 (genus folds, macro over traits and seeds):")
    print(g.to_string())


if __name__ == "__main__":
    main()
