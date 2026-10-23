"""AutoML arm with the protocol mismatches fixed (audit of 2026-09-21, finding 4).

Two differences from src/63_automl.py, both of which the audit measured as worth about +0.015 macro AUPRC on a
probe, against a reported deficit of 0.040:
  * FLAML's internal model selection now uses GROUPED splits (the genus/family group of each training genome), so
    the searcher is blocked exactly like every comparator arm, instead of the default StratifiedKFold that puts
    the same genus on both sides.
  * Class weighting is supplied (sample_weight = n_neg/n_pos on positives), matching class_weight="balanced" in
    RF/ET/DT/L1 and scale_pos_weight in B2's XGBoost.
Also replaces the wall-clock time budget with a trial-count budget, so the result does not depend on machine load.

The outer loop is written here rather than reusing nested_cv because the searcher needs the training rows' group
labels, which nested_cv does not pass through.

usage: python src/75_automl_fixed.py --task da [--max-iter 40] [--jobs 5]
Output: results/trees_xai/<task>/oof/AMG.parquet
"""
import argparse
import importlib
import json
import logging
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

tx = importlib.import_module("47_trees_xai")
from redox_common import LEVELS, best_mcc_threshold, inner_folds  # noqa: E402

ESTIMATORS = ["lgbm", "xgboost", "rf", "extra_tree", "lrl2"]


def fit_predict(Xtr, ytr, gtr, Xte, max_iter, seed=0):
    from flaml import AutoML
    from sklearn.model_selection import GroupKFold
    warnings.filterwarnings("ignore")
    logging.getLogger("flaml").setLevel(logging.ERROR)
    n_pos = int(ytr.sum())
    if n_pos < 4 or n_pos > len(ytr) - 4:
        return np.full(len(Xte), ytr.mean())
    w = np.where(ytr == 1, (len(ytr) - n_pos) / max(n_pos, 1), 1.0)
    am = AutoML()
    n_groups = len(np.unique(gtr))
    kw = {}
    if n_groups >= 3:
        kw = dict(split_type=GroupKFold(n_splits=min(3, n_groups)), groups=gtr)
    try:
        am.fit(np.asarray(Xtr), np.asarray(ytr), task="classification", metric="ap", estimator_list=ESTIMATORS,
               max_iter=max_iter, time_budget=-1, eval_method="cv", n_splits=3, sample_weight=w,
               seed=seed, verbose=0, early_stop=True, **kw)
        p = np.asarray(am.predict_proba(np.asarray(Xte)))
        cls = [float(c) for c in am.classes_]
    except Exception as e:
        # FLAML fails on folds where a grouped split leaves one class (seen for donor_nitrite_oxidation, 12
        # positives). Fall back to the class-weighted logistic regression that is already in its estimator list,
        # so the fold stays inside the arm's own model space rather than dropping the trait.
        from sklearn.linear_model import LogisticRegression
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
        print(f"    fallback (logreg) after {type(e).__name__}", flush=True)
        m = make_pipeline(StandardScaler(), LogisticRegression(class_weight="balanced", max_iter=5000))
        m.fit(np.asarray(Xtr), np.asarray(ytr))
        return m.predict_proba(np.asarray(Xte))[:, list(m.classes_).index(1)]
    if p.ndim == 1:                                   # single-column probability of the positive class
        return p
    if p.shape[1] == 1:                               # degenerate learner that saw one class
        return p[:, 0] if cls[0] == 1 else 1.0 - p[:, 0]
    return p[:, cls.index(1.0)] if 1.0 in cls else np.full(len(Xte), float(cls[0] == 1))


def run_trait(Xv, lab, traits, Y, j, level, max_iter):
    fold = lab[f"fold_{level}"].to_numpy()
    groups = lab[f"group_{level}"].to_numpy()
    y = Y[:, j]
    rows = []
    for k in sorted(np.unique(fold)):
        tr, te = np.where(fold != k)[0], np.where(fold == k)[0]
        ytr = y[tr]
        if ytr.sum() == 0 or ytr.sum() == len(ytr):
            p_te, thr = np.full(len(te), ytr.mean()), 0.5
        else:
            folds = inner_folds(ytr, groups[tr])
            oof = np.full(len(tr), np.nan)
            for itr, iva in folds:                     # same 4-fold grouped inner OOF the other arms use
                yi = ytr[itr]
                oof[iva] = yi.mean() if yi.sum() in (0, len(yi)) else \
                    fit_predict(Xv[tr][itr], yi, groups[tr][itr], Xv[tr][iva], max_iter)
            thr = best_mcc_threshold(ytr, oof)
            p_te = fit_predict(Xv[tr], ytr, groups[tr], Xv[te], max_iter)
        for i, p in zip(te, p_te):
            rows.append(("AMG", level, traits[j], lab.genome_id[i], int(k), y[i], float(p), float(thr),
                         json.dumps({"model": "flaml_grouped", "max_iter": max_iter}), np.nan))
    print(f"[AMG/{level}] {traits[j]} done", flush=True)
    return pd.DataFrame(rows, columns=["arm", "level", "trait", "genome_id", "fold", "y", "prob", "threshold",
                                       "config", "inner_auprc"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=["da", "pw"], default="da")
    ap.add_argument("--max-iter", type=int, default=40)
    ap.add_argument("--jobs", type=int, default=5)
    args = ap.parse_args()
    lab, traits, Y, X = tx.load_task(args.task)
    Xv = X.to_numpy()
    print(f"{args.task}: {Xv.shape[0]} genomes, {len(traits)} traits, {Xv.shape[1]} features, "
          f"{args.max_iter} trials per fit", flush=True)
    def guarded(j, lvl):
        try:
            return run_trait(Xv, lab, traits, Y, j, lvl, args.max_iter)
        except Exception as e:                    # one estimator failing must not lose the whole arm
            print(f"[AMG/{lvl}] {traits[j]} FAILED: {type(e).__name__}: {e}", flush=True)
            return None
    out = [d for d in Parallel(n_jobs=args.jobs)(delayed(guarded)(j, lvl)
                                                 for lvl in LEVELS for j in range(len(traits))) if d is not None]
    df = pd.concat(out, ignore_index=True)
    d = tx.OUT / args.task / "oof"
    d.mkdir(parents=True, exist_ok=True)
    df.to_parquet(d / "AMG.parquet", index=False)
    print(d / "AMG.parquet", df.groupby(["arm", "level"]).size().to_dict())


if __name__ == "__main__":
    main()
