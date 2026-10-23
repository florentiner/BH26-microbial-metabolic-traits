"""Calibration of the released model and a learning curve in training genera (audit round 3, missing analyses 4 and 5).

calibration: reliability of the out-of-fold probabilities (held-out genus and family) and on the external set:
  Brier score, calibration slope/intercept from a logistic recalibration of the logit, and 5-bin reliability tables.
learning curve: macro AUPRC as a function of the number of training genera (25/50/75/100 % of the genera in each
  outer training fold, 3 draws per size, same frozen outer test folds), for the released arm and for the best
  fusion arm, to put a number on "not distinguishable at N = 358".
usage: python src/26_calibration_learning_curve.py [calibration|curve|both] [--arms B2 B2XLE_large]
Output: results/calibration_per_trait.csv, results/calibration_reliability_bins.csv, results/learning_curve.csv
"""
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from redox_common import LEVELS, RES, fit_models, inner_folds, load_features, load_labelled, oof_matrix, safe_ap

EPS = 1e-6
SIZES = (0.25, 0.5, 0.75, 1.0)
DRAWS = 3
SEED = 20260918


def calibration():
    lab, traits, Y = load_labelled()
    rows, bins = [], []
    for level in LEVELS:
        for arm in ("B2", "B2XL", "B2XLE_large"):
            f = RES / "oof" / f"{arm}.parquet"
            if not f.exists():
                continue
            P = oof_matrix(RES / "oof", arm, level, lab, traits)
            for j, t in enumerate(traits):
                y, p = Y[:, j], np.clip(P[:, j], EPS, 1 - EPS)
                lr = LogisticRegression(C=1e6, solver="lbfgs").fit(np.log(p / (1 - p)).reshape(-1, 1), y)
                rows.append(dict(source="internal_cv", split_level=level, arm=arm, trait=t, n=len(y),
                                 prevalence=y.mean(), mean_predicted=p.mean(), brier=float(np.mean((p - y) ** 2)),
                                 calibration_slope=float(lr.coef_[0][0]), calibration_intercept=float(lr.intercept_[0])))
                q = pd.qcut(p, 5, duplicates="drop", labels=False)
                for b in np.unique(q):
                    m = q == b
                    bins.append(dict(source="internal_cv", split_level=level, arm=arm, trait=t, bin=int(b),
                                     n=int(m.sum()), mean_predicted=float(p[m].mean()), observed=float(y[m].mean())))
    ext = RES / "external_validation_predictions.csv"
    if ext.exists():
        e = pd.read_csv(ext)
        for (arm, t), d in e.groupby(["arm", "trait"]):
            if arm.startswith("null"):
                continue
            y, p = d.y.to_numpy(), np.clip(d.score.to_numpy(), EPS, 1 - EPS)
            if not (0 < y.sum() < len(y)):
                continue
            lr = LogisticRegression(C=1e6, solver="lbfgs").fit(np.log(p / (1 - p)).reshape(-1, 1), y)
            rows.append(dict(source="external", split_level="novel_genus", arm=arm, trait=t, n=len(y),
                             prevalence=y.mean(), mean_predicted=p.mean(), brier=float(np.mean((p - y) ** 2)),
                             calibration_slope=float(lr.coef_[0][0]), calibration_intercept=float(lr.intercept_[0])))
    pd.DataFrame(rows).to_csv(RES / "calibration_per_trait.csv", index=False, float_format="%.4f")
    pd.DataFrame(bins).to_csv(RES / "calibration_reliability_bins.csv", index=False, float_format="%.4f")
    r = pd.DataFrame(rows)
    print(r.groupby(["source", "split_level", "arm"])[["brier", "calibration_slope", "calibration_intercept"]]
          .mean().round(3).to_string())


def curve(arms):
    lab, traits, Y = load_labelled()
    b2 = load_features("energy_panel", lab.genome_id)
    b2 = b2.loc[:, b2.std(axis=0) > 0].to_numpy()
    emb = load_features("bacformer_large_mean", lab.genome_id).to_numpy() if "B2XLE_large" in arms else None
    rng = np.random.default_rng(SEED)
    rows = []
    for level in ("genus",):
        fold = lab[f"fold_{level}"].to_numpy()
        groups = lab[f"group_{level}"].to_numpy()
        for frac in SIZES:
            for draw in range(DRAWS if frac < 1.0 else 1):
                P = {a: np.full(Y.shape, np.nan) for a in arms}
                n_gen = []
                for k in np.unique(fold):
                    tr_all, te = np.where(fold != k)[0], np.where(fold == k)[0]
                    g = np.unique(groups[tr_all])
                    keep = rng.choice(g, max(2, int(round(frac * len(g)))), replace=False)
                    tr = tr_all[np.isin(groups[tr_all], keep)]
                    n_gen.append(len(keep))
                    for j in range(len(traits)):
                        y = Y[tr, j]
                        if y.sum() < 2 or y.sum() == len(y):
                            for a in arms:
                                P[a][te, j] = y.mean()
                            continue
                        folds = inner_folds(y, groups[tr])
                        best, best_ap = None, -1
                        for cfg in [{"model": "logreg", "C": c} for c in (0.01, 0.1, 1.0, 10.0)] + \
                                   [{"model": "xgb", "max_depth": d} for d in (2, 4)]:
                            oof = np.zeros(len(tr))
                            for itr, iva in folds:
                                yi = y[itr]
                                oof[iva] = yi.mean() if yi.sum() in (0, len(yi)) else np.mean(
                                    [m.predict_proba(b2[tr][iva])[:, 1] for m in fit_models(cfg, b2[tr][itr], yi)],
                                    axis=0)
                            ap = safe_ap(y, oof)
                            if ap > best_ap:
                                best, best_ap = cfg, ap
                        ko = np.mean([m.predict_proba(b2[te])[:, 1] for m in fit_models(best, b2[tr], y)], axis=0)
                        if "B2" in arms:
                            P["B2"][te, j] = ko
                        if "B2XLE_large" in arms:
                            lrc = {"model": "logreg", "C": best.get("C", 0.1)}
                            xgc = {"model": "xgb", "max_depth": best.get("max_depth", 2)}
                            kox = 0.5 * (np.mean([m.predict_proba(b2[te])[:, 1] for m in fit_models(lrc, b2[tr], y)], axis=0) +
                                         np.mean([m.predict_proba(b2[te])[:, 1] for m in fit_models(xgc, b2[tr], y)], axis=0))
                            pe = np.mean([m.predict_proba(emb[te])[:, 1] for m in fit_models(lrc, emb[tr], y)], axis=0)
                            P["B2XLE_large"][te, j] = 0.5 * (kox + pe)
                for a in arms:
                    rows.append(dict(arm=a, split_level=level, frac_genera=frac, draw=draw,
                                     mean_train_genera=float(np.mean(n_gen)),
                                     macro_auprc=float(np.nanmean([safe_ap(Y[:, j], P[a][:, j])
                                                                   for j in range(len(traits))]))))
                print(rows[-1], flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(RES / "learning_curve.csv", index=False, float_format="%.4f")
    print(out.groupby(["arm", "frac_genera"]).macro_auprc.agg(["mean", "min", "max"]).round(4).to_string())


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "both"
    arms = sys.argv[sys.argv.index("--arms") + 1:] if "--arms" in sys.argv else ["B2", "B2XLE_large"]
    if what in ("calibration", "both"):
        calibration()
    if what in ("curve", "both"):
        curve(arms)
