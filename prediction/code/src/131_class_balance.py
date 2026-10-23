"""Pre-registration CB (docs/preregistration_CB_class_balance.md): other treatments of class imbalance for RedoxKO.

Arms CB_none (no weighting), CB_under (negatives undersampled to 1:1, 5 draws averaged), CB_over (positives
oversampled to 1:1), each with the B2 features, configurations, nested CV and threshold rule. Compared with the
existing B2 out-of-fold predictions (inverse-frequency class weights).

usage: python src/131_class_balance.py [--jobs N]
Output: results/class_balance/oof/<arm>.parquet, results/class_balance/summary.csv
"""
import argparse
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from redox_common import (LEVELS, RES, load_features, load_labelled, nested_cv, oof_matrix, paired_macro_delta,
                          safe_ap)

OUT = RES / "class_balance"
LR_C = [0.01, 0.1, 1.0, 10.0]
CONFIGS = [{"model": "logreg", "C": c} for c in LR_C] + [{"model": "xgb", "max_depth": d} for d in (2, 4)]
ARMS = ("CB_none", "CB_under", "CB_over")


def fit_unweighted(cfg, X, y):
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from xgboost import XGBClassifier
    if cfg["model"] == "logreg":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            return [make_pipeline(StandardScaler(), LogisticRegression(C=cfg["C"], max_iter=5000)).fit(X, y)]
    return [XGBClassifier(n_estimators=300, learning_rate=0.05, max_depth=cfg["max_depth"], subsample=0.8,
                          colsample_bytree=0.5, min_child_weight=1, scale_pos_weight=1.0, n_jobs=1,
                          random_state=1000 + s, verbosity=0).fit(X, y) for s in range(3)]


def predict(models, X):
    return np.mean([m.predict_proba(X)[:, 1] for m in models], axis=0)


def make_fit_predict(arm):
    def fp(cfg, Xtr, ytr, Xte, seed=0):
        pos, neg = np.where(ytr == 1)[0], np.where(ytr == 0)[0]
        if arm == "CB_none":
            return predict(fit_unweighted(cfg, Xtr, ytr), Xte)
        rng = np.random.default_rng(20261020 + seed)
        if arm == "CB_under":
            out = []
            for _ in range(5):
                idx = np.concatenate([pos, rng.choice(neg, size=min(len(pos), len(neg)), replace=False)])
                out.append(predict(fit_unweighted(cfg, Xtr[idx], ytr[idx]), Xte))
            return np.mean(out, axis=0)
        extra = rng.choice(pos, size=max(len(neg) - len(pos), 0), replace=True)
        idx = np.concatenate([np.arange(len(ytr)), extra])
        return predict(fit_unweighted(cfg, Xtr[idx], ytr[idx]), Xte)
    return fp


def macro_mcc(oof):
    from sklearn.metrics import matthews_corrcoef
    return float(np.mean([matthews_corrcoef(g.y, (g.prob >= g.threshold).astype(int)) for _, g in oof.groupby("trait")]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()
    lab, traits, Y = load_labelled()
    X = load_features("energy_panel", lab.genome_id).to_numpy()
    X = X[:, X.std(axis=0) > 0]
    (OUT / "oof").mkdir(parents=True, exist_ok=True)
    for arm in ARMS:
        f = OUT / "oof" / f"{arm}.parquet"
        if f.exists():
            continue
        fp = make_fit_predict(arm)

        def one(level, j):
            m = ~np.isnan(Y[:, j])
            return nested_cv(arm, X[m], lab.loc[m].reset_index(drop=True), [traits[j]], Y[m][:, [j]], CONFIGS, fp,
                             level, log=lambda s: None)

        out = Parallel(n_jobs=args.jobs)(delayed(one)(lvl, j) for lvl in LEVELS for j in range(len(traits)))
        pd.concat(out, ignore_index=True).to_parquet(f)
        print(arm, "done", flush=True)
    prev = np.nanmean(Y, axis=0)
    rare = [j for j in range(len(traits)) if prev[j] < 0.10]
    ref_oof = pd.read_parquet(RES / "oof" / "B2.parquet")
    rows = []
    for level in LEVELS:
        Pb = oof_matrix(RES / "oof", "B2", level, lab, traits)
        for arm in ("B2",) + ARMS:
            d = OUT / "oof" if arm != "B2" else RES / "oof"
            Pa = oof_matrix(d, arm, level, lab, traits)
            oof = (ref_oof if arm == "B2" else pd.read_parquet(d / f"{arm}.parquet"))
            oof = oof[oof.level == level]
            obs = ~np.isnan(Y)
            ap_all = np.nanmean([safe_ap(Y[obs[:, j], j], Pa[obs[:, j], j]) for j in range(len(traits))])
            ap_rare = np.nanmean([safe_ap(Y[obs[:, j], j], Pa[obs[:, j], j]) for j in rare])
            row = dict(arm=arm, level=level, macro_auprc=ap_all, macro_auprc_rare=ap_rare, macro_mcc=macro_mcc(oof))
            if arm != "B2":
                pt, lo, hi = paired_macro_delta(lab, traits, Y, Pa, Pb, level)
                row.update(delta_vs_B2=pt, ci_low=lo, ci_high=hi)
            rows.append(row)
    s = pd.DataFrame(rows)
    s.to_csv(OUT / "summary.csv", index=False, float_format="%.4f")
    print(f"rare categories (<10% positives): {[traits[j] for j in rare]}")
    print(s.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
