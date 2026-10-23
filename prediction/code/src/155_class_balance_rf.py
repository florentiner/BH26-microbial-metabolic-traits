"""Class balancing for the final model RedoxKO (random forest on the 852-KO panel): the comparison of
docs/preregistration_CB_class_balance.md (run there on RedoxKO-LX) repeated for the forest. Not registered.

Arms: RFCB_none (no class weights), RFCB_under (negatives undersampled to 1:1, 5 draws averaged), RFCB_over
(positives oversampled to 1:1), each with the forest's configurations and nested CV; compared with the forest's own
out-of-fold predictions (class weights within each bootstrap sample, results/trees_xai/da/oof/RF.parquet).
usage: python src/155_class_balance_rf.py [--jobs N]
Output: results/class_balance_rf/oof/<arm>.parquet, summary.csv, by_prevalence.csv
"""
import argparse
import importlib
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import LEVELS, RES, load_features, load_labelled, nested_cv, oof_matrix, paired_macro_delta, safe_ap  # noqa

tx = importlib.import_module("47_trees_xai")
tx.TREE_JOBS = 1
OUT = RES / "class_balance_rf"
REF_DIR, REF = RES / "trees_xai/da/oof", "RF"
CONFIGS = tx.CONFIGS["rf"]
ARMS = ("RFCB_none", "RFCB_under", "RFCB_over")


def forest(cfg, X, y):
    from sklearn.ensemble import RandomForestClassifier
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return [RandomForestClassifier(n_estimators=500, max_features="sqrt", min_samples_leaf=cfg["msl"],
                                       class_weight=None, n_jobs=1, random_state=s).fit(X, y) for s in tx.SEEDS]


def predict(models, X):
    return np.mean([m.predict_proba(X)[:, 1] for m in models], axis=0)


def make_fit_predict(arm):
    def fp(cfg, Xtr, ytr, Xte, seed=0):
        pos, neg = np.where(ytr == 1)[0], np.where(ytr == 0)[0]
        if arm == "RFCB_none":
            return predict(forest(cfg, Xtr, ytr), Xte)
        rng = np.random.default_rng(20261021 + seed)
        if arm == "RFCB_under":
            out = []
            for _ in range(5):
                idx = np.concatenate([pos, rng.choice(neg, size=min(len(pos), len(neg)), replace=False)])
                out.append(predict(forest(cfg, Xtr[idx], ytr[idx]), Xte))
            return np.mean(out, axis=0)
        extra = rng.choice(pos, size=max(len(neg) - len(pos), 0), replace=True)
        idx = np.concatenate([np.arange(len(ytr)), extra])
        return predict(forest(cfg, Xtr[idx], ytr[idx]), Xte)
    return fp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()
    lab, traits, Y = load_labelled()
    X = load_features("energy_panel", lab.genome_id)
    X = X.loc[:, X.std(axis=0) > 0].to_numpy(dtype=np.float32)
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
    groups = {"< 10% positive": prev < 0.10, "10-30% positive": (prev >= 0.10) & (prev <= 0.30), "> 30% positive": prev > 0.30}
    obs = ~np.isnan(Y)
    rows, grows = [], []
    for level in LEVELS:
        Pb = oof_matrix(REF_DIR, REF, level, lab, traits)
        apb = np.array([safe_ap(Y[obs[:, j], j], Pb[obs[:, j], j]) for j in range(len(traits))])
        rows.append(dict(arm=REF, level=level, macro_auprc=float(np.nanmean(apb))))
        for g, m in groups.items():
            grows.append(dict(level=level, group=g, n=int(m.sum()), arm=REF, auprc=float(np.nanmean(apb[m]))))
        for arm in ARMS:
            Pa = oof_matrix(OUT / "oof", arm, level, lab, traits)
            apa = np.array([safe_ap(Y[obs[:, j], j], Pa[obs[:, j], j]) for j in range(len(traits))])
            pt, lo, hi = paired_macro_delta(lab, traits, Y, Pa, Pb, level)
            rows.append(dict(arm=arm, level=level, macro_auprc=float(np.nanmean(apa)), delta_vs_RF=pt, ci_low=lo, ci_high=hi,
                             n_lower=int((apa < apb).sum())))
            for g, m in groups.items():
                grows.append(dict(level=level, group=g, n=int(m.sum()), arm=arm, auprc=float(np.nanmean(apa[m])),
                                  delta=float(np.nanmean(apa[m] - apb[m])), n_lower=int((apa[m] < apb[m]).sum())))
    s = pd.DataFrame(rows)
    s.to_csv(OUT / "summary.csv", index=False, float_format="%.4f")
    pd.DataFrame(grows).to_csv(OUT / "by_prevalence.csv", index=False, float_format="%.4f")
    print(s.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(pd.DataFrame(grows).to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
