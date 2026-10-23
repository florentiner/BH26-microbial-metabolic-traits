"""Fold-allocation sensitivity (Supplementary Table S4), random split and learning curve (Fig. 5) for the final model
RedoxKO (random forest on the 852-KO panel), with the same splits, nested CV and genus subsampling as for the former
final model (src/06_run_arms.py --split, src/26_calibration_learning_curve.py curve).
usage: python src/148_rf_sensitivity_curve.py [splits|curve|both]
Output: results/sensitivity/<split>__main/oof/RF.parquet, results/sensitivity_rf.csv, results/learning_curve_rf.csv
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import RES, inner_folds, load_features, load_labelled, nested_cv, safe_ap  # noqa: E402

tx = importlib.import_module("47_trees_xai")
tx.TREE_JOBS = 2
CFG = tx.CONFIGS["rf"]
SPLITS = ["seed1", "seed2", "seed3", "seed4", "seed5", "gtdb", "random"]


def fp(cfg, Xtr, ytr, Xte, seed=0):
    return tx.fit_predict(cfg, Xtr, ytr, Xte)


def splits():
    rows = []
    for split in SPLITS:
        lab, traits, Y = load_labelled(split=split)
        X = load_features("energy_panel", lab.genome_id)
        X = X.loc[:, X.std(axis=0) > 0].to_numpy(dtype=np.float32)
        out = RES / f"sensitivity/{split}__main/oof"
        out.mkdir(parents=True, exist_ok=True)
        levels = ("genus",) if split == "random" else ("genus", "family")

        def one(level, j):
            m = ~np.isnan(Y[:, j])
            return nested_cv("RF", X[m], lab.loc[m].reset_index(drop=True), [traits[j]], Y[m][:, [j]], CFG, fp, level,
                             log=lambda s: None)
        oof = pd.concat(Parallel(n_jobs=5)(delayed(one)(lv, j) for lv in levels for j in range(len(traits))))
        oof.to_parquet(out / "RF.parquet")
        for lv in levels:
            o = oof[oof.level == lv]
            ap = [safe_ap(g.y.to_numpy(), g.prob.to_numpy()) for _, g in o.groupby("trait")]
            rows.append(dict(split=split, level=lv, macro_auprc=float(np.nanmean(ap))))
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(RES / "sensitivity_rf.csv", index=False, float_format="%.4f")


def curve():
    c26 = importlib.import_module("26_calibration_learning_curve")
    lab, traits, Y = load_labelled()
    X = load_features("energy_panel", lab.genome_id)
    X = X.loc[:, X.std(axis=0) > 0].to_numpy(dtype=np.float32)
    rng = np.random.default_rng(c26.SEED)
    fold, groups = lab.fold_genus.to_numpy(), lab.group_genus.to_numpy()
    plans = []
    for frac in c26.SIZES:
        for draw in range(c26.DRAWS if frac < 1.0 else 1):
            parts = []
            for k in np.unique(fold):
                tr_all, te = np.where(fold != k)[0], np.where(fold == k)[0]
                g = np.unique(groups[tr_all])
                keep = rng.choice(g, max(2, int(round(frac * len(g)))), replace=False)
                parts.append((tr_all[np.isin(groups[tr_all], keep)], te, len(keep)))
            plans.append((frac, draw, parts))

    def run(frac, draw, parts, j):
        P = np.full(len(lab), np.nan)
        for tr, te, _ in parts:
            y = Y[tr, j]
            ok = ~np.isnan(y)
            tr, y = tr[ok], y[ok]
            if y.sum() < 2 or y.sum() == len(y):
                P[te] = y.mean()
                continue
            folds = inner_folds(y, groups[tr])
            best, best_ap = None, -1
            for cfg in CFG:
                oof = np.zeros(len(tr))
                for itr, iva in folds:
                    yi = y[itr]
                    oof[iva] = yi.mean() if yi.sum() in (0, len(yi)) else fp(cfg, X[tr][itr], yi, X[tr][iva])
                ap = safe_ap(y, oof)
                if ap > best_ap:
                    best, best_ap = cfg, ap
            P[te] = fp(best, X[tr], y, X[te])
        m = ~np.isnan(Y[:, j])
        return frac, draw, j, safe_ap(Y[m, j], P[m])
    res = Parallel(n_jobs=5)(delayed(run)(f, d, p, j) for f, d, p in plans for j in range(len(traits)))
    r = pd.DataFrame(res, columns=["frac_genera", "draw", "j", "auprc"])
    n_gen = {(f, d): float(np.mean([p[2] for p in parts])) for f, d, parts in plans}
    out = r.groupby(["frac_genera", "draw"]).auprc.mean().reset_index().rename(columns={"auprc": "macro_auprc"})
    out["mean_train_genera"] = [n_gen[(f, d)] for f, d in zip(out.frac_genera, out.draw)]
    out.insert(0, "arm", "RF")
    out.insert(1, "split_level", "genus")
    out.to_csv(RES / "learning_curve_rf.csv", index=False, float_format="%.4f")
    print(out.groupby("frac_genera").macro_auprc.agg(["mean", "min", "max"]).round(4).to_string())


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "both"
    if what in ("splits", "both"):
        splits()
    if what in ("curve", "both"):
        curve()
