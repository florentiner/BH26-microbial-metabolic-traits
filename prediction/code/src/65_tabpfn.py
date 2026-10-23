"""Tabular foundation model arm (docs/preregistration_trees_xai.md, Addendum Q).

TabPFN v2 is a transformer pre-trained on synthetic tabular tasks that predicts in context: it is "fitted" by being
shown the training rows, with no hyper-parameter search. Same KO features and the same frozen folds as every other
arm; the per-trait MCC cut-off comes from one inner grouped split of the outer training fold.

usage: SCIPY_ARRAY_API=1 python src/65_tabpfn.py --task da [--device mps] [--n-estimators 4]
Output: results/trees_xai/<task>/oof/TPFN.parquet
"""
import argparse
import importlib
import json
import time

import numpy as np
import pandas as pd

tx = importlib.import_module("47_trees_xai")
from redox_common import LEVELS, best_mcc_threshold, inner_folds  # noqa: E402


def fit_predict(X_tr, y_tr, X_te, device, n_est, seed=0):
    from tabpfn import TabPFNClassifier
    m = TabPFNClassifier(device=device, n_estimators=n_est, ignore_pretraining_limits=True, random_state=seed)
    m.fit(X_tr, y_tr)
    return m.predict_proba(X_te)[:, 1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=["da", "pw"], default="da")
    ap.add_argument("--device", default="mps")
    ap.add_argument("--n-estimators", type=int, default=4)
    args = ap.parse_args()
    lab, traits, Y, X = tx.load_task(args.task)
    Xv = X.to_numpy().astype(np.float32)
    print(f"{args.task}: {Xv.shape[0]} genomes, {len(traits)} traits, {Xv.shape[1]} features", flush=True)
    cfg = json.dumps({"model": "tabpfn2", "n_estimators": args.n_estimators, "device": args.device})
    rows, t0 = [], time.time()
    for level in LEVELS:
        fold = lab[f"fold_{level}"].to_numpy()
        groups = lab[f"group_{level}"].to_numpy()
        for j, trait in enumerate(traits):
            y = Y[:, j]
            for k in sorted(np.unique(fold)):
                tr, te = np.where(fold != k)[0], np.where(fold == k)[0]
                ytr = y[tr]
                if ytr.sum() == 0 or ytr.sum() == len(ytr):        # degenerate training fold
                    p_te, thr = np.full(len(te), ytr.mean()), 0.5
                else:
                    itr, iva = inner_folds(ytr, groups[tr])[0]      # one inner split, for the cut-off only
                    if 0 < ytr[itr].sum() < len(itr):
                        p_in = fit_predict(Xv[tr][itr], ytr[itr], Xv[tr][iva], args.device, args.n_estimators)
                        thr = best_mcc_threshold(ytr[iva], p_in)
                    else:
                        thr = 0.5
                    p_te = fit_predict(Xv[tr], ytr, Xv[te], args.device, args.n_estimators)
                for i, p in zip(te, p_te):
                    rows.append(("TPFN" + tx.SUFFIX, level, trait, lab.genome_id[i], int(k), y[i], float(p), float(thr),
                                 cfg, np.nan))
            print(f"[TPFN/{level}] {trait} done ({time.time() - t0:.0f}s)", flush=True)
    df = pd.DataFrame(rows, columns=["arm", "level", "trait", "genome_id", "fold", "y", "prob", "threshold",
                                     "config", "inner_auprc"])
    out = tx.oof_dir(args.task)
    out.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out / f"TPFN{tx.SUFFIX}.parquet", index=False)
    print(out / f"TPFN{tx.SUFFIX}.parquet", df.groupby(["arm", "level"]).size().to_dict())


if __name__ == "__main__":
    main()
