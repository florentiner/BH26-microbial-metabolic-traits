"""AutoML arm (docs/preregistration_trees_xai.md, Addendum O): FLAML searches the model space per trait and fold.

Everything outside the searcher is identical to the other arms: frozen genus/family folds, inner grouped CV for the
MCC cut-off, the same KO feature matrix as B2/RF/ET (da) and PB2 (pw). The searcher runs inside the outer training
fold only, so no test genome is ever seen during the search.

usage: python src/63_automl.py --task {da,pw} [--jobs N] [--budget 30] [--inner-budget 15]
Output: results/trees_xai/<task>/oof/AM.parquet   (47_trees_xai.py metrics --task <task> picks it up)
"""
import argparse
import importlib
import logging
import os
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

tx = importlib.import_module("47_trees_xai")
from redox_common import LEVELS, nested_cv  # noqa: E402

ESTIMATORS = ["lgbm", "xgboost", "rf", "extra_tree", "lrl2"]
BUDGET = 30   # seconds per FLAML search, inner and outer alike


def fit_predict(cfg, Xtr, ytr, Xte, seed=0):
    """One FLAML search. cfg carries only the time budget; the searcher owns the model and its hyper-parameters."""
    from flaml import AutoML
    warnings.filterwarnings("ignore")
    logging.getLogger("flaml").setLevel(logging.ERROR)
    n_pos = int(ytr.sum())
    # FLAML needs at least one positive per internal split; with very few positives fall back to the prevalence.
    if n_pos < 4 or n_pos > len(ytr) - 4:
        return np.full(len(Xte), ytr.mean())
    am = AutoML()
    am.fit(np.asarray(Xtr), np.asarray(ytr), task="classification", metric="ap",
           estimator_list=ESTIMATORS, time_budget=cfg["budget"], eval_method="cv", n_splits=3,
           seed=seed, verbose=0, early_stop=True)
    p = am.predict_proba(np.asarray(Xte))
    cls = list(am.classes_)
    return p[:, cls.index(1)] if 1 in cls else np.full(len(Xte), float(cls[0] == 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=["da", "pw"], required=True)
    ap.add_argument("--jobs", type=int, default=5)
    ap.add_argument("--budget", type=int, default=BUDGET)
    args = ap.parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    lab, traits, Y, X = tx.load_task(args.task)
    Xv = X.to_numpy()
    print(f"{args.task}: {len(lab)} genomes, {len(traits)} traits, feature matrix {Xv.shape}", flush=True)

    def one(j, level):
        # a single "config": the budget. nested_cv still runs the inner folds, which set the MCC cut-off.
        cfgs = [{"model": "flaml", "budget": args.budget}]
        return nested_cv("AM" + tx.SUFFIX, Xv, lab, [traits[j]], Y[:, [j]], cfgs, fit_predict, level,
                         log=lambda *a, **k: None)

    out = Parallel(n_jobs=args.jobs)(delayed(one)(j, lvl) for lvl in LEVELS for j in range(len(traits)))
    df = pd.concat(out, ignore_index=True)
    outdir = tx.oof_dir(args.task)
    outdir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(outdir / f"AM{tx.SUFFIX}.parquet", index=False)
    print(outdir / f"AM{tx.SUFFIX}.parquet", df.groupby(["arm", "level"]).size().to_dict())


if __name__ == "__main__":
    main()
