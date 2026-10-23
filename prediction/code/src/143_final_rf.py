"""Final donor and acceptor model RedoxKO: random forest on the 852-KO panel (decision of 21 October 2026, after the
random forest did better than the former final model, now RedoxKO-LX, on every independent test; Supplementary Note 1).

Per category, as for every other final model (src/09_train_final.py): min_samples_leaf (1 or 5) chosen by grouped
(genus) inner CV AUPRC on all labeled genomes, decision threshold = best MCC on the inner out-of-fold scores, Platt
calibrator on the same scores, refit on all 358 labeled genomes (500 trees, balanced subsample, 3 seeds).
Yes/no calls are released per category by the existing rule (src/71_score_validation.reliable_traits): at family level
MCC lower bound > 0.2 and better than the taxonomic baseline at BH q < 0.05, and not useless on the external genomes
(AUPRC lift > 1); computed from the forest's own outer-CV metrics.
Output: results/models/RF/model.joblib, metadata.json; results/models/RELEASED.json (released_arm = RF)
"""
import importlib
import json
import sys
from datetime import date
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import (PROC, RES, FittedModel, best_mcc_threshold, fit_calibrator, inner_folds, load_features,  # noqa
                          load_labelled, safe_ap, sha256)

tx = importlib.import_module("47_trees_xai")
tx.TREE_JOBS = 8


def reliable(traits):
    mt = pd.read_csv(RES / "trees_xai/da/metrics_per_trait.csv")
    mt = mt[mt.arm == "RF"]
    ev = pd.read_csv(RES / "external_validation_missing_rows_predictions.csv")
    ev = ev[ev.arm == "Random forest"]
    ext = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t").set_index("genome_id")
    lift_ok = {}
    for t, g in ev.groupby("trait"):
        if t not in ext.columns:
            continue
        y = ext[t].reindex(g.genome_id).to_numpy(dtype=float)
        m = ~np.isnan(y)
        if y[m].sum() >= 10 and (y[m] == 0).sum() >= 10:
            lift_ok[t] = safe_ap(y[m], g.score.to_numpy()[m]) / y[m].mean() > 1
    out = {}
    for level in ("genus", "family"):
        d = mt[mt.split_level == level]
        ok = set(d.loc[(d.mcc_ci_low > 0.2) & (d.q_vs_null < 0.05) & (d.delta_vs_null > 0), "trait"])
        out[level] = sorted(t for t in ok if lift_ok.get(t, True))
    return out


def main():
    lab, traits, Y = load_labelled()
    X = load_features("energy_panel", lab.genome_id)
    X = X.loc[:, X.std(axis=0) > 0]
    cols = list(X.columns)
    Xa = X.to_numpy(dtype=np.float32)
    groups = lab.group_genus.to_numpy()
    bundle = {}
    for j, t in enumerate(traits):
        y, m = Y[:, j], ~np.isnan(Y[:, j])
        folds = inner_folds(y[m], groups[m])
        best = None
        for ci, cfg in enumerate(tx.CONFIGS["rf"]):
            oof = np.zeros(m.sum())
            for itr, iva in folds:
                oof[iva] = tx.fit_predict(cfg, Xa[m][itr], y[m][itr], Xa[m][iva])
            ap = safe_ap(y[m], oof)
            if best is None or ap > best[0]:
                best = (ap, cfg, oof)
        ap, cfg, oof = best
        model = FittedModel(cols, cfg, tx.fit_all(cfg, Xa[m], y[m]))
        bundle[t] = {"type": "estimator", "model": model, "threshold": best_mcc_threshold(y[m], oof), "config": cfg,
                     "inner_cv_auprc_all_labelled": ap, "calibrator": fit_calibrator(y[m], oof)}
        print(t, cfg, f"thr={bundle[t]['threshold']:.3f}", f"inner AUPRC={ap:.3f}", flush=True)
    out = RES / "models/RF"
    out.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, out / "model.joblib", compress=3)
    rel = reliable(traits)
    meta = {"arm": "RF", "name": "RedoxKO", "trained_on": "358 labelled genomes (all)", "date": str(date.today()),
            "traits": traits, "feature_set": "energy_panel", "learner": "random forest, 500 trees x 3 seeds, "
            "balanced subsample, min_samples_leaf by inner CV", "reliable_traits": rel["family"],
            "reliable_traits_by_level": rel, "kofam_release": "2026-08-25 (profiles.tar.gz Last-Modified)",
            "labels_sha256": sha256(PROC / "labels.tsv"),
            "cv_performance_file": "results/trees_xai/da/metrics_per_trait.csv (arm RF; outer CV)"}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2))
    r = json.loads((RES / "models/RELEASED.json").read_text())
    r["previous_released_arm"] = r.get("released_arm")
    r["released_arm"] = "RF"
    r["released_name"] = "RedoxKO"
    r["switch_note"] = ("21 October 2026: the random forest on the 852-KO panel replaced B2 (now RedoxKO-LX) after it did "
                        "better on the external genomes (+0.043 macro AUPRC [0.018, 0.069], also without FAPROTAX-derived "
                        "labels), in yes/no calls on new species (MCC +0.056 [0.006, 0.115]) and on incomplete MAGs; "
                        "decision made after these tests, reported in Supplementary Note 1.")
    (RES / "models/RELEASED.json").write_text(json.dumps(r, indent=2))
    print("reliable (calls released) at family level:", rel["family"])


if __name__ == "__main__":
    main()
