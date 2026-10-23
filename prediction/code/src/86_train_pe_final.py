"""Save the task-2 deep-learning fusion PE as a deployable model: PB2's saved KO model averaged with an L2 logistic
probe on the Bacformer-large genome vector, exactly as the PE arm was defined in cross-validation
(docs/preregistration_pathways.md, FM addendum). Trained on all 3 378 labelled pathway genomes.

usage: python src/86_train_pe_final.py
Output: results/models/pathways/PE/{model.joblib, metadata.json}
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
from pathway_common import load_pathways  # noqa: E402
from redox_common import FittedFusion, FittedModel, best_mcc_threshold, fit_calibrator, fit_models, inner_folds, safe_ap  # noqa: E402

fm = importlib.import_module("62_pathway_fm")
LR_C = [0.01, 0.1, 1.0, 10.0]


def main():
    lab, traits, Y = load_pathways()
    groups = lab.group_genus.to_numpy()
    emb = fm.embeddings(lab.genome_id)
    emb.index = lab.genome_id.to_numpy()
    E = emb.to_numpy(dtype=np.float32)
    pb2 = joblib.load(ROOT / "results/models/pathways/PB2/model.joblib")
    out = {}
    for j, t in enumerate(traits):
        y = Y[:, j]
        folds = inner_folds(y, groups)
        scores, oofs = {}, {}
        for ci, C in enumerate(LR_C):
            oof = np.zeros(len(y))
            for itr, iva in folds:
                m = fit_models({"model": "logreg", "C": C}, E[itr], y[itr])
                oof[iva] = np.mean([mm.predict_proba(E[iva])[:, 1] for mm in m], axis=0)
            scores[ci], oofs[ci] = safe_ap(y, oof), oof
        bi = max(scores, key=lambda c: (np.nan_to_num(scores[c], nan=-1), -c))
        probe = FittedModel(list(emb.columns), {"model": "logreg", "C": LR_C[bi]},
                            fit_models({"model": "logreg", "C": LR_C[bi]}, E, y))
        fused = FittedFusion([pb2[t]["model"], probe], {"fusion": "mean(PB2, bacformer_large_probe)"})
        # cut-off from the averaged inner out-of-fold scores of the two components
        oof_pb2 = pb2[t].get("oof_inner")
        oof_fused = 0.5 * (oofs[bi] + oof_pb2) if oof_pb2 is not None and len(oof_pb2) == len(y) else oofs[bi]
        out[t] = {"model": fused, "threshold": best_mcc_threshold(y, oof_fused), "config": fused.cfg,
                  "calibrator": fit_calibrator(y, oof_fused)}
        print(f"  {t}: probe C={LR_C[bi]}, inner AUPRC {scores[bi]:.3f}", flush=True)
    d = ROOT / "results/models/pathways/PE"
    d.mkdir(parents=True, exist_ok=True)
    joblib.dump(out, d / "model.joblib")
    json.dump({"arm": "PE", "trained_on": "3378 labelled genomes (all)", "date": str(date.today()),
               "feature_set": "bh26_full_ko_profile (PB2 model) + bacformer_large_mean (logistic probe), averaged",
               "note": "deployment requires a Bacformer-large genome vector; PB3 is the released arm"},
              open(d / "metadata.json", "w"), indent=1)
    print("saved", d)


if __name__ == "__main__":
    main()
