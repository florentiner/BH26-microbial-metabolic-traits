"""Collect the Kaggle fine-tuning results (Addendum N) into the project's formats.

Reads data/kaggle/finetune/job_<ARM>_<level>_<fold>.json (kaggle/bacformer_finetune) and writes
  results/trees_xai/da/oof/<ARM>.parquet   out-of-fold predictions, same columns as every other arm; the per-trait
                                          cut-off is the MCC-maximising threshold on that fold's validation genomes
  results/trees_xai/da/external_<ARM>.csv  scores of the final model on the 562 lineage-novel external genomes

usage: python src/61_finetune_collect.py [pw]      (pw = the pathway fine-tuning of the pathway FM addendum)
"""
import importlib
import json
import sys

import numpy as np
import pandas as pd

from redox_common import best_mcc_threshold

tx = importlib.import_module("47_trees_xai")
TASK = "pw" if "pw" in sys.argv[1:] else "da"
IN = tx.RES.parent / ("data/kaggle/pw_finetune" if TASK == "pw" else "data/kaggle/finetune")
OUT = tx.OUT / TASK
ARMS = ("FTP_BFL2", "FTP_BFL0") if TASK == "pw" else ("FT_BFL2", "FT_BFL0")


def main():
    lab, traits, Y, _ = tx.load_task(TASK)
    ids = list(lab.genome_id)
    for arm in ARMS:
        rows = []
        for level in ("genus", "family"):
            for k in range(5):
                d = json.load(open(IN / f"job_{arm}_{level}_{k}.json"))
                vi, vp = np.array(d["val_idx"]), np.array(d["val_pred"])
                ti, tp = np.array(d["test_idx"]), np.array(d["test_pred"])
                for j, t in enumerate(traits):
                    thr = best_mcc_threshold(Y[vi, j], vp[:, j]) if 0 < Y[vi, j].sum() < len(vi) else 0.5
                    for i, p in zip(ti, tp[:, j]):
                        rows.append((arm, level, t, ids[i], k, Y[i, j], float(p), float(thr),
                                     json.dumps({"model": arm, "best_epoch": d["best_epoch"]}), d["val_auprc"]))
        o = pd.DataFrame(rows, columns=["arm", "level", "trait", "genome_id", "fold", "y", "prob", "threshold",
                                        "config", "inner_auprc"])
        assert o.groupby("level").genome_id.nunique().eq(len(ids)).all(), "every genome must be predicted once"
        o.to_parquet(OUT / "oof" / f"{arm}.parquet", index=False)
        if not (IN / f"job_{arm}_final_-1.json").exists():   # pathway run: cross-validation only, no external set
            print(f"{arm}: {len(o)} out-of-fold rows (no external job)")
            continue
        f = json.load(open(IN / f"job_{arm}_final_-1.json"))
        S = np.array(f["ext_scores"])
        ext = pd.DataFrame([dict(genome_id=g, trait=t, arm=arm, score=float(S[i, j]))
                            for i, g in enumerate(f["ext_ids"]) for j, t in enumerate(traits)])
        ext.to_csv(OUT / f"external_{arm}.csv", index=False, float_format="%.5f")
        print(f"{arm}: {len(o)} out-of-fold rows, {len(f['ext_ids'])} external genomes")


if __name__ == "__main__":
    main()
