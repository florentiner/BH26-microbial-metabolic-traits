"""Table 5 for the final model RedoxKO (random forest): agreement of MAG and isolate calls and change in AUROC, with the
analysis of src/88_real_mags.py unchanged; predictions from the KO panels of the pairs (src/141_mag_rf_vs_redoxko.py).
Output: results/additional/rf_mags/x2_*.csv
"""
import importlib
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
m88 = importlib.import_module("88_real_mags")


def rf_predictions():
    X = (pd.read_parquet(ROOT / "data/additional/mags/ko_panel_pairs.parquet") > 0).astype(np.float32)
    rows = []
    for t, b in joblib.load(ROOT / "results/models/RF/model.joblib").items():
        s = np.asarray(b["model"].predict_proba(X.reindex(columns=b["model"].columns).fillna(0)))
        rows += [dict(genome=g, trait=t, score=v, predicted=int(v >= b["threshold"])) for g, v in zip(X.index, s)]
    return pd.DataFrame(rows)


m88.OUT = ROOT / "results/additional/rf_mags"
m88.OUT.mkdir(parents=True, exist_ok=True)
m88.load_predictions = rf_predictions
m88.analyze()
