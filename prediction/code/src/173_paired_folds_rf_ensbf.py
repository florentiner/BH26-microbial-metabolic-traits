"""Paired macro-AUPRC differences between RedoxKO-ens-BF (B2XLE_large) and RedoxKO (random forest) in every fold
allocation: the folds of the paper, five further random allocations and GTDB-based folds, at genus and family level.
Paired genus- or family-cluster bootstrap (1,000 resamples, seed 20260917), as for all paired contrasts.
Output: results/paired_folds_rf_ensbf.csv
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import RES, load_labelled, oof_matrix, paired_macro_delta  # noqa: E402

ALLOC = {"paper": (None, RES / "oof", RES / "trees_xai/da/oof")}
for s in ["seed1", "seed2", "seed3", "seed4", "seed5", "gtdb"]:
    d = RES / f"sensitivity/{s}__main/oof"
    ALLOC[s] = (s, d, d)

rows = []
for name, (split, d_bf, d_rf) in ALLOC.items():
    lab, traits, Y = load_labelled(split=split) if split else load_labelled()
    for level in ("genus", "family"):
        Pbf = oof_matrix(d_bf, "B2XLE_large", level, lab, traits)
        Prf = oof_matrix(d_rf, "RF", level, lab, traits)
        d, lo, hi = paired_macro_delta(lab, traits, Y, Pbf, Prf, level)
        rows.append(dict(allocation=name, level=level, ensBF_minus_RedoxKO=d, ci_low=lo, ci_high=hi,
                         significant=bool(lo > 0 or hi < 0)))
        print(rows[-1], flush=True)
out = pd.DataFrame(rows)
out.to_csv(RES / "paired_folds_rf_ensbf.csv", index=False, float_format="%.4f")
print(out.round(3).to_string(index=False))
