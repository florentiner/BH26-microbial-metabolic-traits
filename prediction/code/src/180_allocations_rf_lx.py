"""Review round 3: macro AUPRC of RedoxKO (RF) and RedoxKO-LX (B2) in every fold allocation, the paired differences
RF minus B2 and RedoxKO-ens-BF (B2XLE_large) minus RedoxKO-ens (B2XL), and the six-allocation means (paper folds and
five further random allocations). Paired genus- or family-cluster bootstrap (1,000 resamples, seed 20260917).
Output: results/allocations_rf_lx.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import RES, load_labelled, oof_matrix, paired_macro_delta, safe_ap  # noqa: E402

ALLOC = {"paper": (None, RES / "oof", RES / "trees_xai/da/oof")}
for s in ["seed1", "seed2", "seed3", "seed4", "seed5", "gtdb"]:
    d = RES / f"sensitivity/{s}__main/oof"
    ALLOC[s] = (s, d, d)


def macro(P, Y, traits):
    return float(np.mean([safe_ap(Y[:, j], P[:, j]) for j in range(len(traits))]))


rows = []
for name, (split, d_b, d_rf) in ALLOC.items():
    lab, traits, Y = load_labelled(split=split) if split else load_labelled()
    for level in ("genus", "family"):
        P = {"RF": oof_matrix(d_rf, "RF", level, lab, traits), "B2": oof_matrix(d_b, "B2", level, lab, traits),
             "B2XL": oof_matrix(d_b, "B2XL", level, lab, traits), "B2XLE_large": oof_matrix(d_b, "B2XLE_large", level, lab, traits)}
        r = dict(allocation=name, level=level, **{f"macro_{k}": macro(v, Y, traits) for k, v in P.items()})
        r["RF_minus_B2"], r["RF_minus_B2_lo"], r["RF_minus_B2_hi"] = paired_macro_delta(lab, traits, Y, P["RF"], P["B2"], level)
        r["ensBF_minus_ens"], r["ensBF_minus_ens_lo"], r["ensBF_minus_ens_hi"] = paired_macro_delta(
            lab, traits, Y, P["B2XLE_large"], P["B2XL"], level)
        rows.append(r)
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
out = pd.DataFrame(rows)
out.to_csv(RES / "allocations_rf_lx.csv", index=False, float_format="%.4f")
six = out[out.allocation != "gtdb"].groupby("level")[[c for c in out.columns if c.startswith("macro_")]].agg(["mean", "std"])
print(six.round(3).to_string())
