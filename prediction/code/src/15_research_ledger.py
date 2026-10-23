"""Research-loop ledger (docs/research_loop.md): evaluate variants against the reference with the promotion rule.

usage: python src/15_research_ledger.py --reference B2 --variants B2cnt B2tax [--s3]
  --s3 : also compute the mean delta (genus level) over the 5 re-drawn fold allocations (requires results/sensitivity/seed<k>__main/oof/<arm>.parquet for both arms)
Appends/updates rows in results/research_ledger.csv (one row per variant x reference).
"""
import argparse
from datetime import datetime

import numpy as np
import pandas as pd

from redox_common import RES, load_labelled, oof_matrix, paired_macro_delta, safe_ap

HYPOTHESIS = {
    "B2cnt": "R1 KO copy number (log1p counts) instead of presence",
    "B2tax": "R2 energy panel + leakage-free taxonomic prior",
    "LG": "R3 MagicLamp LithoGenie/FeGenie HMM counts alone",
    "B2LG": "R3 energy panel + LithoGenie/FeGenie HMM counts",
    "SEL": "R7 per-trait feature-set selection {B2, LG, B2LG} x models inside inner CV",
    "ENS": "R8 equal-weight average of the same model trained on B2 and on LG features",
    "STACK": "R5 cross-trait stacking: level-2 logistic on 14 level-1 (B2+LG logistic) scores",
    "ENSH": "R11 ENS with multiheme cytochrome (CXXCH motif) features added to the LG group",
    "B3": "R12 full significant-KO profile from TogoMCP BH26 RDF (all KOs, not only the energy panel)",
    "B3ENS": "R12b average of B3 model and LithoGenie model",
    "B2B1r": "R15 control: B2 + marker-subset model (no embeddings)",
    "B2RP": "R15 control: B2 + logistic probe on a random projection of the KO panel (embedding-shaped, no PLM)",
    "B2XL": "R15 control: B2 with logreg+xgb averaged (no second feature block)",
    "B2XLE": "R16 B2XL (logreg+xgb on KO panel) averaged with Bacformer-base probe: PLM beyond learner ensembling",
    "B2XLC": "R16 B2XL averaged with pooled ESM-2 probe",
    "B2E": "R13 average of B2 model and Bacformer-base probe",
    "B2C": "R13 average of B2 model and pooled ESM-2 probe",
    "B2LGE": "R14 average of B2 model, LithoGenie model and Bacformer-base probe",
    "D_nopf": "Q3 clean contrast: ABMIL without the KO bag prefilter (pure attention over all proteins)",
    "C_large": "R6 pooled ESM-C 300M embeddings (protocol arm C with the large encoder)",
    "E_large": "R6 Bacformer-large frozen (protocol arm E, preferred checkpoint)",
    "E_large_MAG": "R6 Bacformer-large MAG checkpoint",
    "B2E_large": "R6 average of B2 model and Bacformer-large probe",
    "B2XLE_large": "R6 B2XL averaged with Bacformer-large probe",
    "C": "R4 pooled ESM-2 35M embeddings (pre-registered arm C)",
    "D": "R4 ABMIL over ESM-2 protein bags (pre-registered arm D)",
    "E": "R4 Bacformer-base frozen (pre-registered arm E)",
    "E_MAG": "R4 Bacformer-base MAG checkpoint",
}


def macro(Y, P):
    obs = ~np.isnan(Y) & ~np.isnan(P)
    return np.nanmean([safe_ap(Y[obs[:, j], j], P[obs[:, j], j]) for j in range(Y.shape[1])])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True)
    ap.add_argument("--variants", nargs="+", required=True)
    ap.add_argument("--s3", action="store_true")
    args = ap.parse_args()
    lab, traits, Y = load_labelled()
    rows = []
    for v in args.variants:
        row = {"date": datetime.now().isoformat(timespec="minutes"), "variant": v, "hypothesis": HYPOTHESIS.get(v, ""),
               "reference": args.reference}
        for level in ("genus", "family"):
            Pv = oof_matrix(RES / "oof", v, level, lab, traits)
            Pr = oof_matrix(RES / "oof", args.reference, level, lab, traits)
            common = ~np.isnan(Pv) & ~np.isnan(Pr)  # e.g. B3 lacks the 70 genomes absent from TogoMCP
            Yc = np.where(common, Y, np.nan)
            d, lo, hi = paired_macro_delta(lab, traits, Yc, np.nan_to_num(Pv), np.nan_to_num(Pr), level)
            row[f"{level}_n_genomes_compared"] = int(common.any(axis=1).sum())
            row.update({f"{level}_macro_auprc": macro(Yc, Pv), f"{level}_reference_macro_auprc": macro(Yc, Pr),
                        f"{level}_delta": d, f"{level}_delta_ci_low": lo, f"{level}_delta_ci_high": hi})
        passes = row["genus_delta_ci_low"] > 0 and row["family_delta_ci_high"] >= 0
        row["passes_genus_family_rule"] = passes
        row["s3_mean_delta_genus"] = np.nan
        if args.s3:  # always informative; promotion still requires `passes`
            ds = []
            for k in range(1, 6):
                labk, _, Yk = load_labelled(f"seed{k}")
                d_ = RES / f"sensitivity/seed{k}__main/oof"
                ds.append(macro(Yk, oof_matrix(d_, v, "genus", labk, traits)) -
                          macro(Yk, oof_matrix(d_, args.reference, "genus", labk, traits)))
            row["s3_mean_delta_genus"] = float(np.mean(ds))
        row["promoted"] = bool(passes and args.s3 and row["s3_mean_delta_genus"] > 0)
        rows.append(row)
        print({k: (round(x, 4) if isinstance(x, float) else x) for k, x in row.items()})
    f = RES / "research_ledger.csv"
    old = pd.read_csv(f) if f.exists() else pd.DataFrame()  # append-only: every evaluation ever run is kept
    pd.concat([old, pd.DataFrame(rows)], ignore_index=True).to_csv(f, index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
