"""Simultaneous inference over every research-loop variant (audit round 2, M1).

For all arms with full-coverage out-of-fold predictions in results/oof (pre-registered + research loop), the paired
cluster-bootstrap distribution of macro-AUPRC(variant) - macro-AUPRC(reference) is computed on shared resamples.
Max-T: each variant's bootstrap deltas are standardised by their bootstrap SD; the 95th percentile of the maximum
|t| over variants gives a family-wise critical value, so the simultaneous 95% CIs (delta +/- c * sd) hold jointly
across all variants tried. Arms with partial genome coverage (B3, B3ENS) are listed but excluded (different genomes).
Output: results/simultaneous_inference.csv, results/variants_tried.csv
"""
import numpy as np
import pandas as pd

from redox_common import LEVELS, RES, load_labelled, oof_matrix, safe_ap

N_BOOT = 10000  # audit round 4 (M-7): 1000 left the max-T critical value within Monte-Carlo noise
SEED = 20260918
NULLS = {"null_prevalence", "null_taxNN"}


def main():
    lab, traits, Y = load_labelled()
    arms = sorted({f.stem for f in (RES / "oof").glob("*.parquet")} - {"nulls"})
    rows, tried = [], []
    for level in LEVELS:
        P = {a: oof_matrix(RES / "oof", a, level, lab, traits) for a in arms}
        full = [a for a in arms if not np.isnan(P[a]).any()]
        partial = [a for a in arms if a not in full]
        if level == "genus":
            tried = [dict(arm=a, full_coverage=a in full) for a in arms]
        groups = lab[f"group_{level}"].to_numpy()
        ug = np.unique(groups)
        members = {g: np.where(groups == g)[0] for g in ug}
        rng = np.random.default_rng(SEED)
        boots = [np.concatenate([members[g] for g in rng.choice(ug, len(ug))]) for _ in range(N_BOOT)]
        macro = {a: np.array([np.nanmean([safe_ap(Y[b, j], P[a][b, j]) for j in range(len(traits))]) for b in boots])
                 for a in full}
        point = {a: np.nanmean([safe_ap(Y[:, j], P[a][:, j]) for j in range(len(traits))]) for a in full}
        for ref in ("B2", "B2XL"):
            if ref not in full:
                continue
            variants = [a for a in full if a != ref and a not in NULLS]
            D = np.stack([macro[a] - macro[ref] for a in variants])  # variants x boots
            sd = D.std(axis=1, ddof=1)
            delta = np.array([point[a] - point[ref] for a in variants])
            T = np.abs((D - delta[:, None]) / sd[:, None])
            c = float(np.percentile(T.max(axis=0), 95))
            for a, d, s, row in zip(variants, delta, sd, D):
                rows.append(dict(split_level=level, reference=ref, variant=a, macro_auprc=point[a], delta=d,
                                 marginal_ci_low=np.percentile(row, 2.5), marginal_ci_high=np.percentile(row, 97.5),
                                 simultaneous_ci_low=d - c * s, simultaneous_ci_high=d + c * s, maxT_critical=c,
                                 n_variants_in_family=len(variants),
                                 simultaneous_significant=bool(d - c * s > 0 or d + c * s < 0),
                                 note="family-wise over the variants of this reference only; the two reference "
                                      "families (B2, B2XL) are not jointly corrected"))
        print(level, "partial-coverage arms excluded:", partial)
    out = pd.DataFrame(rows)
    out.to_csv(RES / "simultaneous_inference.csv", index=False, float_format="%.4f")
    pd.DataFrame(tried).to_csv(RES / "variants_tried.csv", index=False)
    print(out[out.reference == "B2"].sort_values(["split_level", "delta"], ascending=[True, False])
          [["split_level", "variant", "macro_auprc", "delta", "simultaneous_ci_low", "simultaneous_ci_high",
            "simultaneous_significant"]].round(4).to_string())


if __name__ == "__main__":
    main()
