"""Pre-registration BT (docs/preregistration_BT_bacformer_tree_learners.md): score the five learners fitted on the
Bacformer-large genome vector against logistic regression on the same vector (arm E_large).

Primary: macro AUPRC (14 categories) at genus level, paired difference learner - E_large, genus-cluster bootstrap
(redox_common.paired_macro_delta, 1,000 resamples, seed 20260917). Secondary: family level. Intervals are also given
Bonferroni-adjusted for the five comparisons per level (99.5% two-sided).
Output: results/trees_xai/da_bf/bt_summary.csv
"""
import numpy as np
import pandas as pd

from pathway_common import RES
from redox_common import LEVELS, load_labelled, oof_matrix, safe_ap

ARMS = ["RF_BF", "ET_BF", "XGB_BF", "TPFN_BF", "AM_BF"]
BF = RES / "trees_xai/da_bf/oof"
REF = RES / "oof"
N_BOOT, SEED = 1000, 20260917


def macro(P, Y, idx):
    obs = ~np.isnan(Y)
    return np.nanmean([safe_ap(Y[idx, j][obs[idx, j]], P[idx, j][obs[idx, j]]) for j in range(Y.shape[1])])


def main():
    lab, traits, Y = load_labelled()
    rows = []
    for level in LEVELS:
        Pref = oof_matrix(REF, "E_large", level, lab, traits)
        groups = lab[f"group_{level}"].to_numpy()
        ug = np.unique(groups)
        members = {g: np.where(groups == g)[0] for g in ug}
        rng = np.random.default_rng(SEED)
        boots = [np.concatenate([members[g] for g in rng.choice(ug, len(ug))]) for _ in range(N_BOOT)]
        allidx = np.arange(len(lab))
        ref_point = macro(Pref, Y, allidx)
        ref_boot = np.array([macro(Pref, Y, b) for b in boots])
        rows.append(dict(arm="E_large", level=level, macro_auprc=ref_point,
                         ci_low=np.percentile(ref_boot, 2.5), ci_high=np.percentile(ref_boot, 97.5)))
        for arm in ARMS:
            f = BF / f"{arm}.parquet"
            if not f.exists():
                print("missing", f)
                continue
            P = oof_matrix(BF, arm, level, lab, traits)
            point = macro(P, Y, allidx)
            pb = np.array([macro(P, Y, b) for b in boots])
            d = pb - ref_boot
            rows.append(dict(arm=arm, level=level, macro_auprc=point,
                             ci_low=np.percentile(pb, 2.5), ci_high=np.percentile(pb, 97.5),
                             delta_vs_LR=point - ref_point,
                             delta_ci_low=np.percentile(d, 2.5), delta_ci_high=np.percentile(d, 97.5),
                             delta_ci_low_bonf=np.percentile(d, 0.5), delta_ci_high_bonf=np.percentile(d, 99.5)))
    out = pd.DataFrame(rows)
    (RES / "trees_xai/da_bf").mkdir(parents=True, exist_ok=True)
    out.to_csv(RES / "trees_xai/da_bf/bt_summary.csv", index=False, float_format="%.4f")
    print(out.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
