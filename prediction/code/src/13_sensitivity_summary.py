"""Collect primary and sensitivity analyses into summary tables.

Output:
  results/sensitivity_summary_macro.csv     macro AUPRC per analysis x arm x split level
  results/sensitivity_summary_contrasts.csv macro paired contrasts per analysis
  results/sensitivity_fold_seeds.csv        S3: across the frozen partition + 5 re-drawn allocations, per arm and
                                            level: mean/min/max macro AUPRC, rank stability, genus-family gap
"""
import pandas as pd

from redox_common import RES

ANALYSES = {
    "primary": RES,
    "S1_block_nan": RES / "sensitivity/frozen__block_nan",
    "S2_gtdb": RES / "sensitivity/gtdb__main",
    "S4_quality_ok": RES / "sensitivity/frozen__main__quality_ok",
    **{f"S3_seed{k}": RES / f"sensitivity/seed{k}__main" for k in range(1, 6)},
}


def main():
    macro, contr = [], []
    for name, d in ANALYSES.items():
        if not (d / "metrics_macro.csv").exists():
            print("missing", name)
            continue
        m = pd.read_csv(d / "metrics_macro.csv")
        m.insert(0, "analysis", name)
        macro.append(m)
        c = pd.read_csv(d / "paired_comparisons.csv")
        c = c[c.trait == "MACRO"].drop(columns=["trait", "q_bh"])
        c.insert(0, "analysis", name)
        contr.append(c)
    macro = pd.concat(macro, ignore_index=True)
    contr = pd.concat(contr, ignore_index=True)
    macro.to_csv(RES / "sensitivity_summary_macro.csv", index=False, float_format="%.4f")
    contr.to_csv(RES / "sensitivity_summary_contrasts.csv", index=False, float_format="%.4f")

    s3 = macro[macro.analysis.str.match(r"primary|S3_seed")].copy()
    s3["rank"] = s3.groupby(["analysis", "split_level"]).macro_auprc.rank(ascending=False)
    agg = s3.groupby(["arm", "split_level"]).agg(n_partitions=("analysis", "nunique"),
                                                 macro_auprc_mean=("macro_auprc", "mean"),
                                                 macro_auprc_min=("macro_auprc", "min"),
                                                 macro_auprc_max=("macro_auprc", "max"),
                                                 macro_auprc_sd=("macro_auprc", "std"),
                                                 rank_min=("rank", "min"), rank_max=("rank", "max")).reset_index()
    gap = s3.pivot_table(index=["analysis", "arm"], columns="split_level", values="macro_auprc").reset_index()
    gap["genus_minus_family"] = gap["genus"] - gap["family"]
    g = gap.groupby("arm").genus_minus_family.agg(["mean", "min", "max",
                                                   lambda x: int((x > 0).sum())]).reset_index()
    g.columns = ["arm", "genus_minus_family_mean", "genus_minus_family_min", "genus_minus_family_max",
                 "n_partitions_genus_gt_family"]
    agg = agg.merge(g, on="arm", how="left")
    agg.to_csv(RES / "sensitivity_fold_seeds.csv", index=False, float_format="%.4f")
    print(macro.pivot_table(index=["arm", "split_level"], columns="analysis", values="macro_auprc").round(3)
          .to_string())
    print(agg.to_string())


if __name__ == "__main__":
    main()
