"""Put the two tasks on the same footing (pre-registered prediction PP5).

Raw AUPRC cannot be compared across tasks with different prevalences: a trait present in 1 % of genomes has a
prevalence floor of 0.01, one present in 50 % has a floor of 0.5. The comparable quantities are the lift over
prevalence, the gain over the taxonomic group-mean null, and MCC (threshold-based, prevalence-aware).

Output: results/task_comparison.csv + a printed summary.
usage: python src/39_task_comparison.py
"""
import pandas as pd

from pathway_common import RES

RELEASED_DA = "B2"  # released donor/acceptor arm


def summarise(df, arm, task, level="genus"):
    d = df[(df.arm == arm) & (df.split_level == level)]
    return dict(task=task, arm=arm, split_level=level, n_traits=len(d),
                mean_prevalence=d.prevalence.mean(),
                macro_auprc=d.auprc.mean(),
                macro_auprc_lift=(d.auprc / d.prevalence).mean(),
                macro_mcc=d.mcc.mean(),
                macro_delta_vs_taxonomic_null=d.delta_vs_null.mean(),
                n_beating_null_q05=int(((d.delta_vs_null_q_bh < 0.05) & (d.delta_vs_null > 0)).sum()),
                frac_beating_null=float(((d.delta_vs_null_q_bh < 0.05) & (d.delta_vs_null > 0)).mean()))


def main():
    da = pd.read_csv(RES / "metrics_per_trait.csv")
    pw = pd.read_csv(RES / "pathways/metrics_per_trait.csv")
    pw_arm = (pw[(pw.split_level == "genus") & pw.arm.str.startswith("PB")]
              .groupby("arm").auprc.mean().idxmax())
    rows = []
    for level in ("genus", "family"):
        rows.append(summarise(da, RELEASED_DA, "donors/acceptors", level))
        rows.append(summarise(pw, pw_arm, "pathways", level))
        rows.append(summarise(pw, "PA", "pathways (KEGG module rule)", level))
    out = pd.DataFrame(rows)
    out.to_csv(RES / "task_comparison.csv", index=False, float_format="%.4f")
    print(out.round(3).to_string(index=False))

    g = out[out.split_level == "genus"].set_index("task")
    pp5 = ("supported" if g.at["pathways", "macro_auprc"] > g.at["donors/acceptors", "macro_auprc"]
           else "not supported")
    pp5_lift = ("supported" if g.at["pathways", "macro_auprc_lift"] > g.at["donors/acceptors", "macro_auprc_lift"]
                else "not supported")
    print("\nPP5 (pathways easier than donor/acceptor categories):")
    print(f"  on raw macro AUPRC : {pp5} "
          f"({g.at['pathways', 'macro_auprc']:.3f} vs {g.at['donors/acceptors', 'macro_auprc']:.3f}; "
          f"mean prevalence {g.at['pathways', 'mean_prevalence']:.3f} vs "
          f"{g.at['donors/acceptors', 'mean_prevalence']:.3f} - not a like-for-like comparison)")
    print(f"  on lift over prevalence : {pp5_lift} "
          f"({g.at['pathways', 'macro_auprc_lift']:.1f}x vs {g.at['donors/acceptors', 'macro_auprc_lift']:.1f}x)")
    print(f"  on macro MCC : {g.at['pathways', 'macro_mcc']:.3f} vs "
          f"{g.at['donors/acceptors', 'macro_mcc']:.3f}")
    print(f"  traits beating the taxonomic null: "
          f"{g.at['pathways', 'n_beating_null_q05']:.0f}/{g.at['pathways', 'n_traits']:.0f} vs "
          f"{g.at['donors/acceptors', 'n_beating_null_q05']:.0f}/{g.at['donors/acceptors', 'n_traits']:.0f}")
    print("\nThe primary, prevalence-free statement is the MCC and the gain over the taxonomic null; raw macro "
          "AUPRC is reported for both tasks but is not comparable between them.")


if __name__ == "__main__":
    main()
