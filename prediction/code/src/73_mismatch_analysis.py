"""Classify every mismatch in the prospective validation against the organism's own marker genes.

For each wrong call the genome's canonical marker groups (redox_common.RULE_GROUPS, the same mapping the rule
baseline uses) are read from the KO panel computed for the validation genomes. Classes:

  feature_gap     the category has no marker KOs in the panel at all (iron reduction, most of "other"), so no
                  KO-based model can see the chemistry - a missed positive here is a limit of the feature space
  model_error     the genome's marker status agrees with the published result and the model still called the
                  opposite: a genuine error
  genotype_gap    a false positive where the genome does carry a complete marker group - the capability is
                  encoded but the organism did not use it under the tested condition
  marker_absent   a missed positive where no marker group is complete: the published capability is not visible in
                  this genome's gene content (divergent enzymes, or annotation recall)
  unreliable_cat  the category is not in this arm's released reliable-trait set, so the model emits no call there

usage: python src/73_mismatch_analysis.py
Output: validation_new_organisms/mismatch178_<arm>.csv and a printed summary
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from redox_common import RULE_GROUPS  # noqa: E402

V = Path(__file__).resolve().parents[1] / "validation_new_organisms"


def marker_state(X, acc, cat):
    """(n complete groups, n groups, joined names of the complete ones) for this category in this genome."""
    groups = RULE_GROUPS.get(cat, {})
    if not groups:
        return 0, 0, ""
    done = []
    for name, kos in groups.items():
        present = [k for k in kos if k in X.columns and X.loc[acc, k] > 0]
        if present and len(present) == len([k for k in kos if k in X.columns]):
            done.append(name)
    return len(done), len(groups), ",".join(done)


def classify(row, X):
    if not row.deployed:
        return "unreliable_cat"
    n_done, n_groups, _ = marker_state(X, row.accession, row.category)
    if n_groups == 0:
        return "feature_gap"
    if row.expected == 1 and row.call == 0:
        return "marker_absent" if n_done == 0 else "model_error"
    if row.expected == 0 and row.call == 1:
        return "genotype_gap" if n_done > 0 else "model_error"
    return "correct"


def main():
    X = pd.read_parquet(V / "ko_panel_validation.parquet")
    for arm in ("B2", "B2XLE_large"):
        d = pd.read_csv(V / f"scored178_{arm}.csv")
        d = d[d.accession.isin(X.index)]
        d["marker_groups_complete"] = [marker_state(X, r.accession, r.category)[2] for r in d.itertuples()]
        d["verdict"] = [classify(r, X) if r.agree == 0 else "correct" for r in d.itertuples()]
        d.to_csv(V / f"mismatch178_{arm}.csv", index=False)
        m = d[d.agree == 0]
        print(f"\n== {arm}: {len(m)} mismatches of {len(d)} calls")
        print(m.verdict.value_counts().to_string())
        prim = m[m.blind & ~m.category.isin({"donor_other", "acceptor_other"})]
        print(f"   primary set only ({len(prim)} mismatches):")
        print(prim.verdict.value_counts().to_string())
        print("   by category (primary):")
        print(prim.groupby(["category", "verdict"]).size().to_string())


if __name__ == "__main__":
    main()
