"""X5 (docs/preregistration_additional_experiments.md): a stratified random sample of 50 of the 255 primary growth
tests (seed 20261014; at least 3 per category where available) exported for two independent human reviewers.
Analysis (agreement with the LLM extraction, Cohen's kappa) is `python src/91_label_check_sheet.py analyze` and runs
only when both reviewer columns are completely filled.

usage: python src/91_label_check_sheet.py [analyze]
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/additional"
SHEET = OUT / "x5_label_check_sheet.csv"
SEED = 20261014
N = 50


def make():
    calls = pd.read_csv(ROOT / "validation_new_organisms/scored178_B2_balanced100.csv")
    calls = calls[calls.blind & ~calls.category.isin(["donor_other", "acceptor_other"])].reset_index(drop=True)
    ev = json.loads((ROOT / "validation_new_organisms/validation_set.json").read_text())
    rng = np.random.default_rng(SEED)
    pick = []
    for c, d in calls.groupby("category"):
        pick += list(rng.choice(d.index, min(3, len(d)), replace=False))
    rest = calls.index.difference(pick)
    pick += list(rng.choice(rest, N - len(pick), replace=False))
    S = calls.loc[sorted(pick)].copy()

    by_acc = {o.get("refseq_accession"): o for o in ev}

    def lookup(acc, cat):
        o = by_acc.get(acc, {})
        c = o.get("categories", {}).get(cat, {})
        return c.get("quote", ""), o.get("description_url", "")
    q = [lookup(a, c) for a, c in zip(S.accession, S.category)]
    S["quote"], S["source"] = [x[0] for x in q], [x[1] for x in q]
    S["llm_label"] = S.published
    S = S[["accession", "organism", "paper", "category", "quote", "source", "llm_label"]]
    S["reviewer1_label (tested_positive/tested_negative/not_tested)"] = ""
    S["reviewer1_note"] = ""
    S["reviewer2_label (tested_positive/tested_negative/not_tested)"] = ""
    S["reviewer2_note"] = ""
    S.to_csv(SHEET, index=False)
    print(len(S), "calls written to", SHEET, "| quotes found:", int((S.quote != "").sum()))


def analyze():
    from sklearn.metrics import cohen_kappa_score
    S = pd.read_csv(SHEET).fillna("")
    r1, r2 = S.columns[7], S.columns[9]
    if (S[r1] == "").any() or (S[r2] == "").any():
        raise SystemExit("both reviewer columns must be filled before analysis")
    print("reviewer1 vs reviewer2 kappa:", round(cohen_kappa_score(S[r1], S[r2]), 3))
    agree = S[S[r1] == S[r2]]
    print("calls where reviewers agree:", len(agree), "| LLM extraction matches them:",
          round((agree.llm_label == agree[r1]).mean(), 3))


if __name__ == "__main__":
    analyze() if len(sys.argv) > 1 and sys.argv[1] == "analyze" else make()
