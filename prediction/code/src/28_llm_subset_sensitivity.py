"""Seed sensitivity of the 96-genome LLM-baseline subset (integrity audit, M4).

The subset is drawn by a seeded, label-based rule (>=5 positives per trait where possible). This checks how much the
realised draw matters: the same rule is re-run over many seeds and the models' macro AUPRC on each subset is measured
from the existing out-of-fold predictions (no refitting - the models never see the subset).
Output: results/qc/llm_subset_seed_sensitivity.json
"""
import importlib
import json

import numpy as np
import pandas as pd

from redox_common import RES, load_labelled, oof_matrix, safe_ap

b2mod = importlib.import_module("21_baseline2_llm")
N_SEEDS = 200
ARMS = ["B2", "B2XL", "B2E", "A", "null_taxNN"]


def main():
    lab, traits, Y = load_labelled()
    P = {a: oof_matrix(RES / "oof", a, "genus", lab, traits) for a in ARMS}

    def macro(idx, a):
        return float(np.nanmean([safe_ap(Y[idx, j], P[a][idx, j]) for j in range(len(traits))]))

    shipped = b2mod.pick_subset(lab, traits, Y)
    rows = []
    orig = b2mod.SEED
    for s in range(N_SEEDS):
        b2mod.SEED = 100000 + s
        try:
            idx = b2mod.pick_subset(lab, traits, Y)
        except Exception:
            continue
        rows.append({a: macro(np.array(idx), a) for a in ARMS})
    b2mod.SEED = orig
    d = pd.DataFrame(rows)
    full = {a: macro(np.arange(len(lab)), a) for a in ARMS}
    ship = {a: macro(np.array(shipped), a) for a in ARMS}
    res = {"n_seeds": len(d), "full_set_358": full, "shipped_subset_96": ship,
           "resampled_subsets": {a: {"mean": float(d[a].mean()), "sd": float(d[a].std()),
                                     "p05": float(d[a].quantile(0.05)), "p95": float(d[a].quantile(0.95)),
                                     "percentile_of_shipped": float((d[a] < ship[a]).mean() * 100)} for a in ARMS},
           "note": "LLM-vs-model deltas are genome-paired and valid for the shipped subset; this quantifies how "
                   "representative that subset is of the full labelled set."}
    (RES / "qc").mkdir(parents=True, exist_ok=True)
    (RES / "qc/llm_subset_seed_sensitivity.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
