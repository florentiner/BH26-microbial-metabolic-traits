"""Upper bound on "combine everything": exhaustive search over averages of existing out-of-fold predictions.

Every subset of a candidate pool of arms is scored by averaging their out-of-fold probabilities (rank-normalised
per trait so that arms on different probability scales mix fairly). Because the best subset is chosen on the same
out-of-fold data it is scored on, the winning number is OPTIMISTIC - it is an upper bound on what ensembling the
arms we already have could give, not an estimate of what a new ensemble would score. To expose the optimism, the
subset is selected on genus-level folds and its family-level score is reported without re-selection, and vice
versa.

usage: python src/80_ensemble_search.py --task da|pw
Output: results/trees_xai/<task>/ensemble_search.csv
"""
import argparse
import importlib
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

tx = importlib.import_module("47_trees_xai")
from redox_common import safe_ap  # noqa: E402

RES = tx.RES
POOL = {"da": ["B2XL", "B2", "RF", "ET", "NN", "E_large", "C_large", "LG", "TPFN", "ST", "D_nopf", "B1"],
        "pw": ["PB3", "PB2", "ET", "RF", "PB4", "PE_BFL", "PE_MLP", "PB1", "L1"]}
SRC = {"da": [RES / "oof", RES / "trees_xai/da/oof"], "pw": [RES / "pathways/oof", RES / "trees_xai/pw/oof"]}


def load(task, arm):
    for d in SRC[task]:
        f = d / f"{arm}.parquet"
        if f.exists():
            x = pd.read_parquet(f)
            return x[x.arm == arm] if "arm" in x.columns else x
    raise FileNotFoundError(arm)


def matrices(task):
    lab, traits, Y, _ = tx.load_task(task)
    gid = {g: i for i, g in enumerate(lab.genome_id)}
    P = {}
    for arm in POOL[task]:
        try:
            d = load(task, arm)
        except FileNotFoundError:
            continue
        for level in ("genus", "family"):
            s = d[d.level == level]
            M = np.full(Y.shape, np.nan)
            for j, t in enumerate(traits):
                st = s[s.trait == t]
                M[st.genome_id.map(gid).to_numpy(), j] = st.prob.to_numpy()
            if np.isnan(M).any():
                continue
            # rank-normalise within trait so probability scales are comparable across arms
            R = np.column_stack([rankdata(M[:, j]) / len(M) for j in range(M.shape[1])])
            P.setdefault(level, {})[arm] = R
    return lab, traits, Y, P


def macro(Y, M, groups=None):
    return float(np.mean([safe_ap(Y[:, j], M[:, j]) for j in range(Y.shape[1]) if 0 < Y[:, j].sum() < len(Y)]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=["da", "pw"], required=True)
    ap.add_argument("--max-size", type=int, default=5)
    args = ap.parse_args()
    lab, traits, Y, P = matrices(args.task)
    arms = sorted(set(P["genus"]) & set(P["family"]))
    print(f"{args.task}: {len(arms)} arms with complete out-of-fold predictions at both levels: {arms}")
    rows = []
    for k in range(1, args.max_size + 1):
        for combo in itertools.combinations(arms, k):
            g = macro(Y, np.mean([P["genus"][a] for a in combo], axis=0))
            f = macro(Y, np.mean([P["family"][a] for a in combo], axis=0))
            rows.append(dict(size=k, arms="+".join(combo), genus=g, family=f))
    d = pd.DataFrame(rows).sort_values("genus", ascending=False)
    out = RES / f"trees_xai/{args.task}/ensemble_search.csv"
    d.to_csv(out, index=False, float_format="%.4f")
    single = d[d["size"] == 1].set_index("arms")
    print("\nsingle arms (rank-normalised probabilities):")
    print(single[["genus", "family"]].sort_values("genus", ascending=False).to_string())
    print("\nbest 8 combinations selected on GENUS, with their unselected FAMILY score:")
    print(d.head(8).to_string(index=False))
    df = d.sort_values("family", ascending=False)
    print("\nbest 5 selected on FAMILY, with their unselected GENUS score:")
    print(df.head(5).to_string(index=False))
    best_single_g = single.genus.max()
    print(f"\nupper bound from ensembling: genus {d.genus.max():.4f} vs best single {best_single_g:.4f} "
          f"(+{d.genus.max() - best_single_g:.4f}); family {df.family.max():.4f} vs best single "
          f"{single.family.max():.4f} (+{df.family.max() - single.family.max():.4f})")
    print(f"saved {out} ({len(d)} combinations)")


if __name__ == "__main__":
    main()
