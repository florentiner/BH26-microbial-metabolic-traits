"""Foundation-model arms for the pathway task (docs/preregistration_pathways.md, FM addendum).

Vectors: ESM-C 300M per protein -> Bacformer-large (masked-complete-genomes @ 84f85c37), mean over proteins,
computed on Kaggle GPUs (data/kaggle/pw_embed_1..4) plus the genomes already embedded for the donor/acceptor task.

Arms (same frozen folds, nested CV and thresholds as PB2):
  PE_BFL : L2 logistic probe on the 960-number genome vector (the foundation model without gene content)
  PE_MLP : one multi-label MLP (20 outputs) on the same vector - the PB4 network on embeddings instead of KOs
  PE     : 1/2 PB2 + 1/2 PE_BFL, averaged out of fold (no extra training)
  FTP_FUSE : 1/2 PB2 + 1/2 FTP_BFL2 (fine-tuning addendum), also a pure out-of-fold average

usage: python src/62_pathway_fm.py {PE_BFL,PE_MLP,PE,FTP_FUSE} [--jobs N]
Output: results/trees_xai/pw/oof/<arm>.parquet (where 47_trees_xai.py metrics --task pw reads the arms)
"""
import argparse
import importlib
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

a32 = importlib.import_module("32_pathway_arms")
tx = importlib.import_module("47_trees_xai")
from pathway_common import RES, load_pathways  # noqa: E402
from redox_common import LEVELS, nested_cv  # noqa: E402

ROOT = RES.parent
OUT = tx.OUT / "pw" / "oof"
EMB = "bacformer_large_mean"
FUSIONS = {"PE": ("PB2", "PE_BFL"), "FTP_FUSE": ("PB2", "FTP_BFL2")}


def embeddings(genome_ids):
    """Bacformer-large genome vectors for every pathway genome, from the Kaggle kernels and the earlier runs."""
    srcs = [ROOT / f"data/kaggle/pw_embed_{i}/genome_features_{EMB}.parquet" for i in (1, 2, 3, 4)]
    srcs += [ROOT / f"data/processed/genome_features_{EMB}.parquet", ROOT / f"data/kaggle/large_external/genome_features_{EMB}.parquet"]
    parts = [pd.read_parquet(p) for p in srcs if Path(p).exists()]
    X = pd.concat(parts)
    X = X[~X.index.duplicated()]
    missing = set(genome_ids) - set(X.index)
    assert not missing, f"{len(missing)} genomes without a vector, e.g. {sorted(missing)[:5]}"
    return X.loc[list(genome_ids)].astype(np.float32)


def oof_path(arm):
    for p in (OUT / f"{arm}.parquet", RES / f"pathways/oof/{arm}.parquet"):
        if p.exists():
            return p
    raise SystemExit(f"{arm}: out-of-fold predictions not found")


def fuse(arm):
    """Average two arms' out-of-fold probabilities; thresholds re-picked on the same rows are not available out of
    fold, so the fused threshold is the mean of the two arms' fold thresholds (as for B2XL on donors/acceptors)."""
    a, b = FUSIONS[arm]
    A, B = (pd.read_parquet(oof_path(x)) for x in (a, b))
    key = ["level", "trait", "genome_id", "fold"]
    m = A.merge(B, on=key, suffixes=("_a", "_b"))
    assert len(m) == len(A) == len(B), f"{arm}: {len(A)} vs {len(B)} rows, merged {len(m)}"
    assert (m.y_a == m.y_b).all()
    return pd.DataFrame({"arm": arm, "level": m.level, "trait": m.trait, "genome_id": m.genome_id, "fold": m.fold,
                         "y": m.y_a, "prob": (m.prob_a + m.prob_b) / 2,
                         "threshold": (m.threshold_a + m.threshold_b) / 2,
                         "config": f"mean({a},{b})", "inner_auprc": m.inner_auprc_a})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arm", choices=["PE_BFL", "PE_MLP", "PE", "FTP_FUSE"])
    ap.add_argument("--jobs", type=int, default=6)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.arm in FUSIONS:
        df = fuse(args.arm)
    else:
        lab, traits, Y = load_pathways()
        X = embeddings(lab.genome_id).to_numpy()
        X = X[:, X.std(axis=0) > 0]
        print(f"{len(lab)} genomes, {len(traits)} traits, feature matrix {X.shape}", flush=True)
        if args.arm == "PE_MLP":
            out = [a32.run_mlp(args.arm, X, lab, traits, Y, lvl) for lvl in LEVELS]
        else:
            cfgs = [c for c in a32.configs() if c["model"] == "logreg"]   # probe only, as pre-registered
            out = Parallel(n_jobs=args.jobs)(
                delayed(nested_cv)(
                    args.arm, X, lab, [traits[j]], Y[:, [j]], cfgs, a32.fit_predict, lvl)
                for lvl in LEVELS for j in range(len(traits)))
        df = pd.concat(out, ignore_index=True)
    df.to_parquet(OUT / f"{args.arm}.parquet", index=False)
    print(OUT / f"{args.arm}.parquet", df.groupby(["arm", "level"]).size().to_dict())


if __name__ == "__main__":
    main()
