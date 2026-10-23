"""MAG-robustness: how pathway prediction degrades on incomplete proteomes (docs/preregistration_pathways.md).

The application is metagenome-assembled genomes, which are partial. Proteins are randomly removed from each
genome to leave 90 / 70 / 50 % of the proteome (3 draws, fixed seeds), KO features are recomputed from the
surviving proteins only, and the models selected on the *complete* training fold are applied to the degraded
test-fold genomes. Training always uses complete genomes: the degradation is a property of the query, not of
the training set, which is what deployment on MAGs looks like.

usage: python src/34_pathway_mag_robustness.py [--arm PB2] [--level genus] [--jobs 6]
Output: results/pathways/mag_robustness.csv (+ _macro.csv)
"""
import argparse
import json

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from pathway_common import PATH, RES, completeness_from_matrix, load_pathways, module_completeness
from redox_common import fit_models, safe_ap

FRACTIONS = (1.0, 0.9, 0.7, 0.5)
DRAWS = 3
SEED = 20260918


def protein_ko_table():
    """Long table of significant (genome, protein, KO) hits: the atoms that subsampling removes."""
    h = pd.read_parquet(PATH / "togo_kofam.parquet", columns=["genome_id", "protein", "KO", "significant"])
    return h[h.significant].drop_duplicates(["genome_id", "protein", "KO"])


def degraded_matrix(hits, genome_ids, columns, frac, seed):
    """KO presence/absence after keeping a random `frac` of each genome's proteins."""
    if frac >= 1.0:
        keep = hits
    else:
        rng = np.random.default_rng(seed)
        prot = hits[["genome_id", "protein"]].drop_duplicates()
        # hash-free but reproducible: one uniform draw per protein, thresholded at frac
        mask = rng.random(len(prot)) < frac
        keep = hits.merge(prot[mask], on=["genome_id", "protein"])
    m = pd.crosstab(keep.genome_id, keep.KO).clip(upper=1).astype(np.float32)
    return m.reindex(index=genome_ids, columns=columns).fillna(0).to_numpy(dtype=np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="PB2")
    ap.add_argument("--level", default="genus")
    ap.add_argument("--jobs", type=int, default=6)
    args = ap.parse_args()
    lab, traits, Y = load_pathways()
    oof = pd.read_parquet(RES / f"pathways/oof/{args.arm}.parquet")
    oof = oof[oof.level == args.level]

    hits = protein_ko_table()
    full = pd.read_csv(PATH / "ko_matrix.tsv", sep="\t", index_col=0)
    columns = list(full.columns)
    Xfull = full.reindex(lab.genome_id).fillna(0).to_numpy(dtype=np.float32)
    if args.arm == "PB3":
        MC = module_completeness(lab.genome_id).to_numpy()
        Xfull = np.hstack([Xfull, MC])
    fold = lab[f"fold_{args.level}"].to_numpy()

    def one(frac, draw):
        Xd = degraded_matrix(hits, lab.genome_id, columns, frac, SEED + 1000 * draw + int(frac * 100))
        if args.arm == "PB3":  # modules recomputed from the surviving KOs, not carried over from the complete genome
            Xd = np.hstack([Xd, completeness_from_matrix(Xd, columns)])
        rows = []
        for j, t in enumerate(traits):
            dt = oof[oof.trait == t]
            for k in sorted(dt.fold.unique()):
                cfg = json.loads(dt[dt.fold == k].config.iloc[0])
                if cfg.get("model") not in ("logreg", "xgb"):
                    continue
                tr, te = np.where(fold != k)[0], np.where(fold == k)[0]
                ms = fit_models(cfg, Xfull[tr], Y[tr, j])
                p = np.mean([m.predict_proba(Xd[te])[:, 1] for m in ms], axis=0)
                rows.append(dict(arm=args.arm, level=args.level, trait=t, fold=int(k), completeness=frac,
                                 draw=draw, n_pos=int(Y[te, j].sum()), auprc=safe_ap(Y[te, j], p)))
        return rows

    jobs = [(f, d) for f in FRACTIONS for d in (range(DRAWS) if f < 1.0 else [0])]
    out = Parallel(n_jobs=args.jobs)(delayed(one)(f, d) for f, d in jobs)
    df = pd.DataFrame([r for rs in out for r in rs])
    (RES / "pathways").mkdir(parents=True, exist_ok=True)
    df.to_csv(RES / "pathways/mag_robustness.csv", index=False, float_format="%.4f")
    # pooled over folds within a draw, then averaged over draws, then macro-averaged over traits
    per_trait = df.groupby(["trait", "completeness", "draw"]).auprc.mean().reset_index()
    macro = per_trait.groupby(["completeness", "draw"]).auprc.mean().reset_index()
    macro = macro.groupby("completeness").auprc.agg(["mean", "std"]).reset_index()
    macro.columns = ["completeness", "macro_auprc", "macro_auprc_sd_over_draws"]
    macro["delta_vs_complete"] = macro.macro_auprc - macro.macro_auprc[macro.completeness == 1.0].iloc[0]
    macro.to_csv(RES / "pathways/mag_robustness_macro.csv", index=False, float_format="%.4f")
    print(macro.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
