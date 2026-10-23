"""End-to-end neural pipeline (docs/preregistration_trees_xai.md, Addendum M):
genome -> Pyrodigal proteins -> frozen protein foundation model -> genome vector -> MLP head trained per trait.

Only the MLP head is trained; the foundation-model vectors were computed once (Kaggle GPU kernels, public proteomes).
Arms: FM_ESMC_MLP (ESM-C 300M mean+max), FM_BFL_MLP (Bacformer-large mean), FM_KO_BFL_MLP (852-KO panel +
Bacformer-large, one network over both).

usage: python src/59_fm_mlp.py cv --arm FM_BFL_MLP        # nested genus/family CV -> results/trees_xai/da/oof/<ARM>.parquet
       python src/59_fm_mlp.py external --arm FM_BFL_MLP  # refit on all 358, score the 562 lineage-novel genomes
       python src/59_fm_mlp.py summary                    # external table vs B2 (paired genus bootstrap)
"""
import argparse
import importlib
import os

import numpy as np
import pandas as pd

from redox_common import LEVELS, best_mcc_threshold, inner_folds, load_features, nested_cv, safe_ap

tx = importlib.import_module("47_trees_xai")
ROOT = tx.RES.parent
KAGGLE_EXT = ROOT / "data/kaggle/large_external"
OUT = tx.OUT / "da"
EXT_TRAITS = ["acceptor_sulfate", "acceptor_sulfur", "acceptor_CO2_methanogen_acetogen", "acceptor_iron",
              "acceptor_nitrate", "acceptor_oxygen", "donor_hydrogen_oxidation", "donor_sulfur_oxidation",
              "donor_iron_oxidation"]   # the 9 traits of results/external_validation_macro.csv (subset "all")
CONFIGS = [{"model": "fm_mlp", "hidden": h, "wd": w} for h in (64, 256) for w in (1e-3, 1e-1)]
ARMS = {"FM_ESMC_MLP": ["esmc_mean", "esmc_max"], "FM_BFL_MLP": ["bacformer_large_mean"],
        "FM_KO_BFL_MLP": ["energy_panel", "bacformer_large_mean"]}


def features(arm, ids_int, ids_ext=None):
    """Internal (and external) feature matrices for an arm, columns aligned; constant KO columns dropped."""
    Xi, Xe = [], []
    for name in ARMS[arm]:
        fi = load_features(name, ids_int)
        if name == "energy_panel":
            fi = fi.loc[:, fi.std(axis=0) > 0]
        Xi.append(fi)
        if ids_ext is not None:
            if name == "energy_panel":
                fe = pd.read_csv(ROOT / "data/external/ko_matrix_panel.tsv", sep="\t", index_col=0).reindex(ids_ext)
            else:
                fe = pd.read_parquet(KAGGLE_EXT / f"genome_features_{name}.parquet").reindex(ids_ext)
            fe = fe[fi.columns]
            assert not fe.isna().any().any(), f"external {name} features missing"
            Xe.append(fe)
    Xi = pd.concat(Xi, axis=1).to_numpy(dtype=np.float32)
    return Xi, (pd.concat(Xe, axis=1).to_numpy(dtype=np.float32) if ids_ext is not None else None)


def fit_predict(cfg, Xtr, ytr, Xte, seed=0):
    """MLP head (tx.NNClassifier: standardise -> Linear -> ReLU -> Dropout 0.3 -> Linear), 3 seeds averaged."""
    import torch
    torch.set_num_threads(int(os.environ.get("FM_THREADS", "2")))
    return np.mean([tx.NNClassifier(wd=cfg["wd"], seed=s + seed, hidden=cfg["hidden"]).fit(Xtr, ytr)
                    .predict_proba(Xte)[:, 1] for s in tx.SEEDS], axis=0)


def cmd_cv(arm, jobs):
    from joblib import Parallel, delayed
    lab, traits, Y, _ = tx.load_task("da")
    X, _ = features(arm, lab.genome_id)
    print(f"[{arm}] {len(lab)} genomes, {X.shape[1]} features, {len(traits)} traits", flush=True)
    frames = Parallel(n_jobs=jobs)(delayed(nested_cv)(arm, X, lab, [traits[j]], Y[:, [j]], CONFIGS, fit_predict, lvl)
                                   for lvl in LEVELS for j in range(len(traits)))
    (OUT / "oof").mkdir(parents=True, exist_ok=True)
    pd.concat(frames).to_parquet(OUT / "oof" / f"{arm}.parquet", index=False)
    print(f"[{arm}] written oof/{arm}.parquet", flush=True)


def cmd_external(arm):
    lab, traits, Y, _ = tx.load_task("da")
    ext = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t").merge(
        pd.read_csv(ROOT / "data/external/lineage.tsv", sep="\t"), on="genome_id", validate="1:1")
    ext = ext[ext.novel].reset_index(drop=True)
    Xi, Xe = features(arm, lab.genome_id, ext.genome_id)
    groups = lab.group_genus.to_numpy()
    rows = []
    for j, t in enumerate(traits):
        y = Y[:, j]
        folds = inner_folds(y, groups)
        best = None
        for ci, cfg in enumerate(CONFIGS):          # configuration chosen in inner grouped CV on all 358 genomes
            oof = np.full(len(y), np.nan)
            for itr, iva in folds:
                oof[iva] = fit_predict(cfg, Xi[itr], y[itr], Xi[iva], seed=ci)
            ap = safe_ap(y, oof)
            if best is None or ap > best[0]:
                best = (ap, ci, cfg, oof)
        _, ci, cfg, oof = best
        p = fit_predict(cfg, Xi, y, Xe, seed=ci)
        thr = best_mcc_threshold(y, oof)
        rows += [dict(genome_id=g, trait=t, arm=arm, score=float(s), threshold=thr) for g, s in zip(ext.genome_id, p)]
        print(f"[{arm}/external] {t} config {cfg}", flush=True)
    pd.DataFrame(rows).to_csv(OUT / f"external_{arm}.csv", index=False, float_format="%.5f")


def cmd_summary():
    ext = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t").merge(
        pd.read_csv(ROOT / "data/external/lineage.tsv", sep="\t"), on="genome_id", validate="1:1")
    ext = ext[ext.novel].set_index("genome_id")
    b2 = pd.read_csv(tx.RES / "external_validation_predictions.csv")
    b2 = b2[b2.arm == "B2"].pivot(index="genome_id", columns="trait", values="score")
    genus = ext.ncbi_genus.fillna("NA")
    ug = genus.unique()
    members = {g: np.where(genus.to_numpy() == g)[0] for g in ug}
    rng = np.random.default_rng(20260917)
    boots = [np.concatenate([members[g] for g in rng.choice(ug, len(ug))]) for _ in range(1000)]
    out = []
    arms = sorted(f.stem[len("external_"):] for f in OUT.glob("external_*.csv")
                  if not f.stem.endswith("_per_trait") and "summary" not in f.stem)   # FM MLP arms + fine-tuned arms
    for arm in arms:
        f = OUT / f"external_{arm}.csv"
        s = pd.read_csv(f).pivot(index="genome_id", columns="trait", values="score").reindex(ext.index)
        per, diffs = [], np.zeros(len(boots))
        for t in EXT_TRAITS:
            y = ext[t].to_numpy(dtype=float)
            ok = ~np.isnan(y)
            a, b = s[t].to_numpy(), b2[t].reindex(ext.index).to_numpy()
            per.append(dict(arm=arm, trait=t, n_pos=int(np.nansum(y)), auprc=safe_ap(y[ok], a[ok]),
                            auprc_B2=safe_ap(y[ok], b[ok])))
            for i, bi in enumerate(boots):
                m = bi[ok[bi]]
                diffs[i] += (safe_ap(y[m], a[m]) - safe_ap(y[m], b[m])) / len(EXT_TRAITS)
        per = pd.DataFrame(per)
        per.to_csv(OUT / f"external_{arm}_per_trait.csv", index=False, float_format="%.4f")
        out.append(dict(arm=arm, macro_auprc=per.auprc.mean(), macro_auprc_B2=per.auprc_B2.mean(),
                        delta_vs_B2=per.auprc.mean() - per.auprc_B2.mean(),
                        delta_ci_low=np.nanpercentile(diffs, 2.5), delta_ci_high=np.nanpercentile(diffs, 97.5)))
    d = pd.DataFrame(out)
    d.to_csv(OUT / "external_fm_mlp_summary.csv", index=False, float_format="%.4f")
    print(d.round(3).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["cv", "external", "summary"])
    ap.add_argument("--arm", choices=list(ARMS))
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    {"cv": lambda: cmd_cv(a.arm, a.jobs), "external": lambda: cmd_external(a.arm), "summary": cmd_summary}[a.cmd]()
