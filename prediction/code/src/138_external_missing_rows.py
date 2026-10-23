"""External macro AUPRC for the Table 1 rows that had none: frequency baseline, curated marker rules, ESM-2 and
random forest (attention-based MIL needs per-protein ESM-2 embeddings of the external genomes; see below).

Same external genomes, categories, nitrification rule and bootstrap as src/22_external_eval.py. Models are fitted on
all 358 labeled genomes with their cross-validation choices and applied unchanged:
  frequency : training prevalence of each category (constant score)
  rules     : arm A, highest share of present KOs over the marker groups of a category (no fitting)
  ESM-2     : arm C, logistic regression on mean- or max-pooled ESM-2 genome vectors (pooling and C by inner CV)
  RF        : random forest on the 852-KO panel (47_trees_xai config grid, min_samples_leaf by inner CV)
usage: python src/138_external_missing_rows.py
Output: results/external_validation_missing_rows.csv
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import (RULE_GROUPS, best_mcc_threshold, fit_models, inner_folds, load_features, load_labelled,  # noqa
                          safe_ap)

ev = importlib.import_module("22_external_eval")
tx = importlib.import_module("47_trees_xai")
tx.TREE_JOBS = 8


def select(blocks, y, groups, fp):
    """blocks: {name: (Xint, Xext, configs)}; config and block chosen by inner grouped CV AUPRC."""
    folds = inner_folds(y, groups)
    best = None
    for name, (Xi, Xe, cfgs) in blocks.items():
        for cfg in cfgs:
            oof = np.zeros(len(y))
            for itr, iva in folds:
                oof[iva] = fp(cfg, Xi[itr], y[itr], Xi[iva])
            ap = safe_ap(y, oof)
            if best is None or ap > best[0]:
                best = (ap, name, cfg)
    _, name, cfg = best
    Xi, Xe, _ = blocks[name]
    return fp(cfg, Xi, y, Xe)


def lr_fp(cfg, Xtr, ytr, Xte):
    return np.mean([m.predict_proba(Xte)[:, 1] for m in fit_models(cfg, Xtr, ytr)], axis=0)


def rf_fp(cfg, Xtr, ytr, Xte):
    return tx.fit_predict(cfg, Xtr, ytr, Xte)


def main():
    lab, traits, Y = load_labelled()
    ext = pd.read_csv(ev.EXT / "labels.tsv", sep="\t").merge(pd.read_csv(ev.EXT / "lineage.tsv", sep="\t"),
                                                             on="genome_id", validate="1:1")
    ext = ext[ext.novel].reset_index(drop=True)
    groups = lab.group_genus.to_numpy()
    b2_int = load_features("energy_panel", lab.genome_id)
    keep = b2_int.columns[b2_int.std(axis=0) > 0]
    b2_ext = pd.read_csv(ev.EXT / "ko_matrix_panel.tsv", sep="\t", index_col=0).reindex(ext.genome_id)
    es = {p: (load_features(f"esm2_{p}", lab.genome_id),
              pd.read_parquet(ev.KAGGLE_EXT / f"genome_features_esm2_{p}.parquet").reindex(ext.genome_id))
          for p in ("mean", "max")}
    lr_cfgs = [{"model": "logreg", "C": c} for c in (0.01, 0.1, 1.0, 10.0)]
    scores = {a: {} for a in ("Frequency baseline", "Curated marker rules", "ESM-2", "Random forest")}
    for j, t in enumerate(traits):
        y = Y[:, j]
        scores["Frequency baseline"][t] = np.full(len(ext), np.nanmean(y))
        g = RULE_GROUPS[t]
        scores["Curated marker rules"][t] = np.max([b2_ext.reindex(columns=k).fillna(0).mean(axis=1).to_numpy()
                                                    for k in g.values()], axis=0) if g else np.zeros(len(ext))
        scores["ESM-2"][t] = select({p: (a.to_numpy(), b[a.columns].to_numpy(dtype=np.float32), lr_cfgs)
                                     for p, (a, b) in es.items()}, y, groups, lr_fp)
        scores["Random forest"][t] = select({"ko": (b2_int[keep].to_numpy(), b2_ext[keep].to_numpy(dtype=np.float32),
                                                    tx.CONFIGS["rf"])}, y, groups, rf_fp)
        print("fitted", t, flush=True)

    def s_for(arm, et):
        if et == "nitrification":
            return np.maximum(scores[arm]["donor_ammonia_oxidation"], scores[arm]["donor_nitrite_oxidation"])
        return scores[arm][et]
    genus = ext.ncbi_genus.fillna("NA").to_numpy()
    ug = np.unique(genus)
    members = {g: np.where(genus == g)[0] for g in ug}
    rng = np.random.default_rng(ev.SEED)
    boots = [np.concatenate([members[g] for g in rng.choice(ug, len(ug))]) for _ in range(ev.N_BOOT)]
    rows = []
    for arm in scores:
        aps, bs = [], []
        for et in ev.TRAITS_EXT:
            y_all = ext[et].to_numpy(dtype=float)
            m = ~np.isnan(y_all)
            if np.nansum(y_all[m]) < 10 or (y_all[m] == 0).sum() < 10:
                continue
            s = s_for(arm, et)
            aps.append(safe_ap(y_all[m], s[m]))
            bs.append([safe_ap(y_all[b][m[b]], s[b][m[b]]) for b in boots])
        bm = np.nanmean(np.array(bs), axis=0)
        rows.append(dict(arm=arm, n_traits=len(aps), macro_auprc=np.mean(aps), ci_low=np.nanpercentile(bm, 2.5),
                         ci_high=np.nanpercentile(bm, 97.5)))
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results/external_validation_missing_rows.csv", index=False, float_format="%.4f")
    pd.DataFrame([dict(genome_id=g, trait=et, arm=arm, score=s) for arm in scores for et in ev.TRAITS_EXT
                  for g, s in zip(ext.genome_id, s_for(arm, et))]).to_csv(
        ROOT / "results/external_validation_missing_rows_predictions.csv", index=False, float_format="%.6f")
    # random forest vs RedoxKO: paired genus bootstrap, all genomes and labels not derived from FAPROTAX
    b2 = pd.read_csv(ROOT / "results/external_validation_predictions.csv")
    b2 = b2[b2.arm == "B2"].pivot_table(index="genome_id", columns="trait", values="score")
    nofap = ~ext.faprotax_source.astype(bool).to_numpy()
    prs = []
    for sname, smask in (("all", np.ones(len(ext), bool)), ("no_faprotax_source", nofap)):
        d_pt, d_b = [], []
        for et in ev.TRAITS_EXT:
            y_all = ext[et].to_numpy(dtype=float)
            m = smask & ~np.isnan(y_all)
            if np.nansum(y_all[m]) < 10 or (y_all[m] == 0).sum() < 10:
                continue
            s_rf, s_b2 = s_for("Random forest", et), b2[et].reindex(ext.genome_id).to_numpy()
            d_pt.append(safe_ap(y_all[m], s_rf[m]) - safe_ap(y_all[m], s_b2[m]))
            d_b.append([safe_ap(y_all[b][m[b]], s_rf[b][m[b]]) - safe_ap(y_all[b][m[b]], s_b2[b][m[b]]) for b in boots])
        bm = np.nanmean(np.array(d_b), axis=0)
        prs.append(dict(subset=sname, n_traits=len(d_pt), rf_minus_redoxko=np.mean(d_pt), ci_low=np.nanpercentile(bm, 2.5),
                        ci_high=np.nanpercentile(bm, 97.5)))
    pr = pd.DataFrame(prs)
    pr.to_csv(ROOT / "results/external_validation_rf_vs_redoxko.csv", index=False, float_format="%.4f")
    print(pr.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(out.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
