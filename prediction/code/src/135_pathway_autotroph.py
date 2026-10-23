"""Pathway models trained only on autotroph-like genomes, and their metrics on autotroph-like test genomes.

Autotroph-like: at least one carbon-fixation module or CO2-to-methane module >= 80% complete
(src/133_autotroph_subset.py CFIX), from the full KO profile.
Arms, each on the autotroph-like pathway genomes and the traits with >= 25 positives among them:
  EnergyKO-M_auto : the EnergyKO-M procedure (KO profile + module completeness, logistic regression / XGBoost, nested CV)
  ET_auto         : extra-trees on the KO profile (47_trees_xai config grid), the best learner in the full-data CV
Compared with the final EnergyKO-M on the same genomes and traits:
  cross-validation: macro AUPRC at genus and family level (frozen folds restricted to these genomes), paired
                    genus/family-cluster bootstrap of the difference;
  new species:      final models refitted on the autotroph-like genomes, macro AUROC and pooled MCC on the 203 growth
                    experiments (all species, and autotroph-like species only).
usage: python src/135_pathway_autotroph.py
Output: results/additional/pathway_autotroph/{cv.csv,new_species.csv}
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pathway_common import completeness_from_matrix, ko_matrix, load_pathways, module_completeness  # noqa: E402
from redox_common import PROC, nested_cv, safe_ap  # noqa: E402

arms = importlib.import_module("32_pathway_arms")
tf = importlib.import_module("35_pathway_train_final")
tx = importlib.import_module("47_trees_xai")
tx.TREE_JOBS = 4
a133 = importlib.import_module("133_autotroph_subset")
nov = importlib.import_module("84_novel_pathway_comparison")
OUT = ROOT / "results/additional/pathway_autotroph"
SEED, B = 20261021, 1000


def auto_mask(X):
    md = pd.read_csv(PROC / "kegg_modules_energy.tsv", sep="\t")
    md = md[md.module.isin(a133.CFIX)]
    return completeness_from_matrix(X.to_numpy(), X.columns, md).max(axis=1) >= 0.8


def et_fp(cfg, Xtr, ytr, Xte, seed=0):
    return tx.fit_predict(cfg, Xtr, ytr, Xte)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lab, traits, Y = load_pathways()
    K = ko_matrix(lab.genome_id)
    keep = auto_mask(K)
    labA, YA = lab[keep].reset_index(drop=True), Y[keep]
    use = [t for j, t in enumerate(traits) if YA[:, j].sum() >= 25]
    J = [traits.index(t) for t in use]
    print(f"autotroph-like pathway genomes {len(labA)} of {len(lab)}; traits {len(use)}: {use}", flush=True)
    XM = pd.concat([K[keep], module_completeness(labA.genome_id).add_prefix("mod_").set_index(K[keep].index)], axis=1)
    XM = XM.loc[:, XM.std(axis=0) > 0]
    XK = K[keep]
    prev = XK.mean()
    XK = XK.loc[:, (prev >= 0.01) & (prev <= 0.99)]
    XK = XK.loc[:, XK.std(axis=0) > 0]

    # ---------- cross-validation
    oof = {}
    for name, X, cfgs, fp, levels in [("EnergyKO-M_auto", XM, arms.configs(), arms.fit_predict, ("genus", "family")),
                                      ("ET_auto", XK, tx.CONFIGS["et"], et_fp, ("genus",))]:
        Xa = X.to_numpy(dtype=np.float32)
        outs = Parallel(n_jobs=4)(delayed(nested_cv)(name, Xa, labA, [t], YA[:, [j]], cfgs, fp, lvl, log=lambda s: None)
                                  for lvl in levels for t, j in zip(use, J))
        oof[name] = pd.concat(outs).groupby(["level", "trait", "genome_id"]).prob.mean()
        print(name, "CV done", flush=True)
    full = pd.read_parquet(ROOT / "results/pathways/oof/PB3.parquet")
    full = full[full.genome_id.isin(labA.genome_id)].groupby(["level", "trait", "genome_id"]).prob.mean()
    ids = labA.genome_id.to_numpy()
    rows = []
    for lvl in ("genus", "family"):
        g = labA[f"group_{lvl}"].to_numpy()
        ug = np.unique(g)
        rng = np.random.default_rng(SEED)
        boots = [np.concatenate([np.where(g == u)[0] for u in rng.choice(ug, len(ug))]) for _ in range(B)]

        def macro(src, idx):
            return np.nanmean([safe_ap(YA[idx, j], src.loc[(lvl, t)].reindex(ids[idx]).to_numpy()) for t, j in zip(use, J)
                               if 0 < YA[idx, j].sum() < len(idx)])
        allidx = np.arange(len(ids))
        ref = macro(full, allidx)
        rows.append(dict(level=lvl, model="EnergyKO-M (final, all genomes)", macro_auprc=ref))
        for name in oof:
            if (lvl, use[0]) not in oof[name].index:
                continue
            v = macro(oof[name], allidx)
            d = [macro(oof[name], i) - macro(full, i) for i in boots]
            rows.append(dict(level=lvl, model=name, macro_auprc=v, diff=v - ref, ci_low=np.percentile(d, 2.5),
                             ci_high=np.percentile(d, 97.5)))
    cv = pd.DataFrame(rows)
    cv.to_csv(OUT / "cv.csv", index=False, float_format="%.4f")
    print(cv.to_string(index=False, float_format=lambda v: f"{v:.3f}"), flush=True)

    # ---------- new species: refit on the autotroph-like genomes, score the 99 species
    accs = sorted(pd.read_csv(ROOT / "validation_new_organisms/novel100_pathway_extra_scores.csv").accession.unique())
    d = pd.read_csv(ROOT / "validation_new_organisms/novel100_pathway_extra_scores.csv")
    d = d[d.trait.isin(use)].reset_index(drop=True)
    Xv = nov.features(accs)
    Kv = Xv[[c for c in Xv.columns if c.startswith("K") and c[1:].isdigit()]]
    d["auto"] = d.accession.map(pd.Series(auto_mask(Kv.reindex(columns=K.columns).fillna(0)), index=Kv.index))
    sc, th = {}, {}
    for t, j in zip(use, J):
        if t not in set(d.trait):
            continue
        m, thr, cfg, inner, cal = tf.select_and_fit(XM, YA[:, j], labA.group_genus.to_numpy())
        Xa = Xv.reindex(columns=m.columns).fillna(0.0)
        sc[t], th[t] = pd.Series(np.asarray(m.predict_proba(Xa)), index=Xv.index), float(thr)
    d["score_auto"] = [sc[t][a] for a, t in zip(d.accession, d.trait)]
    d["call_auto"] = [int(sc[t][a] >= th[t]) for a, t in zip(d.accession, d.trait)]
    et_lab = labA.copy()
    et_sc, et_th = a129_extra_trees(et_lab, use, YA[:, J], Kv)
    d["score_ETauto"] = [et_sc[t][a] for a, t in zip(d.accession, d.trait)]
    d["call_ETauto"] = [int(et_sc[t][a] >= et_th[t]) for a, t in zip(d.accession, d.trait)]
    res = []
    for sub, s in [("all species", d), ("autotroph-like species", d[d.auto])]:
        for name, col in [("EnergyKO-M (final, all genomes)", "PB3"), ("EnergyKO-M_auto", "auto"), ("ET_auto", "ETauto")]:
            au = [roc_auc_score(g.expected, g[f"score_{col}"]) for _, g in s.groupby("trait") if 0 < g.expected.sum() < len(g)]
            res.append(dict(subset=sub, model=name, species=s.accession.nunique(), tests=len(s), traits=len(au),
                            macro_auroc=np.mean(au), mcc=matthews_corrcoef(s.expected, s[f"call_{col}"])))
    ns = pd.DataFrame(res)
    ns.to_csv(OUT / "new_species.csv", index=False, float_format="%.4f")
    d.to_csv(OUT / "new_species_scores.csv", index=False)
    print(ns.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


def a129_extra_trees(labA, use, YA, Kv):
    m129 = importlib.import_module("129_pathway_new_species_extra")
    m129.tx.TREE_JOBS = 4
    return m129.extra_trees(labA, use, YA, Kv)


if __name__ == "__main__":
    main()
