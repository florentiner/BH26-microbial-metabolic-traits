"""Table 6 gaps: score extra-trees, the taxonomic baseline and the KEGG module rule on the pathway new-species test.

Same 99 species, 203 growth experiments and seven traits as src/84_novel_pathway_comparison.py (its saved
novel100_pathway_scores.csv supplies the tests). Each comparator is fitted on all 3,378 labelled pathway genomes,
with the choices it made in cross-validation:
* KEGG module rule (arm PA): score = highest completeness among the trait's mapped modules
  (pathway_common.TRAIT_MODULES); threshold chosen by MCC on the training genomes among >0, 0.5, 0.8, 1.0.
* Taxonomic baseline (null_taxNN): mean training label of the genomes sharing the lowest rank with the test
  genome (genus -> ... -> domain -> prevalence), lineages of the new species from data/additional/
  lineage_new_species.tsv; threshold = best MCC on the genus-level out-of-fold baseline scores.
* Extra-trees on the KO profile (47_trees_xai.py config grid, msl 1 or 5, 500 trees, 3 seeds, same prevalence
  filter): config by grouped (genus) inner CV AUPRC, threshold by MCC on the inner out-of-fold scores, refit on all.
The fine-tuned Bacformer model is not scored (GPU fine-tuning, no saved final model).

usage: python src/129_pathway_new_species_extra.py
Output: validation_new_organisms/novel100_pathway_extra.csv (+ _scores.csv)
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pathway_common import RES, TRAIT_MODULES, ko_matrix, load_pathways, module_completeness  # noqa: E402
from redox_common import best_mcc_threshold, inner_folds, safe_ap  # noqa: E402

V = ROOT / "validation_new_organisms"
TAX = ["domain", "phylum", "class", "order", "family", "genus"]
cmp84 = importlib.import_module("84_novel_pathway_comparison")
tx = importlib.import_module("47_trees_xai")
tx.TREE_JOBS = 10


def module_rule(lab, traits, Y, Xnew_mod):
    mc = module_completeness(lab.genome_id)
    out, thr = {}, {}
    for t in traits:
        j = traits.index(t)
        mods = [m for m in TRAIT_MODULES.get(t, []) if m in mc.columns]
        s_tr = mc[mods].to_numpy().max(axis=1) if mods else np.zeros(len(lab))
        cols = [f"mod_{m}" for m in mods if f"mod_{m}" in Xnew_mod.columns]
        out[t] = Xnew_mod[cols].max(axis=1) if cols else pd.Series(0.0, index=Xnew_mod.index)
        cands = {">0": 1e-6, "0.5": 0.5, "0.8": 0.8, "1.0": 1.0}
        m = {n: cmp84.mcc(Y[:, j], (s_tr >= v).astype(int)) for n, v in cands.items()}
        thr[t] = cands[max(m, key=lambda n: (m[n], -cands[n]))]
    return out, thr


def taxonomic(lab, traits, Y, lin):
    nulls = pd.read_parquet(RES / "pathways/oof/nulls.parquet")
    nulls = nulls[(nulls.arm == "null_taxNN") & (nulls.level == "genus")]
    out, thr = {}, {}
    for t in traits:
        y = pd.Series(Y[:, traits.index(t)])
        means = {r: y.groupby(lab[r].fillna("NA").to_numpy()).mean() for r in TAX}
        s = {}
        for acc, row in lin.iterrows():
            p = y.mean()
            for r in reversed(TAX):
                v = row.get(r)
                if pd.notna(v) and v in means[r].index:
                    p = means[r][v]
                    break
            s[acc] = p
        out[t] = pd.Series(s)
        n = nulls[nulls.trait == t]
        thr[t] = best_mcc_threshold(n.y.to_numpy(), n.prob.to_numpy())
    return out, thr


def extra_trees(lab, traits, Y, Xnew):
    X = ko_matrix(lab.genome_id)
    prev = X.mean()
    X = X.loc[:, (prev >= 0.01) & (prev <= 0.99)]
    X = X.loc[:, X.std(axis=0) > 0]
    cols = list(X.columns)
    Xa = X.to_numpy(dtype=np.float32)
    Xn = Xnew.reindex(columns=cols).fillna(0.0).to_numpy(dtype=np.float32)
    groups = lab.group_genus.to_numpy()
    out, thr, cfgs = {}, {}, {}
    for t in traits:
        y = Y[:, traits.index(t)].astype(int)
        folds = inner_folds(y, groups)
        best = None
        for cfg in tx.CONFIGS["et"]:
            oof = np.zeros(len(y))
            for itr, iva in folds:
                oof[iva] = tx.fit_predict(cfg, Xa[itr], y[itr], Xa[iva])
            ap = safe_ap(y, oof)
            if best is None or ap > best[0]:
                best = (ap, cfg, oof)
        _, cfg, oof = best
        thr[t], cfgs[t] = best_mcc_threshold(y, oof), cfg
        out[t] = pd.Series(tx.fit_predict(cfg, Xa, y, Xn), index=Xnew.index)
        print(t, cfg, f"inner AUPRC {best[0]:.3f}", flush=True)
    return out, thr


def main():
    d = pd.read_csv(V / "novel100_pathway_scores.csv")[["accession", "paper", "trait", "expected", "score_PB3",
                                                       "call_PB3"]]
    accs = sorted(d.accession.unique())
    traits_all = load_pathways()[1]
    lab, _, Y = load_pathways()
    traits = sorted(d.trait.unique())
    Yt = Y[:, [traits_all.index(t) for t in traits]]
    Xf = cmp84.features(accs)
    Xko = Xf[[c for c in Xf.columns if c.startswith("K") and c[1:].isdigit()]]
    Xmod = Xf[[c for c in Xf.columns if c.startswith("mod_")]]
    lin = pd.concat([pd.read_csv(ROOT / f"data/additional/{f}.tsv", sep="\t")
                     for f in ("lineage_new_species", "lineage_pathway_new_species_extra")]).set_index("accession")
    missing = set(accs) - set(lin.index)
    print(f"{len(accs)} species, {len(d)} tests, {len(traits)} traits; lineage missing for {len(missing)}")
    arms = {"rule": module_rule(lab, traits, Yt, Xmod), "tax": taxonomic(lab, traits, Yt, lin.reindex(accs))}
    prev = V / "novel100_pathway_extra_scores.csv"
    if "--reuse-et" in sys.argv and prev.exists():   # extra-trees is slow and does not use lineages
        p = pd.read_csv(prev).set_index(["accession", "trait"])
        d["score_ET"] = [p.score_ET[(a, t)] for a, t in zip(d.accession, d.trait)]
        d["call_ET"] = [p.call_ET[(a, t)] for a, t in zip(d.accession, d.trait)]
    else:
        arms["ET"] = extra_trees(lab, traits, Yt, Xko)
    for a, (sc, th) in arms.items():
        d[f"score_{a}"] = [sc[t].get(acc, np.nan) for acc, t in zip(d.accession, d.trait)]
        d[f"call_{a}"] = [int(s >= th[t]) for s, t in zip(d[f"score_{a}"], d.trait)]
    d.to_csv(V / "novel100_pathway_extra_scores.csv", index=False)
    rows = []
    for a in ["PB3", "ET", "tax", "rule"]:
        au, ap, nt = cmp84.ranking(d, f"score_{a}")
        lo, hi = cmp84.boot(d, lambda x: cmp84.ranking(x, f"score_{a}")[0])
        rows.append(dict(model=a, auroc=au, auroc_lo=lo, auroc_hi=hi, auprc=ap, traits=nt,
                         mcc=cmp84.mcc(d.expected, d[f"call_{a}"])))
    t = pd.DataFrame(rows)
    t.to_csv(V / "novel100_pathway_extra.csv", index=False, float_format="%.4f")
    print(t.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
