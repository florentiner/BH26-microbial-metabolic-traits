"""External validation on novel-genus genomes with Madin et al. labels (docs/external_validation.md).

Audit round 2 fixes: novelty is lineage-based (src/18b_external_lineage.py: accession, current NCBI genus, GTDB
genus); the taxonomic null uses current NCBI lineages (family -> order -> class -> phylum -> domain); paired
genus-cluster bootstrap contrasts are reported; B2E is included when external Bacformer-base features exist.

Models are fit on all 358 internal genomes (config and threshold by inner grouped CV, as in 09) and applied unchanged:
  B2   : energy-panel KO model (pre-registered)
  ENS  : average of the same model on B2 and LithoGenie features (research loop R8)
  B2XL : B2 KO panel with logistic regression and XGBoost averaged (research loop R15)
  B2E  : average of B2 model and Bacformer-base logistic probe (research loop R13)
  B2XLE_large / B2E_large : the same fusions with Bacformer-large (the protocol's preferred encoder)
  null_taxNN : internal label mean at the lowest rank shared with the external genome
Subsets: all novel genomes; species whose Madin data_source excludes FAPROTAX (taxonomy-inferred labels);
nitrate_strict (acceptor_nitrate positives only from denitrification / nitrate_reduction_to_ammonia).
Output: results/external_validation_per_trait.csv, results/external_validation_macro.csv,
        results/external_validation_contrasts.csv, results/external_validation_predictions.csv
"""
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from redox_common import ROOT, best_mcc_threshold, fit_models, inner_folds, load_features, load_labelled, safe_ap

EXT = ROOT / "data/external"
RES = ROOT / "results"
KAGGLE_EXT = ROOT / "data/kaggle/base_external"
KAGGLE_EXT_LARGE = ROOT / "data/kaggle/large_external"
TRAITS_EXT = ["acceptor_sulfate", "acceptor_sulfur", "acceptor_CO2_methanogen_acetogen", "acceptor_iron",
              "acceptor_nitrate", "acceptor_oxygen", "donor_hydrogen_oxidation", "donor_sulfur_oxidation",
              "donor_iron_oxidation", "nitrification"]
NULL_RANKS = ["family", "order", "class", "phylum", "domain"]
N_BOOT = 1000
SEED = 20260917


def ko_configs():
    return [{"model": "logreg", "C": c} for c in (0.01, 0.1, 1.0, 10.0)] + \
           [{"model": "xgb", "max_depth": d} for d in (2, 4)]


def predict(models, X):
    return np.mean([m.predict_proba(X)[:, 1] for m in models], axis=0)


def fit_select(blocks, y, groups):
    """blocks: list of (Xtrain, Xext, kind); kind 'ko' uses the config, 'emb' the logistic probe with the same C."""
    folds = inner_folds(y, groups)
    best = None

    def cfg_for(cfg, kind):
        if kind == "ko":
            return cfg
        if kind == "xgb":  # B2XL second learner (same definition as 06_run_arms.py)
            return {"model": "xgb", "max_depth": cfg.get("max_depth", 2)}
        return {"model": "logreg", "C": cfg.get("C", 0.1)}  # "emb" probe and B2XL "lr" learner

    w = np.array([b[3] if len(b) > 3 else 1.0 for b in blocks], dtype=float)
    w = w / w.sum()  # block weights reproduce the internal definitions (e.g. B2XLE = 0.25 lr + 0.25 xgb + 0.5 emb)

    for cfg in ko_configs():
        oof = np.zeros(len(y))
        for itr, iva in folds:
            oof[iva] = np.sum([wi * predict(fit_models(cfg_for(cfg, b[2]), b[0][itr], y[itr]), b[0][iva])
                               for wi, b in zip(w, blocks)], axis=0)
        ap = safe_ap(y, oof)
        if best is None or ap > best[0]:
            best = (ap, cfg, oof)
    ap, cfg, oof = best
    ext = np.sum([wi * predict(fit_models(cfg_for(cfg, b[2]), b[0], y), b[1]) for wi, b in zip(w, blocks)], axis=0)
    return ext, best_mcc_threshold(y, oof)


def main():
    lab, traits, Y = load_labelled()
    ext = pd.read_csv(EXT / "labels.tsv", sep="\t")
    lin = pd.read_csv(EXT / "lineage.tsv", sep="\t")
    ext = ext.merge(lin, on="genome_id", how="inner", validate="1:1")
    n_before = len(ext)
    ext = ext[ext.novel].reset_index(drop=True)
    print(f"novel external genomes: {len(ext)} (removed {n_before - len(ext)} non-novel)")
    groups = lab.group_genus.to_numpy()

    b2_int = load_features("energy_panel", lab.genome_id)
    keep = b2_int.columns[b2_int.std(axis=0) > 0]
    b2_ext = pd.read_csv(EXT / "ko_matrix_panel.tsv", sep="\t", index_col=0).reindex(ext.genome_id)[keep]
    lg_int = load_features("lithogenie", lab.genome_id)
    lg_keep = lg_int.columns[lg_int.std(axis=0) > 0]
    lg_ext = pd.read_parquet(EXT / "genome_features_lithogenie.parquet").reindex(ext.genome_id)[lg_keep]
    blocks = {"B2": [(b2_int[keep].to_numpy(), b2_ext.to_numpy(dtype=np.float32), "ko")]}
    blocks["ENS"] = blocks["B2"] + [(lg_int[lg_keep].to_numpy(), lg_ext.to_numpy(dtype=np.float32), "ko")]
    b2x = (b2_int[keep].to_numpy(), b2_ext.to_numpy(dtype=np.float32))
    blocks["B2XL"] = [(*b2x, "lr"), (*b2x, "xgb")]  # learner-ensembling control
    bf_file = KAGGLE_EXT / "genome_features_bacformer_mean.parquet"
    if bf_file.exists():
        bf_int = load_features("bacformer_mean", lab.genome_id)
        bf_ext = pd.read_parquet(bf_file).reindex(ext.genome_id)[bf_int.columns]
        assert not bf_ext.isna().any().any(), "external Bacformer features missing for some genomes"
        bfx = (bf_int.to_numpy(), bf_ext.to_numpy(dtype=np.float32))
        blocks["B2E"] = blocks["B2"] + [(*bfx, "emb")]
        blocks["B2XLE"] = [(*b2x, "lr", 0.25), (*b2x, "xgb", 0.25), (*bfx, "emb", 0.5)]
        blocks["B2LGE"] = blocks["ENS"] + [(*bfx, "emb")]
        es_file = KAGGLE_EXT / "genome_features_esm2_mean.parquet"
        if es_file.exists():
            es_int = load_features("esm2_mean", lab.genome_id)
            es_ext = pd.read_parquet(es_file).reindex(ext.genome_id)[es_int.columns]
            blocks["B2XLC"] = [(*b2x, "lr", 0.25), (*b2x, "xgb", 0.25),
                               (es_int.to_numpy(), es_ext.to_numpy(dtype=np.float32), "emb", 0.5)]
    else:
        print("external Bacformer-base features not available yet: B2E skipped")
    bfl_file = KAGGLE_EXT_LARGE / "genome_features_bacformer_large_mean.parquet"
    if bfl_file.exists():  # protocol's preferred encoder: the best internal variants
        bl_int = load_features("bacformer_large_mean", lab.genome_id)
        bl_ext = pd.read_parquet(bfl_file).reindex(ext.genome_id)[bl_int.columns]
        assert not bl_ext.isna().any().any(), "external Bacformer-large features missing for some genomes"
        blx = (bl_int.to_numpy(), bl_ext.to_numpy(dtype=np.float32))
        blocks["B2E_large"] = blocks["B2"] + [(*blx, "emb")]
        blocks["B2XLE_large"] = [(*b2x, "lr", 0.25), (*b2x, "xgb", 0.25), (*blx, "emb", 0.5)]
        es_l = KAGGLE_EXT_LARGE / "genome_features_esmc_mean.parquet"
        if es_l.exists():
            ec_int = load_features("esmc_mean", lab.genome_id)
            ec_ext = pd.read_parquet(es_l).reindex(ext.genome_id)[ec_int.columns]
            blocks["C_large"] = [(ec_int.to_numpy(), ec_ext.to_numpy(dtype=np.float32), "emb")]
            blocks["E_large"] = [(*blx, "emb")]

    scores = {a: {} for a in blocks}
    thr = {a: {} for a in blocks}
    scores["null_taxNN"] = {}
    for j, t in enumerate(traits):
        y = Y[:, j]
        for arm, bl in blocks.items():
            scores[arm][t], thr[arm][t] = fit_select(bl, y, groups)
        prior = np.full(len(ext), y.mean())
        for i, r in ext.iterrows():
            for rank in NULL_RANKS:
                v = r[f"ncbi_{rank}"]
                m = lab[rank].to_numpy() == v
                if pd.notna(v) and m.any():
                    prior[i] = y[m].mean()
                    break
        scores["null_taxNN"][t] = prior
        print("fitted", t, flush=True)
    used = []
    for _, r in ext.iterrows():
        used.append(next((rk for rk in NULL_RANKS if pd.notna(r[f"ncbi_{rk}"]) and (lab[rk] == r[f"ncbi_{rk}"]).any()),
                         "global"))
    print("taxonomic null rank used:", pd.Series(used).value_counts().to_dict())

    def s_for(arm, et):
        if et == "nitrification":
            return np.maximum(scores[arm]["donor_ammonia_oxidation"], scores[arm]["donor_nitrite_oxidation"])
        return scores[arm][et]

    genus = ext.ncbi_genus.fillna("NA").to_numpy()
    ug = np.unique(genus)
    members = {g: np.where(genus == g)[0] for g in ug}
    rng = np.random.default_rng(SEED)
    boots = [np.concatenate([members[g] for g in rng.choice(ug, len(ug))]) for _ in range(N_BOOT)]
    labels = {et: ext[et].to_numpy(dtype=float) for et in TRAITS_EXT}
    strict = ext.pathways.fillna("").str.contains("denitrification|nitrate_reduction_to_ammonia")
    labels_strict = labels.copy()
    labels_strict["acceptor_nitrate"] = np.where(strict, 1.0, np.where(ext.acceptor_nitrate == 1, np.nan, 0.0))
    subsets = {"all": (np.ones(len(ext), bool), labels),
               "no_faprotax_source": (~ext.faprotax_source.astype(bool).to_numpy(), labels),
               "nitrate_strict": (np.ones(len(ext), bool), labels_strict)}
    arms = list(scores)
    rows, contrasts, preds = [], [], []
    for sname, (smask, lbl) in subsets.items():
        B = {a: [] for a in arms}
        evaluable = []
        for et in TRAITS_EXT:
            if sname == "nitrate_strict" and et != "acceptor_nitrate":
                continue
            y_all = lbl[et]
            m = smask & ~np.isnan(y_all)
            n_pos, n_neg = int(np.nansum(y_all[m])), int((y_all[m] == 0).sum())
            n_pos_fap = int(np.nansum(y_all[m & ext.faprotax_source.astype(bool).to_numpy()]))
            if n_pos < 10 or n_neg < 10:
                rows.append(dict(subset=sname, trait=et, n_pos=n_pos, n_neg=n_neg, evaluable=False))
                continue
            evaluable.append(et)
            for arm in arms:
                s = s_for(arm, et)
                y, p = y_all[m], s[m]
                b = np.array([safe_ap(y_all[bi][m[bi]], s[bi][m[bi]]) for bi in boots])
                B[arm].append(b)
                row = dict(subset=sname, trait=et, arm=arm, evaluable=True, n=int(m.sum()), n_pos=n_pos, n_neg=n_neg,
                           n_pos_from_faprotax_species=n_pos_fap, prevalence=y.mean(), auprc=safe_ap(y, p),
                           auprc_ci_low=np.nanpercentile(b, 2.5), auprc_ci_high=np.nanpercentile(b, 97.5),
                           auroc=roc_auc_score(y, p))
                row["auprc_lift"] = row["auprc"] / row["prevalence"]
                if arm in thr and et != "nitrification":
                    row["mcc"] = matthews_corrcoef(y, p >= thr[arm][et])
                rows.append(row)
                if sname == "all":
                    preds += [dict(genome_id=g, trait=et, arm=arm, y=yy, score=pp)
                              for g, yy, pp in zip(ext.genome_id[m], y, p)]
        pt = pd.DataFrame([r for r in rows if r["subset"] == sname and r.get("evaluable")])
        for a, bref in (("B2", "null_taxNN"), ("ENS", "B2"), ("B2XL", "B2"), ("B2XL", "null_taxNN"), ("B2E", "B2"),
                        ("B2E", "B2XL"), ("B2E", "null_taxNN"), ("B2XLE", "B2"), ("B2XLC", "B2"), ("B2LGE", "B2"),
                        ("B2XLE_large", "B2"), ("B2XLE_large", "null_taxNN"), ("B2E_large", "B2"),
                        ("E_large", "B2"), ("E_large", "null_taxNN"), ("C_large", "B2")):
            if a not in B or bref not in B or not evaluable:
                continue
            dm = np.nanmean(np.stack(B[a]), axis=0) - np.nanmean(np.stack(B[bref]), axis=0)
            point = pt[pt.arm == a].auprc.mean() - pt[pt.arm == bref].auprc.mean()
            contrasts.append(dict(subset=sname, contrast=f"{a} - {bref}", trait="MACRO", n_traits=len(evaluable),
                                  delta_auprc=point, ci_low=np.nanpercentile(dm, 2.5),
                                  ci_high=np.nanpercentile(dm, 97.5), boot_frac_le0=float(np.nanmean(dm <= 0))))
            for k, et in enumerate(evaluable):
                d = B[a][k] - B[bref][k]
                contrasts.append(dict(subset=sname, contrast=f"{a} - {bref}", trait=et, n_traits=1,
                                      delta_auprc=pt[(pt.arm == a) & (pt.trait == et)].auprc.iloc[0] -
                                      pt[(pt.arm == bref) & (pt.trait == et)].auprc.iloc[0],
                                      ci_low=np.nanpercentile(d, 2.5), ci_high=np.nanpercentile(d, 97.5),
                                      boot_frac_le0=float(np.nanmean(d <= 0))))
    df = pd.DataFrame(rows)
    df.to_csv(RES / "external_validation_per_trait.csv", index=False, float_format="%.4f")
    cdf = pd.DataFrame(contrasts)
    # multiplicity control over the external contrast family (audit round 3, M3): two-sided bootstrap p from the
    # paired delta distribution, BH-adjusted within each (subset, level) family
    cdf["p_boot"] = cdf.apply(lambda r: min(1.0, 2 * min(r.boot_frac_le0, 1 - r.boot_frac_le0)), axis=1)
    for scope in cdf.subset.unique():
        for lvl in ("MACRO", "trait"):
            m_ = (cdf.subset == scope) & ((cdf.trait == "MACRO") if lvl == "MACRO" else (cdf.trait != "MACRO"))
            p_ = cdf.loc[m_, "p_boot"].to_numpy()
            order = np.argsort(p_)
            ranked = np.minimum.accumulate((p_[order] * len(p_) / np.arange(1, len(p_) + 1))[::-1])[::-1]
            q = np.empty_like(p_)
            q[order] = np.minimum(ranked, 1)
            cdf.loc[m_, "q_bh"] = q
    cdf.drop(columns=["boot_frac_le0"]).to_csv(RES / "external_validation_contrasts.csv", index=False,
                                               float_format="%.4f")
    pd.DataFrame(preds).to_csv(RES / "external_validation_predictions.csv", index=False, float_format="%.4f")
    ev = df[df.evaluable == True]  # noqa: E712
    mac = ev.groupby(["subset", "arm"]).agg(n_traits=("trait", "nunique"), macro_auprc=("auprc", "mean"),
                                            macro_auroc=("auroc", "mean"), macro_lift=("auprc_lift", "mean"),
                                            macro_mcc=("mcc", "mean")).reset_index()
    mac.to_csv(RES / "external_validation_macro.csv", index=False, float_format="%.4f")
    print(ev.pivot_table(index=["subset", "trait"], columns="arm", values="auprc").round(3).to_string())
    print(mac.round(4).to_string())
    c = pd.DataFrame(contrasts)
    print(c[c.trait == "MACRO"].round(4).to_string())


if __name__ == "__main__":
    main()
