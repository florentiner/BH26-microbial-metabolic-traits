"""Refit every sequence-derivable arm on all 358 labelled genomes and save the model files.

For each arm and trait: hyper-parameters and the decision threshold are selected by grouped
(genus) CV on all labelled genomes with the same procedure used inside the outer folds, then the
model is refitted on all labelled genomes. Performance claims come from the outer CV
(results/metrics_*.csv), never from these refits.

Released arm (docs/preregistration.md): highest genus-level macro AUPRC among arms whose features
are computable from a genome sequence (A, B1r, B2, C, D, E). B1 consumes the provided marker table
and is not deployable on new genomes, so it is evaluated but not eligible.

Output: results/models/<arm>/model.joblib (dict trait -> fitted estimator/threshold/config) and
        results/models/<arm>/metadata.json; results/models/RELEASED.json
"""
import json
import sys
from datetime import date

import joblib
import numpy as np
import pandas as pd

from redox_common import (PROC, RES, ROOT, FittedFusion, FittedModel, best_mcc_threshold, fit_calibrator, fit_models,
                          inner_folds, load_features, load_labelled, safe_ap, sha256)

sys.path.insert(0, str(ROOT / "src"))
arms_mod = __import__("06_run_arms")
# docs/preregistration.md: sequence-derivable arms (B1 needs the provided marker table, so it is not deployable).
# The protocol's preferred encoders are the large ones, so arms C/E include their large variants (integrity audit M5).
PREREGISTERED = ["A", "B1r", "B2", "C", "D", "D_nopf", "E", "E_MAG", "C_large", "E_large", "E_large_MAG"]
RESEARCH_LOOP = ["B2E", "B2XL", "B2C", "B2XLE", "B2XLC", "B2LGE", "B2E_large", "B2XLE_large", "ENS"]
ELIGIBLE = PREREGISTERED + RESEARCH_LOOP
FEATURES = {"A": "rules on markers_recomputed", "B1r": "markers_recomputed", "B2": "energy_panel",
            "C": "esm2_mean | esm2_max (selected per trait)", "D": "ESM-2 protein bags (ABMIL)",
            "E": "bacformer_mean", "E_MAG": "bacformer_mag_mean",
            "B2E": "energy_panel + bacformer_mean (equal-weight fusion; research loop R13)",
            "B2XL": "energy_panel; logistic regression + XGBoost averaged (research loop R15)",
            "C_large": "esmc_mean | esmc_max (selected per trait)", "E_large": "bacformer_large_mean",
            "E_large_MAG": "bacformer_large_mag_mean", "D_nopf": "ESM-2 protein bags (ABMIL, no KO prefilter)"}


def select_and_fit(Xdf, y, groups, configs, column_fn):
    folds = inner_folds(y, groups)
    scores, oofs = {}, {}
    for ci, cfg in enumerate(configs):
        cols = column_fn(cfg)
        X = Xdf[cols].to_numpy(dtype=np.float32)
        oof = np.zeros(len(y))
        for itr, iva in folds:
            oof[iva] = np.mean([m.predict_proba(X[iva])[:, 1] for m in fit_models(cfg, X[itr], y[itr])], axis=0)
        scores[ci], oofs[ci] = safe_ap(y, oof), oof
    bi = max(scores, key=lambda c: (np.nan_to_num(scores[c], nan=-1), -c))
    cfg = configs[bi]
    cols = column_fn(cfg)
    model = FittedModel(cols, cfg, fit_models(cfg, Xdf[cols].to_numpy(dtype=np.float32), y))
    return model, best_mcc_threshold(y, oofs[bi]), cfg, scores[bi], fit_calibrator(y, oofs[bi])


def select_and_fit_xl(b2, y, groups):
    """B2XL: logistic regression (C from config, 0.1 for xgb configs) + XGBoost (depth from config, 2 for logreg
    configs) on the KO panel, averaged; identical to 06_run_arms.py arm B2XL."""
    folds = inner_folds(y, groups)
    X = b2.to_numpy(dtype=np.float32)

    def pred(ms, Z):
        return np.mean([m.predict_proba(Z)[:, 1] for m in ms], axis=0)

    configs = arms_mod.ko_configs()
    scores, oofs = {}, {}
    for ci, cfg in enumerate(configs):
        lr, xg = {"model": "logreg", "C": cfg.get("C", 0.1)}, {"model": "xgb", "max_depth": cfg.get("max_depth", 2)}
        oof = np.zeros(len(y))
        for itr, iva in folds:
            oof[iva] = 0.5 * (pred(fit_models(lr, X[itr], y[itr]), X[iva]) + pred(fit_models(xg, X[itr], y[itr]), X[iva]))
        scores[ci], oofs[ci] = safe_ap(y, oof), oof
    bi = max(scores, key=lambda c: (np.nan_to_num(scores[c], nan=-1), -c))
    cfg = configs[bi]
    lr, xg = {"model": "logreg", "C": cfg.get("C", 0.1)}, {"model": "xgb", "max_depth": cfg.get("max_depth", 2)}
    parts = [FittedModel(list(b2.columns), lr, fit_models(lr, X, y)),
             FittedModel(list(b2.columns), xg, fit_models(xg, X, y))]
    return FittedFusion(parts, cfg), best_mcc_threshold(y, oofs[bi]), cfg, scores[bi], fit_calibrator(y, oofs[bi])


def select_and_fit_xle_large(b2, emb, y, groups):
    """B2XLE_large (B2XL + Bacformer-large): ½ [logistic regression (C) + XGBoost (depth) on the KO panel, averaged] +
    ½ logistic probe (same C) on the Bacformer-large genome vector; identical to 06_run_arms.py arm B2XLE_large.
    Weights 0.25 / 0.25 / 0.5 are reproduced by nesting two equal-weight fusions."""
    folds = inner_folds(y, groups)
    Xb, Xe = b2.to_numpy(dtype=np.float32), emb.to_numpy(dtype=np.float32)

    def pred(ms, X):
        return np.mean([m.predict_proba(X)[:, 1] for m in ms], axis=0)

    configs = arms_mod.ko_configs()
    scores, oofs = {}, {}
    for ci, cfg in enumerate(configs):
        lr, xg = {"model": "logreg", "C": cfg.get("C", 0.1)}, {"model": "xgb", "max_depth": cfg.get("max_depth", 2)}
        oof = np.zeros(len(y))
        for itr, iva in folds:
            ko = 0.5 * (pred(fit_models(lr, Xb[itr], y[itr]), Xb[iva]) + pred(fit_models(xg, Xb[itr], y[itr]), Xb[iva]))
            oof[iva] = 0.5 * (ko + pred(fit_models(lr, Xe[itr], y[itr]), Xe[iva]))
        scores[ci], oofs[ci] = safe_ap(y, oof), oof
    bi = max(scores, key=lambda c: (np.nan_to_num(scores[c], nan=-1), -c))
    cfg = configs[bi]
    lr, xg = {"model": "logreg", "C": cfg.get("C", 0.1)}, {"model": "xgb", "max_depth": cfg.get("max_depth", 2)}
    ko_part = FittedFusion([FittedModel(list(b2.columns), lr, fit_models(lr, Xb, y)),
                            FittedModel(list(b2.columns), xg, fit_models(xg, Xb, y))], cfg)
    emb_part = FittedModel(list(emb.columns), lr, fit_models(lr, Xe, y))
    return FittedFusion([ko_part, emb_part], cfg), best_mcc_threshold(y, oofs[bi]), cfg, scores[bi], \
        fit_calibrator(y, oofs[bi])


def select_and_fit_fusion(b2, emb, y, groups):
    """B2E: config from ko_configs applies to the KO block; the embedding probe is logistic with the same C (or 0.1)."""
    folds = inner_folds(y, groups)
    Xb, Xe = b2.to_numpy(dtype=np.float32), emb.to_numpy(dtype=np.float32)

    def pred(ms, X):
        return np.mean([m.predict_proba(X)[:, 1] for m in ms], axis=0)

    scores, oofs = {}, {}
    configs = arms_mod.ko_configs()
    for ci, cfg in enumerate(configs):
        ecfg = {"model": "logreg", "C": cfg.get("C", 0.1)}
        oof = np.zeros(len(y))
        for itr, iva in folds:
            oof[iva] = 0.5 * (pred(fit_models(cfg, Xb[itr], y[itr]), Xb[iva]) +
                              pred(fit_models(ecfg, Xe[itr], y[itr]), Xe[iva]))
        scores[ci], oofs[ci] = safe_ap(y, oof), oof
    bi = max(scores, key=lambda c: (np.nan_to_num(scores[c], nan=-1), -c))
    cfg = configs[bi]
    ecfg = {"model": "logreg", "C": cfg.get("C", 0.1)}
    parts = [FittedModel(list(b2.columns), cfg, fit_models(cfg, Xb, y)),
             FittedModel(list(emb.columns), ecfg, fit_models(ecfg, Xe, y))]
    return FittedFusion(parts, cfg), best_mcc_threshold(y, oofs[bi]), cfg, scores[bi], fit_calibrator(y, oofs[bi])


def main():
    # --arms restricts the refit to a subset. XGBoost (tabular arms) and PyTorch (arm D) each load their own OpenMP
    # runtime; on macOS loading them in one process deadlocks, so arm D is refitted in its own invocation.
    only = sys.argv[sys.argv.index("--arms") + 1].split(",") if "--arms" in sys.argv else None
    if "--release-only" in sys.argv:  # recompute the release decision from existing metrics; refit nothing
        only = []
    lab, traits, Y = load_labelled()
    groups = lab.group_genus.to_numpy()
    macro = pd.read_csv(RES / "metrics_macro.csv")
    g = macro[(macro.split_level == "genus") & macro.arm.isin(ELIGIBLE)].sort_values("macro_auprc", ascending=False)
    # released arm (docs/deviations.md #12, #16): research-loop variants that improve on the pre-registered arm in
    # internal CV must also confirm on the external novel-genus set (external macro-AUPRC delta vs B2 not negative in
    # any subset). B2XL and ENS failed that test (external ENS - B2 = -0.008; B2XL - B2 = -0.007 overall and
    # -0.050 [-0.080, -0.014] without FAPROTAX-sourced labels), so the pre-registered B2 is released.
    extc = RES / "external_validation_contrasts.csv"
    ext = pd.read_csv(extc) if extc.exists() else pd.DataFrame(columns=["subset", "contrast", "trait", "ci_high"])
    ext = ext[ext.trait == "MACRO"]

    def externally_worse(arm):
        """Arms this one is significantly worse than on the external novel-genus set, in any reported subset."""
        e = ext[ext.contrast.str.startswith(f"{arm} - ")]
        return sorted({c.split(" - ")[1] for c, hi in zip(e.contrast, e.ci_high)
                       if hi < 0 and c.split(" - ")[1] != "null_taxNN"})

    # Released arm. Deviation #16 required a research-loop variant to confirm externally before displacing the
    # pre-registered arm. The same test is applied symmetrically to the pre-registered arms themselves (#21): the
    # released arm is the highest-ranked pre-registered arm on internal genus-level macro AUPRC that is NOT
    # significantly worse than another pre-registered arm on the independently curated external set. Without this,
    # a 0.001 internal difference - far inside the bootstrap CI - could displace a model that is 0.07 better on
    # independent labels.
    skipped = {}
    prereg_best = None
    for arm in g[g.arm.isin(PREREGISTERED)].arm:
        worse_than = [a for a in externally_worse(arm) if a in PREREGISTERED]
        if worse_than:
            skipped[arm] = worse_than
            print(f"  {arm}: internally ranked above {prereg_best or 'all'} but significantly worse externally "
                  f"than {worse_than}; not eligible for release", flush=True)
            continue
        prereg_best = arm
        break
    confirmed = []
    ec = ext[ext.contrast.str.endswith(f"- {prereg_best}")]
    for arm in set(a.split(" - ")[0] for a in ec.contrast):
        rows_ = ec[ec.contrast == f"{arm} - {prereg_best}"]
        if len(rows_) and (rows_.ci_high > 0).all() and (rows_.delta_auprc > 0).all():
            confirmed.append(arm)
    released = confirmed[0] if confirmed else prereg_best
    print("genus-level macro AUPRC (eligible arms):\n", g[["arm", "macro_auprc", "macro_auprc_ci_low",
                                                          "macro_auprc_ci_high"]].to_string())
    # research-loop fusions without a refit branch (B2C, B2XLE, B2XLC, B2LGE, ENS) are evaluated in the outer CV
    # and reported in results/research_ledger.csv, but are not refitted: none is released, and each is statistically
    # tied with the released arm, so a deployable artefact would not be used. Skipped explicitly, never silently.
    HAS_BRANCH = set(FEATURES) | {"A", "D", "D_nopf", "B2XL", "B2E", "C", "C_large", "B2XLE_large"}
    for arm in [a for a in ELIGIBLE if a in set(macro.arm) and (only is None or a in only)]:
        if arm not in HAS_BRANCH:
            print(f"skipping {arm}: no refit branch (CV-only research-loop variant)", flush=True)
            continue
        out = RES / "models" / arm
        out.mkdir(parents=True, exist_ok=True)
        bundle = {}
        if arm == "A":
            S = arms_mod.rule_scores(lab, traits)
            from sklearn.metrics import matthews_corrcoef
            for j, t in enumerate(traits):
                cands = {">0": 1e-6, "0.5": 0.5, "1.0": 1.0}
                mcc = {n: matthews_corrcoef(Y[:, j], S[:, j] >= v) for n, v in cands.items()}
                best = max(mcc, key=lambda n: (mcc[n], -cands[n]))
                bundle[t] = {"type": "rule", "groups": arms_mod.RULE_GROUPS[t], "threshold": cands[best]}
        elif arm == "D":
            d_mod = __import__("06b_arm_D_abmil")
            import torch
            bags = d_mod.load_bags(lab)
            all_idx = np.arange(len(lab))
            folds = inner_folds(Y[:, int(np.argmin(Y.sum(0)))], groups)
            score = {}
            for ci, cfg in enumerate(d_mod.CONFIGS):  # re-select the config on all labelled genomes
                oof = np.zeros_like(Y)
                for itr, iva in folds:
                    oof[iva], _ = d_mod.train_predict(bags[cfg["prefilter"]], Y, itr, iva, seed=ci)
                score[ci] = np.nanmean([safe_ap(Y[:, j], oof[:, j]) for j in range(len(traits))])
            cfg = d_mod.CONFIGS[max(score, key=lambda c: (score[c], -c))]
            thr = d_mod.ensemble_thresholds(bags[cfg["prefilter"]], Y, all_idx, folds)
            states = [d_mod.train_model(bags[cfg["prefilter"]], Y, all_idx, 100 + s).state_dict()
                      for s in range(d_mod.N_SEEDS_TEST)]
            torch.save({"state_dicts": states, "config": cfg, "traits": traits, "thresholds": thr,
                        "inner_cv_macro_auprc_all_labelled": score}, out / "model.pt")
            bundle = None
        elif arm == "B2XL":
            b2 = load_features("energy_panel", lab.genome_id)
            b2 = b2.loc[:, b2.std(axis=0) > 0]
            for j, t in enumerate(traits):
                model, thr, cfg, inner, cal = select_and_fit_xl(b2, Y[:, j], groups)
                bundle[t] = {"type": "estimator", "model": model, "threshold": thr, "config": cfg,
                             "inner_cv_auprc_all_labelled": inner, "calibrator": cal}
                print(arm, t, cfg, f"thr={thr:.3f}", f"inner AUPRC={inner:.3f}", flush=True)
        elif arm == "B2XLE_large":   # the paper's best internal model, saved fully trained for publication
            b2 = load_features("energy_panel", lab.genome_id)
            b2 = b2.loc[:, b2.std(axis=0) > 0]
            emb = load_features("bacformer_large_mean", lab.genome_id)
            for j, t in enumerate(traits):
                model, thr, cfg, inner, cal = select_and_fit_xle_large(b2, emb, Y[:, j], groups)
                bundle[t] = {"type": "estimator", "model": model, "threshold": thr, "config": cfg,
                             "inner_cv_auprc_all_labelled": inner, "calibrator": cal}
                print(arm, t, cfg, f"thr={thr:.3f}", f"inner AUPRC={inner:.3f}", flush=True)
        elif arm == "B2E":
            b2 = load_features("energy_panel", lab.genome_id)
            b2 = b2.loc[:, b2.std(axis=0) > 0]
            emb = load_features("bacformer_mean", lab.genome_id)
            for j, t in enumerate(traits):
                model, thr, cfg, inner, cal = select_and_fit_fusion(b2, emb, Y[:, j], groups)
                bundle[t] = {"type": "estimator", "model": model, "threshold": thr, "config": cfg,
                             "inner_cv_auprc_all_labelled": inner, "calibrator": cal}
                print(arm, t, cfg, f"thr={thr:.3f}", f"inner AUPRC={inner:.3f}", flush=True)
        else:
            if arm in ("C", "C_large"):
                pfx = "esm2" if arm == "C" else "esmc"
                Xdf = pd.concat([load_features(f"{pfx}_mean", lab.genome_id), load_features(f"{pfx}_max", lab.genome_id)],
                                axis=1)
                configs = [{"model": "logreg", "C": c, "pooling": p} for p in ("mean", "max") for c in arms_mod.LR_C]

                def column_fn(cfg, _cols=list(Xdf.columns), _pfx=pfx):
                    return [c for c in _cols if c.startswith(f"{_pfx}_{cfg['pooling']}_")]
            else:
                Xdf = load_features(FEATURES[arm], lab.genome_id)  # B1r, B2, E, E_MAG
                Xdf = Xdf.loc[:, Xdf.std(axis=0) > 0]
                configs = arms_mod.ko_configs() if arm in ("B1r", "B2") else arms_mod.lr_configs()

                def column_fn(cfg, _cols=list(Xdf.columns)):
                    return _cols
            for j, t in enumerate(traits):
                model, thr, cfg, inner, cal = select_and_fit(Xdf, Y[:, j], groups, configs, column_fn)
                bundle[t] = {"type": "estimator", "model": model, "threshold": thr, "config": cfg,
                             "inner_cv_auprc_all_labelled": inner, "calibrator": cal}
                print(arm, t, cfg, f"thr={thr:.3f}", f"inner AUPRC={inner:.3f}", flush=True)
        if bundle is not None:
            joblib.dump(bundle, out / "model.joblib")
        meta = {"arm": arm, "trained_on": "358 labelled genomes (all)", "date": str(date.today()),
                "traits": traits, "feature_set": FEATURES.get(arm, "energy_panel + bacformer_large_mean"),
                "kofam_release": "2026-08-25 (profiles.tar.gz Last-Modified)",
                "esm2_model": "facebook/esm2_t12_35M_UR50D@6fbf070e65b0b7291e7bbcd451118c216cff79d8",
                "bacformer_model": ("macwiatrak/bacformer-large-masked-complete-genomes@84f85c37a05e55559142b842c1569e265183b554"
                                    " (ESM-C 300M = Synthyra/ESMplusplus_small@f840ce88 protein vectors)"
                                    if arm in ("B2XLE_large", "E_large", "B2E_large") else
                                    "macwiatrak/bacformer-masked-complete-genomes@d1c461f1db75902e6502e5ba7588ea1a4ea96b29"),
                "labels_sha256": sha256(PROC / "labels.tsv"),
                "cv_performance_file": "results/metrics_per_trait.csv (outer CV; refit has no held-out estimate)"}
        (out / "metadata.json").write_text(json.dumps(meta, indent=2))
    ties = tied_arms(lab, traits, Y, released, [a for a in g.arm if a != released])
    (RES / "models" / "RELEASED.json").write_text(json.dumps(
        {"released_arm": released,
         "selection_rule": "highest-ranked pre-registered arm on internal genus-level macro AUPRC that is not "
                           "significantly worse than another pre-registered arm on the external novel-genus set "
                           "(docs/deviations.md #16, #21); a research-loop variant displaces it only if it also "
                           "confirms externally.",
         "internally_tied_top_set": sorted(g.arm.head(6)),
         "externally_confirmed_variants": confirmed,
         "preregistered_rule_choice": prereg_best,
         "excluded_by_external_check": skipped,
         "macro_auprc_genus": g.set_index("arm").macro_auprc.round(4).to_dict(),
         "paired_macro_delta_released_minus_arm_genus": ties,  # positive = released arm is ahead
         "tied_with_released": [a for a, v in ties.items() if v["ci_low"] <= 0 <= v["ci_high"]],
         "note": "arms whose 95% cluster-bootstrap CI of the macro-AUPRC difference includes 0 are statistically "
                 "tied with the released arm; the choice among tied arms is not a claim of superiority"},
        indent=2))
    print("released arm:", released)


def tied_arms(lab, traits, Y, ref, others, n_boot=1000, seed=20260917):
    """Paired genus-cluster bootstrap of macro AUPRC (released minus other) on the primary OOF predictions."""
    oof = pd.concat([pd.read_parquet(RES / "oof" / f"{a}.parquet") for a in [ref] + others], ignore_index=True)
    oof = oof[oof.level == "genus"]
    gid = {g: i for i, g in enumerate(lab.genome_id)}
    P = {}
    for a, d in oof.groupby("arm"):
        m = np.full(Y.shape, np.nan)
        for j, t in enumerate(traits):
            dt = d[d.trait == t]
            m[dt.genome_id.map(gid).to_numpy(), j] = dt.prob
        P[a] = m
    groups = lab.group_genus.to_numpy()
    ug = np.unique(groups)
    members = {g_: np.where(groups == g_)[0] for g_ in ug}
    rng = np.random.default_rng(seed)
    boots = [np.concatenate([members[g_] for g_ in rng.choice(ug, len(ug))]) for _ in range(n_boot)]

    def macro(a, idx):
        return np.nanmean([safe_ap(Y[idx, j], P[a][idx, j]) for j in range(len(traits))])

    out = {}
    for a in others:
        d = np.array([macro(ref, b) - macro(a, b) for b in boots])
        out[a] = {"delta": round(float(macro(ref, np.arange(len(Y))) - macro(a, np.arange(len(Y)))), 4),
                  "ci_low": round(float(np.percentile(d, 2.5)), 4), "ci_high": round(float(np.percentile(d, 97.5)), 4)}
    return out


if __name__ == "__main__":
    main()
