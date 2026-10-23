"""Run nested grouped CV for one arm (or the nulls) at both split levels.

usage: python src/06_run_arms.py {nulls,A,B1,B1r,B2,C,E,E_MAG} [--jobs N] [--split frozen|seed<k>|gtdb]
                                 [--label-mode main|block_nan]
Output: results/oof/<arm>.parquet for the primary analysis (split=frozen, label-mode=main), otherwise
        results/sensitivity/<split>__<label-mode>/oof/<arm>.parquet
"""
import argparse
import json

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from redox_common import (LEVELS, PROC, RES, RULE_GROUPS, TAX, fit_models, taxonomic_prior, load_features, load_labelled, nested_cv, per_trait,
                          prevalence_null, taxonomic_nn_null)

LR_C = [0.01, 0.1, 1.0, 10.0]
ARM_FEATURES = {"B1": "markers_provided", "B1r": "markers_recomputed", "B2": "energy_panel",
                "B2cnt": "energy_panel_log1p_counts",  # research loop R1
                "LG": "lithogenie", "B2LG": ("energy_panel", "lithogenie"),  # research loop R3
                "B2M": ("energy_panel", "metabolic_hmm")}  # hybrid A, docs/preregistration_hybrid_AB.md


def lr_configs():
    return [{"model": "logreg", "C": c} for c in LR_C]


def ko_configs():
    return lr_configs() + [{"model": "xgb", "max_depth": d} for d in (2, 4)]


def fit_predict(cfg, Xtr, ytr, Xte, seed=0):
    return np.mean([m.predict_proba(Xte)[:, 1] for m in fit_models(cfg, Xtr, ytr)], axis=0)


def rule_scores(lab, traits):
    """Per-genome rule score per trait = max over marker groups of the fraction of group KOs present."""
    X = load_features("markers_recomputed", lab.genome_id)
    S = np.zeros((len(lab), len(traits)), dtype=np.float32)
    for j, t in enumerate(traits):
        g = RULE_GROUPS[t]
        if g:
            S[:, j] = np.max([X[kos].mean(axis=1).to_numpy() for kos in g.values()], axis=0)
    return S


def run_rules(lab, traits, Y, level, S):
    """Arm A: no fitting. AUPRC uses the raw rule score; theta in {>0, 0.5, 1.0} maximising MCC on the outer training fold."""
    from sklearn.metrics import matthews_corrcoef

    fold = lab[f"fold_{level}"].to_numpy()
    rows = []
    for j, t in enumerate(traits):
        for k in sorted(np.unique(fold)):
            tr, te = np.where(fold != k)[0], np.where(fold == k)[0]
            ytr, s_tr = Y[tr, j], S[tr, j]
            # threshold chosen on training genomes only (no fitting, so inner folds = pooled training set)
            cands = {">0": 1e-6, "0.5": 0.5, "1.0": 1.0}
            mcc = {name: matthews_corrcoef(ytr, s_tr >= v) for name, v in cands.items()}
            best = max(mcc, key=lambda n: (mcc[n], -cands[n]))
            for i in te:
                rows.append(("A", level, t, lab.genome_id[i], int(k), Y[i, j], float(S[i, j]), cands[best],
                             json.dumps({"model": "rule", "theta": best}), np.nan))
    return pd.DataFrame(rows, columns=["arm", "level", "trait", "genome_id", "fold", "y", "prob", "threshold",
                                       "config", "inner_auprc"])


def oof_dir(split, label_mode):
    if split == "frozen" and label_mode == "main":
        return RES / "oof"
    return RES / "sensitivity" / f"{split}__{label_mode}" / "oof"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arm")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--split", default="frozen")
    ap.add_argument("--label-mode", default="main")
    args = ap.parse_args()
    lab, traits, Y = load_labelled(args.split, args.label_mode)
    out_dir = oof_dir(args.split, args.label_mode)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.arm == "nulls":
        out = [per_trait(f, lab, traits, Y, lvl) for lvl in LEVELS for f in (prevalence_null, taxonomic_nn_null)]
    elif args.arm == "A":
        S_all = rule_scores(lab, traits)
        out = []
        for lvl in LEVELS:
            for j, t in enumerate(traits):
                m = ~np.isnan(Y[:, j])
                out.append(run_rules(lab.loc[m].reset_index(drop=True), [t], Y[m][:, [j]], lvl, S_all[m][:, [j]]))
    else:
        if args.arm in ARM_FEATURES:
            names = ARM_FEATURES[args.arm]
            names = (names,) if isinstance(names, str) else names
            X = np.hstack([load_features(n, lab.genome_id).to_numpy() for n in names])
            X = X[:, X.std(axis=0) > 0]  # drop all-constant columns (label-free, unsupervised)
            configs, fp = ko_configs(), fit_predict
        elif args.arm == "C":
            # pooling is a hyper-parameter: both poolings side by side, the config selects the columns
            mean = load_features("esm2_mean", lab.genome_id).to_numpy()
            dim = mean.shape[1]
            X = np.hstack([mean, load_features("esm2_max", lab.genome_id).to_numpy()])
            configs = [{"model": "logreg", "C": c, "pooling": p} for p in ("mean", "max") for c in LR_C]

            def fp(cfg, Xtr, ytr, Xte, seed=0):
                sl = slice(0, dim) if cfg["pooling"] == "mean" else slice(dim, 2 * dim)
                return fit_predict({"model": "logreg", "C": cfg["C"]}, Xtr[:, sl], ytr, Xte[:, sl], seed)
        elif args.arm == "SEL":  # research loop R7: feature set {B2, LG, B2+LG} chosen per trait by inner CV
            b2 = load_features("energy_panel", lab.genome_id).to_numpy()
            lg = load_features("lithogenie", lab.genome_id).to_numpy()
            b2, lg = b2[:, b2.std(axis=0) > 0], lg[:, lg.std(axis=0) > 0]
            nb = b2.shape[1]
            X = np.hstack([b2, lg])
            cols = {"B2": slice(0, nb), "LG": slice(nb, X.shape[1]), "B2LG": slice(0, X.shape[1])}
            configs = [{**c, "features": f} for f in cols for c in ko_configs()]

            def fp(cfg, Xtr, ytr, Xte, seed=0):
                sl = cols[cfg["features"]]
                inner = {k: v for k, v in cfg.items() if k != "features"}
                return fit_predict(inner, Xtr[:, sl], ytr, Xte[:, sl], seed)
        elif args.arm in ("ENS", "ENSH"):  # R8: average of the same model on B2 and on LG features;
            # R11 (ENSH): second group = LithoGenie + multiheme cytochrome motif features
            b2 = load_features("energy_panel", lab.genome_id).to_numpy()
            lg = load_features("lithogenie", lab.genome_id).to_numpy()
            if args.arm == "ENSH":
                lg = np.hstack([lg, load_features("heme", lab.genome_id).to_numpy()])
            b2, lg = b2[:, b2.std(axis=0) > 0], lg[:, lg.std(axis=0) > 0]
            nb = b2.shape[1]
            X = np.hstack([b2, lg])
            configs = ko_configs()

            def fp(cfg, Xtr, ytr, Xte, seed=0):
                return 0.5 * (fit_predict(cfg, Xtr[:, :nb], ytr, Xte[:, :nb], seed) +
                              fit_predict(cfg, Xtr[:, nb:], ytr, Xte[:, nb:], seed))
        elif args.arm == "STACK":  # research loop R5: cross-trait stacking on B2 + LG features
            # level 1: fixed logistic model (C=0.1) per trait on B2 and on LG, averaged (as in ENS);
            # level 2: logistic regression on the 14 level-1 scores. Training genomes' level-1 scores are inner
            # out-of-fold; labels of other traits are looked up only for training rows (row-id column, not a feature).
            b2 = load_features("energy_panel", lab.genome_id).to_numpy()
            lg = load_features("lithogenie", lab.genome_id).to_numpy()
            b2, lg = b2[:, b2.std(axis=0) > 0], lg[:, lg.std(axis=0) > 0]
            nb, nl = b2.shape[1], lg.shape[1]
            X = np.hstack([b2, lg, np.arange(len(lab))[:, None]])  # last column = row id
            Yall = np.nan_to_num(Y, nan=0.0)
            obs_all = ~np.isnan(Y)
            L1 = {"model": "logreg", "C": 0.1}
            configs = [{"model": "stack", "C2": c} for c in (0.1, 1.0)]

            def level1(Xa, ya, Xb):
                ok = ~np.isnan(ya)
                if ya[ok].sum() == 0 or ya[ok].sum() == ok.sum():
                    return np.full(len(Xb), ya[ok].mean() if ok.any() else 0.0)
                return 0.5 * (fit_predict(L1, Xa[ok][:, :nb], ya[ok], Xb[:, :nb]) +
                              fit_predict(L1, Xa[ok][:, nb:nb + nl], ya[ok], Xb[:, nb:nb + nl]))

            def fp(cfg, Xtr, ytr, Xte, seed=0):
                rid = Xtr[:, -1].astype(int)
                Ytr = np.where(obs_all[rid], Yall[rid], np.nan)
                k = min(4, len(Xtr))
                folds = np.arange(len(Xtr)) % k
                Z_tr = np.zeros((len(Xtr), Y.shape[1]))
                Z_te = np.zeros((len(Xte), Y.shape[1]))
                for t in range(Y.shape[1]):
                    for f in range(k):
                        m = folds == f
                        Z_tr[m, t] = level1(Xtr[~m], Ytr[~m, t], Xtr[m])
                    Z_te[:, t] = level1(Xtr, Ytr[:, t], Xte)
                return fit_predict({"model": "logreg", "C": cfg["C2"]}, Z_tr, ytr, Z_te, seed)
        elif args.arm in ("B3", "B3ENS"):  # research loop R12: full significant-KO profile from TogoMCP (BH26 RDF)
            # 70/586 genomes are absent from BH26: their labels are treated as unobserved for this arm only;
            # comparisons with other arms are restricted to the common genomes (15_research_ledger.py)
            tm = pd.read_csv(PROC / "togo_ko_matrix_internal.tsv", sep="\t", index_col=0)
            log = pd.read_csv(PROC / "togo_fetch_log_internal.tsv", sep="\t")
            have = lab.genome_id.isin(set(log.loc[log.status == "ok", "genome_id"])).to_numpy()
            Y = Y.copy()
            Y[~have] = np.nan
            X = tm.reindex(lab.genome_id).fillna(0).to_numpy(dtype=np.float32)
            X = X[:, X[have].std(axis=0) > 0]
            configs = ko_configs()
            if args.arm == "B3":
                fp = fit_predict
            else:  # R12b B3ENS: average of B3 model and LG model (as ENS)
                lg = load_features("lithogenie", lab.genome_id).to_numpy()
                lg = lg[:, lg.std(axis=0) > 0]
                nb = X.shape[1]
                X = np.hstack([X, lg])

                def fp(cfg, Xtr, ytr, Xte, seed=0):
                    return 0.5 * (fit_predict(cfg, Xtr[:, :nb], ytr, Xte[:, :nb], seed) +
                                  fit_predict(cfg, Xtr[:, nb:], ytr, Xte[:, nb:], seed))
        elif args.arm in ("B2B1r", "B2RP", "B2XL"):
            # audit round 2 (M2): matched ensemble controls for the B2E fusion claim. Any equal-weight average of two
            # models may help at this N, so the PLM claim is only supported if B2E beats these:
            #   B2B1r: B2 + marker-subset model (no embeddings, same annotation source)
            #   B2RP : B2 + logistic probe on a random projection of the KO panel to 480 dims (embedding-shaped, no PLM)
            #   B2XL : B2 model with two learners averaged (logreg + xgb), no second feature block
            b2 = load_features("energy_panel", lab.genome_id).to_numpy()
            b2 = b2[:, b2.std(axis=0) > 0]
            nb = b2.shape[1]
            if args.arm == "B2B1r":
                other = load_features("markers_recomputed", lab.genome_id).to_numpy()
                other = other[:, other.std(axis=0) > 0]
            elif args.arm == "B2RP":
                rng = np.random.default_rng(20260918)  # fixed, label-free projection
                other = (b2 - b2.mean(0)) @ rng.standard_normal((nb, 480)) / np.sqrt(nb)
            else:
                other = np.zeros((len(lab), 0), dtype=np.float32)
            X = np.hstack([b2, other])
            configs = ko_configs()

            def fp(cfg, Xtr, ytr, Xte, seed=0):
                if args.arm == "B2XL":
                    a = fit_predict({"model": "logreg", "C": cfg.get("C", 0.1)}, Xtr[:, :nb], ytr, Xte[:, :nb], seed)
                    b = fit_predict({"model": "xgb", "max_depth": cfg.get("max_depth", 2)}, Xtr[:, :nb], ytr,
                                    Xte[:, :nb], seed)
                    return 0.5 * (a + b)
                ocfg = {"model": "logreg", "C": cfg.get("C", 0.1)}
                return 0.5 * (fit_predict(cfg, Xtr[:, :nb], ytr, Xte[:, :nb], seed) +
                              fit_predict(ocfg, Xtr[:, nb:], ytr, Xte[:, nb:], seed))
        elif args.arm in ("B2XLE", "B2XLC"):  # R16: learner-averaged B2 (B2XL) + PLM probe; tests PLM beyond ensembling
            b2 = load_features("energy_panel", lab.genome_id).to_numpy()
            b2 = b2[:, b2.std(axis=0) > 0]
            emb = load_features("bacformer_mean" if args.arm == "B2XLE" else "esm2_mean", lab.genome_id).to_numpy()
            nb = b2.shape[1]
            X = np.hstack([b2, emb])
            configs = ko_configs()

            def fp(cfg, Xtr, ytr, Xte, seed=0):
                c, d = cfg.get("C", 0.1), cfg.get("max_depth", 2)
                ko = 0.5 * (fit_predict({"model": "logreg", "C": c}, Xtr[:, :nb], ytr, Xte[:, :nb], seed) +
                            fit_predict({"model": "xgb", "max_depth": d}, Xtr[:, :nb], ytr, Xte[:, :nb], seed))
                return 0.5 * (ko + fit_predict({"model": "logreg", "C": c}, Xtr[:, nb:], ytr, Xte[:, nb:], seed))
        elif args.arm == "B2LGE":  # R14: average of B2 model, LithoGenie model and Bacformer-base probe
            b2 = load_features("energy_panel", lab.genome_id).to_numpy()
            lg = load_features("lithogenie", lab.genome_id).to_numpy()
            b2, lg = b2[:, b2.std(axis=0) > 0], lg[:, lg.std(axis=0) > 0]
            emb = load_features("bacformer_mean", lab.genome_id).to_numpy()
            n1, n2 = b2.shape[1], b2.shape[1] + lg.shape[1]
            X = np.hstack([b2, lg, emb])
            configs = ko_configs()

            def fp(cfg, Xtr, ytr, Xte, seed=0):
                ecfg = {"model": "logreg", "C": cfg.get("C", 0.1)}
                return (fit_predict(cfg, Xtr[:, :n1], ytr, Xte[:, :n1], seed) +
                        fit_predict(cfg, Xtr[:, n1:n2], ytr, Xte[:, n1:n2], seed) +
                        fit_predict(ecfg, Xtr[:, n2:], ytr, Xte[:, n2:], seed)) / 3
        elif args.arm in ("B2E_large", "B2XLE_large"):  # fusion with the protocol's preferred large encoder
            b2 = load_features("energy_panel", lab.genome_id).to_numpy()
            b2 = b2[:, b2.std(axis=0) > 0]
            emb = load_features("bacformer_large_mean", lab.genome_id).to_numpy()
            nb = b2.shape[1]
            X = np.hstack([b2, emb])
            configs = ko_configs()

            def fp(cfg, Xtr, ytr, Xte, seed=0):
                ecfg = {"model": "logreg", "C": cfg.get("C", 0.1)}
                if args.arm == "B2E_large":
                    ko = fit_predict(cfg, Xtr[:, :nb], ytr, Xte[:, :nb], seed)
                else:
                    ko = 0.5 * (fit_predict({"model": "logreg", "C": cfg.get("C", 0.1)}, Xtr[:, :nb], ytr,
                                            Xte[:, :nb], seed) +
                                fit_predict({"model": "xgb", "max_depth": cfg.get("max_depth", 2)}, Xtr[:, :nb], ytr,
                                            Xte[:, :nb], seed))
                return 0.5 * (ko + fit_predict(ecfg, Xtr[:, nb:], ytr, Xte[:, nb:], seed))
        elif args.arm in ("B2E", "B2C"):  # R13: average of the B2 model and a Bacformer-base (E) / ESM-2 (C) probe
            b2 = load_features("energy_panel", lab.genome_id).to_numpy()
            b2 = b2[:, b2.std(axis=0) > 0]
            emb = load_features("bacformer_mean" if args.arm == "B2E" else "esm2_mean", lab.genome_id).to_numpy()
            nb = b2.shape[1]
            X = np.hstack([b2, emb])
            configs = ko_configs()

            def fp(cfg, Xtr, ytr, Xte, seed=0):
                # the configuration applies to the KO block; the embedding block always uses the logistic probe
                # with the same C (or C=0.1 for XGBoost configurations): fixed before running
                ecfg = {"model": "logreg", "C": cfg.get("C", 0.1)}
                return 0.5 * (fit_predict(cfg, Xtr[:, :nb], ytr, Xte[:, :nb], seed) +
                              fit_predict(ecfg, Xtr[:, nb:], ytr, Xte[:, nb:], seed))
        elif args.arm == "B2tax":  # research loop R2: energy panel + leakage-free taxonomic prior
            X = load_features("energy_panel", lab.genome_id).to_numpy()
            X = X[:, X.std(axis=0) > 0]
            codes = np.stack([pd.factorize(lab[r])[0] for r in TAX], axis=1).astype(float)  # NaN -> -1
            n_feat = X.shape[1]
            X = np.hstack([X, codes])
            configs = None  # built per level below

            def fp(cfg, Xtr, ytr, Xte, seed=0):
                ex = TAX.index(cfg["exclude"])
                ctr, cte = Xtr[:, n_feat:].astype(int), Xte[:, n_feat:].astype(int)
                ptr = taxonomic_prior(ctr, ytr, ctr, ex, same_set=True)
                pte = taxonomic_prior(ctr, ytr, cte, ex)
                inner = {k: v for k, v in cfg.items() if k != "exclude"}
                return fit_predict(inner, np.hstack([Xtr[:, :n_feat], ptr[:, None]]), ytr,
                                   np.hstack([Xte[:, :n_feat], pte[:, None]]), seed)
        elif args.arm in ("E", "E_MAG", "E_large", "E_large_MAG"):
            feat = {"E": "bacformer_mean", "E_MAG": "bacformer_mag_mean", "E_large": "bacformer_large_mean",
                    "E_large_MAG": "bacformer_large_mag_mean"}[args.arm]
            X = load_features(feat, lab.genome_id).to_numpy()
            configs, fp = lr_configs(), fit_predict
        elif args.arm == "C_large":  # protocol arm C with the ESM-C 300M encoder (Bacformer-large's base model)
            mean = load_features("esmc_mean", lab.genome_id).to_numpy()
            dim = mean.shape[1]
            X = np.hstack([mean, load_features("esmc_max", lab.genome_id).to_numpy()])
            configs = [{"model": "logreg", "C": c, "pooling": p} for p in ("mean", "max") for c in LR_C]

            def fp(cfg, Xtr, ytr, Xte, seed=0):
                sl = slice(0, dim) if cfg["pooling"] == "mean" else slice(dim, 2 * dim)
                return fit_predict({"model": "logreg", "C": cfg["C"]}, Xtr[:, sl], ytr, Xte[:, sl], seed)
        else:
            raise SystemExit(f"unknown arm {args.arm}")

        def one(level, j):
            m = ~np.isnan(Y[:, j])
            cfgs = configs if configs is not None else [{**c, "exclude": level} for c in ko_configs()]
            return nested_cv(args.arm, X[m], lab.loc[m].reset_index(drop=True), [traits[j]], Y[m][:, [j]], cfgs,
                             fp, level)

        out = Parallel(n_jobs=args.jobs)(delayed(one)(lvl, j) for lvl in LEVELS for j in range(len(traits)))
    df = pd.concat(out, ignore_index=True)
    df.to_parquet(out_dir / f"{args.arm}.parquet", index=False)
    print(out_dir, df.groupby(["arm", "level"]).size().to_dict())


if __name__ == "__main__":
    main()
