"""Shared data loading, nested grouped cross-validation, nulls and metrics."""
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data/processed"
RES = ROOT / "results"
LEVELS = ("genus", "family")
N_INNER = 4
INNER_SEED = 17
TAX = ["domain", "phylum", "class", "order", "family", "genus"]

# Marker-gene groups per trait for the rule arm (A). Blocks follow the column order of the
# provided marker table (donor block, then acceptor block). Fixed before scoring. Markers that
# are not specific to the trait are excluded because a max-over-groups rule would fire on almost every genome:
# cynS (cyanase), rbcL/prkB (CO2 fixation, not an electron donor), nirBD K00362/K00363
# (assimilatory), hdrABC, Rnf, pta (widespread outside methanogens/acetogens); cooS is kept only under CO oxidation.
RULE_GROUPS = {
    "donor_sulfur_oxidation": {
        "sqr": ["K17218"], "fccAB": ["K17229", "K17230"],
        "soxABXYZ": ["K17222", "K17223", "K17224", "K17225", "K17226", "K17227"],
        "dsrAB": ["K11180", "K11181"], "aprAB_sat": ["K00394", "K00395", "K00958"], "tsdA": ["K19713"],
    },
    "donor_ammonia_oxidation": {"amoABC": ["K10944", "K10945", "K10946"], "hao": ["K10535"]},
    "donor_nitrite_oxidation": {"nxrAB": ["K00370", "K00371"]},
    "donor_hydrogen_oxidation": {"NiFe_hyaB_hybC": ["K06281", "K05922"], "FeFe_hydA": ["K00532", "K00533"]},
    "donor_CO_oxidation": {"coxLMS": ["K03520", "K03519", "K03518"], "cooS": ["K00198"]},
    "donor_iron_oxidation": {},
    "donor_other": {
        "fdhF_fdoG": ["K00123", "K00122"], "pmoABC": ["K10944", "K10945", "K10946"], "mmoX": ["K16157"],
        "mxaFI": ["K14028", "K14029", "K23995"],
    },
    "acceptor_oxygen": {
        "coxABC": ["K02274", "K02275", "K02276"], "ccoNOQP": ["K00404", "K00405", "K00406"],
        "cydAB": ["K00425", "K00426"],
    },
    "acceptor_nitrate": {
        "narGHJI": ["K00370", "K00371", "K00372", "K00373", "K00374"], "napAB": ["K02567", "K02568"],
        "nirK_nirS": ["K00368", "K15864"], "norBC": ["K04561", "K02305"], "nosZ": ["K00376"],
        "nrfAH": ["K03385", "K15876"],
    },
    "acceptor_sulfate": {"sat": ["K00958"], "aprAB": ["K00394", "K00395"], "dsrAB": ["K11180", "K11181"]},
    "acceptor_sulfur": {"asrABC": ["K16950", "K16951", "K16952"], "sreABC": ["K17219"], "psrABC": ["K08352"]},
    "acceptor_CO2_methanogen_acetogen": {
        "mcrABG": ["K00399", "K00401", "K00402"], "acsB": ["K14138"],
    },
    "acceptor_iron": {},
    "acceptor_other": {"hzsABC": ["K20932", "K20933", "K20934", "K20935"]},
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_labelled(split="frozen", label_mode="main"):
    """Labelled genomes, trait list and label matrix Y.

    split      : "frozen" (primary, hash-checked NCBI-taxonomy folds) | "seed<k>" (same NCBI groups, fold
                 allocation re-drawn with seed k) | "gtdb" (groups and taxonomy from GTDB R232)
    label_mode : "main" (closed world: 0 = negative) | "block_nan" (a genome's donor or acceptor block that is
                 entirely 0 is treated as unobserved -> NaN)
    Y may contain NaN only in block_nan mode.
    """
    labels = pd.read_csv(PROC / "labels.tsv", sep="\t")
    splits = pd.read_csv(PROC / "splits.tsv", sep="\t")
    expected = (PROC / "splits.sha256").read_text().split()[0]
    assert sha256(PROC / "splits.tsv") == expected, "splits.tsv changed after freezing"
    traits = [c for c in labels.columns if c.startswith(("donor_", "acceptor_"))]
    lab = splits.merge(labels.drop(columns=["genus", "family"]), on="genome_id", how="left", validate="1:1")
    assert lab.labelled.all() and not lab[traits].isna().any().any()
    Y = lab[traits].to_numpy(dtype=float)
    if split != "frozen":
        lab = resplit(lab, Y, split)
    if label_mode == "block_nan":
        for prefix in ("donor_", "acceptor_"):
            cols = [j for j, t in enumerate(traits) if t.startswith(prefix)]
            empty = Y[:, cols].sum(axis=1) == 0
            Y[np.ix_(empty, cols)] = np.nan
    elif label_mode == "curated_D":  # docs/preregistration_label_recuration_D.md: NaN = unknown after re-curation
        cur = pd.read_csv(PROC / "labels_curated_D.tsv", sep="\t").set_index("genome_id")
        Y = cur.loc[lab.genome_id, traits].to_numpy(dtype=float)
    elif label_mode != "main":
        raise ValueError(label_mode)
    return lab, traits, Y


def resplit(lab, Y, split):
    """Sensitivity partitions (not used for the primary results). Written to data/processed/splits_<split>.tsv."""
    import importlib
    ls = importlib.import_module("05_labels_splits")
    lab = lab.copy()
    if split.startswith("seed"):
        seed = int(split[4:])
    elif split == "random":  # context only: no phylogenetic grouping (each genome is its own group)
        seed = ls.SEED
        lab["group_genus"] = lab.genome_id
        lab["group_family"] = lab.genome_id
    elif split == "gtdb":
        seed = ls.SEED
        g = pd.read_csv(PROC / "gtdb_taxonomy.tsv", sep="\t").set_index("genome_id")
        in_gtdb = lab.genome_id.map(g.gtdb_genus).notna()
        for r in TAX:  # genomes absent from GTDB keep NCBI names, prefixed so they never merge with GTDB taxa
            gt = lab.genome_id.map(g[f"gtdb_{r}"])
            lab[r] = np.where(in_gtdb, gt, ("ncbi:" + lab[r].astype(str)).where(lab[r].notna(), np.nan))
            lab[r] = lab[r].replace("", np.nan)
        lab["group_genus"] = lab.genus.astype(str)
        lab["group_family"] = np.where(lab.family.isna(), "genus:" + lab.genus.astype(str),
                                       "family:" + lab.family.astype(str))
    else:
        raise ValueError(split)
    lab["fold_genus"] = ls.greedy_group_folds(lab.group_genus.to_numpy(), Y, ls.N_FOLDS, seed)
    lab["fold_family"] = ls.greedy_group_folds(lab.group_family.to_numpy(), Y, ls.N_FOLDS, seed)
    out = PROC / f"splits_{split}.tsv"
    lab[["genome_id", "genus", "family", "group_genus", "group_family", "fold_genus", "fold_family"]].to_csv(
        out, sep="\t", index=False)
    digest = sha256(out)  # deterministic given the seed; hashed so sensitivity splits are tamper-evident too
    hashes = PROC / "splits_sensitivity.sha256"
    lines = [ln for ln in (hashes.read_text().splitlines() if hashes.exists() else []) if not ln.endswith(out.name)]
    prev = [ln for ln in lines if ln.split("  ")[-1] == out.name]
    assert not prev or prev[0].startswith(digest), f"{out.name} changed since first written"
    hashes.write_text("\n".join(sorted(lines + [f"{digest}  {out.name}"])) + "\n")
    return lab


def per_trait(fn, lab, traits, Y, *args, **kw):
    """Run fn(lab_sub, [trait], Y_sub, ...) per trait on genomes whose label for that trait is observed."""
    out = []
    for j, t in enumerate(traits):
        m = ~np.isnan(Y[:, j])
        out.append(fn(lab.loc[m].reset_index(drop=True), [t], Y[m][:, [j]], *args, **kw))
    return pd.concat(out, ignore_index=True)


def load_features(name, genome_ids):
    """Genome x feature matrix (DataFrame indexed by genome_id, rows in the given order)."""
    if name == "markers_provided":
        df = pd.read_csv(ROOT / "ko_presence_586genomes.tsv", sep="\t", encoding="utf-8-sig").set_index("genome_id")
        cols = [c for c in df.columns if re.match(r"K\d{5} ", c)]
        X = df[cols]
        X.columns = [c.split()[0] for c in cols]
        X = X.T.groupby(level=0).max().T  # duplicated KO columns are identical (verified)
    elif name in ("markers_recomputed", "energy_panel"):
        X = pd.read_csv(PROC / "ko_matrix_panel.tsv", sep="\t", index_col=0)
        if name == "markers_recomputed":
            panel = pd.read_csv(PROC / "ko_panel.tsv", sep="\t")
            X = X[sorted(panel.loc[panel.is_marker, "KO"])]
    elif name == "energy_panel_log1p_counts":
        X = np.log1p(pd.read_csv(PROC / "ko_counts_panel.tsv", sep="\t", index_col=0))
    elif (PROC / f"genome_features_{name}.parquet").exists():  # esm2_*, bacformer_*, esmc_*, lithogenie, ...
        X = pd.read_parquet(PROC / f"genome_features_{name}.parquet")
    else:
        raise KeyError(name)
    missing = set(genome_ids) - set(X.index)
    assert not missing, f"{name}: missing genomes {sorted(missing)[:5]}"
    return X.loc[list(genome_ids)].astype(np.float32)


def inner_folds(y, groups, seed=INNER_SEED):
    """Grouped, class-stratified inner folds; falls back to fewer splits for rare traits."""
    n_pos_groups = len(set(groups[y == 1]))
    n = min(N_INNER, max(2, n_pos_groups))
    return list(StratifiedGroupKFold(n_splits=n, shuffle=True, random_state=seed).split(np.zeros(len(y)), y, groups))


def safe_ap(y, p):
    return average_precision_score(y, p) if 0 < y.sum() < len(y) else np.nan


def best_mcc_threshold(y, p):
    """Threshold maximising MCC on (inner out-of-fold) predictions."""
    cands = np.unique(p)
    if len(cands) > 200:
        cands = np.quantile(p, np.linspace(0, 1, 201))
    best_t, best_m = 0.5, -np.inf
    for t in cands:
        m = matthews_corrcoef(y, p >= t)
        if m > best_m:
            best_t, best_m = float(t), m
    return best_t


def nested_cv(arm, X, lab, traits, Y, configs, fit_predict, level, log=print):
    """Outer frozen folds x inner grouped CV for config + threshold selection.

    fit_predict(config, Xtr, ytr, Xte, seed) -> probabilities for Xte.
    Returns a long DataFrame of out-of-fold predictions for this arm and level.
    """
    fold = lab[f"fold_{level}"].to_numpy()
    groups = lab[f"group_{level}"].to_numpy()
    rows = []
    for j, trait in enumerate(traits):
        y = Y[:, j]
        for k in sorted(np.unique(fold)):
            tr, te = np.where(fold != k)[0], np.where(fold == k)[0]
            ytr = y[tr]
            if ytr.sum() == 0 or ytr.sum() == len(ytr):  # degenerate training fold
                p_te = np.full(len(te), ytr.mean())
                best, thr, inner_score = "constant", 0.5, np.nan
            else:
                folds = inner_folds(ytr, groups[tr])
                scores = {}
                oof_by_cfg = {}
                for ci, cfg in enumerate(configs):
                    oof = np.full(len(tr), np.nan)
                    for itr, iva in folds:
                        yi = ytr[itr]
                        if yi.sum() == 0 or yi.sum() == len(yi):
                            oof[iva] = yi.mean()
                        else:
                            oof[iva] = fit_predict(cfg, X[tr][itr], yi, X[tr][iva], seed=ci)
                    oof_by_cfg[ci] = oof
                    scores[ci] = safe_ap(ytr, oof)
                bi = max(scores, key=lambda c: (np.nan_to_num(scores[c], nan=-1), -c))
                best, inner_score = configs[bi], scores[bi]
                thr = best_mcc_threshold(ytr, oof_by_cfg[bi])
                p_te = fit_predict(best, X[tr], ytr, X[te], seed=bi)
            for i, p in zip(te, p_te):
                rows.append((arm, level, trait, lab.genome_id[i], int(k), y[i], float(p), float(thr),
                             json.dumps(best, sort_keys=True), inner_score))
        log(f"[{arm}/{level}] {trait} done")
    return pd.DataFrame(rows, columns=["arm", "level", "trait", "genome_id", "fold", "y", "prob", "threshold",
                                       "config", "inner_auprc"])


def taxonomic_nn_null(lab, traits, Y, level):
    """Taxonomic group-mean null: mean training label of the genomes sharing the lowest taxonomic rank with the
    test genome (genus -> family -> order -> class -> phylum -> domain -> global prevalence). Stands in for a
    phylogenetic nearest-neighbour null (no tree). Arm id kept as "null_taxNN"."""
    fold = lab[f"fold_{level}"].to_numpy()
    rows = []
    for k in sorted(np.unique(fold)):
        tr, te = fold != k, fold == k
        for j, trait in enumerate(traits):
            ytr = pd.Series(Y[tr, j])
            prior = ytr.mean()
            means = {r: ytr.groupby(lab.loc[tr, r].fillna("NA").to_numpy()).mean() for r in TAX}
            for i in np.where(te)[0]:
                p, used = prior, "global"
                for r in reversed(TAX):  # genus first
                    v = lab.at[i, r]
                    if pd.notna(v) and v in means[r].index:
                        p, used = means[r][v], r
                        break
                rows.append(("null_taxNN", level, trait, lab.genome_id[i], int(k), Y[i, j], float(p), np.nan,
                             json.dumps({"rank": used}), np.nan))
    return pd.DataFrame(rows, columns=["arm", "level", "trait", "genome_id", "fold", "y", "prob", "threshold",
                                       "config", "inner_auprc"])


def prevalence_null(lab, traits, Y, level):
    fold = lab[f"fold_{level}"].to_numpy()
    rows = []
    for k in sorted(np.unique(fold)):
        tr, te = fold != k, fold == k
        for j, trait in enumerate(traits):
            p = Y[tr, j].mean()
            for i in np.where(te)[0]:
                rows.append(("null_prevalence", level, trait, lab.genome_id[i], int(k), Y[i, j], float(p), np.nan,
                             "{}", np.nan))
    return pd.DataFrame(rows, columns=["arm", "level", "trait", "genome_id", "fold", "y", "prob", "threshold",
                                       "config", "inner_auprc"])


def point_metrics(y, p, thr):
    out = {"auprc": safe_ap(y, p),
           "auroc": roc_auc_score(y, p) if 0 < y.sum() < len(y) else np.nan,
           "brier": float(np.mean((p - y) ** 2))}
    if np.all(np.isnan(thr)):
        out["mcc"] = np.nan
        out["f1"] = np.nan
    else:
        yhat = p >= thr
        tp = np.sum(yhat & (y == 1))
        out["mcc"] = matthews_corrcoef(y, yhat)
        out["f1"] = 2 * tp / max(yhat.sum() + y.sum(), 1)
    return out


def fit_models(cfg, X, y, n_jobs=1):
    """Fit the estimator(s) for a tabular config; xgb is a 3-seed ensemble. Shared by CV and final refit."""
    import warnings

    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from xgboost import XGBClassifier
    if cfg["model"] == "logreg":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            return [make_pipeline(StandardScaler(), LogisticRegression(C=cfg["C"], class_weight="balanced",
                                                                       max_iter=5000)).fit(X, y)]
    if cfg["model"] == "xgb":
        return [XGBClassifier(n_estimators=300, learning_rate=0.05, max_depth=cfg["max_depth"], subsample=0.8,
                              colsample_bytree=0.5, min_child_weight=1, scale_pos_weight=(len(y) - y.sum()) / y.sum(),
                              n_jobs=n_jobs, random_state=1000 + s, verbosity=0).fit(X, y) for s in range(3)]
    raise ValueError(cfg)


class FittedModel:
    """Released per-trait model: feature column list + fitted estimator(s) (probabilities averaged)."""

    def __init__(self, columns, cfg, models):
        self.columns, self.cfg, self.models = columns, cfg, models

    def predict_proba(self, Xdf):
        missing = set(self.columns) - set(Xdf.columns)
        assert not missing, f"feature columns missing at prediction time: {sorted(missing)[:5]}"
        X = Xdf[self.columns].to_numpy(dtype=np.float32)
        return np.mean([m.predict_proba(X)[:, 1] for m in self.models], axis=0)


def taxonomic_prior(codes_ref, y_ref, codes_query, exclude_col, same_set=False):
    """Leakage-free taxonomic prior (research loop R2).

    codes_*: integer rank codes, columns = TAX order (domain..genus), -1 = missing.
    For each query genome: mean y of reference genomes sharing the lowest possible rank, never using reference
    genomes in the query's own `exclude_col` group (genus or family: the held-out rank). When query and reference
    are the same set (training rows), this also excludes the genome itself, so a training genome's prior never
    sees its own label. Falls back to the reference prevalence.
    """
    prior = np.full(len(codes_query), y_ref.mean())
    for i, q in enumerate(codes_query):
        allowed = codes_ref[:, exclude_col] != q[exclude_col] if q[exclude_col] >= 0 else np.ones(len(codes_ref), bool)
        if same_set:
            allowed[i] = False
        for r in range(exclude_col - 1, -1, -1):  # ranks above the held-out rank, lowest first
            if q[r] < 0:
                continue
            m = allowed & (codes_ref[:, r] == q[r])
            if m.any():
                prior[i] = y_ref[m].mean()
                break
    return prior


def oof_matrix(oof_dir, arm, level, lab, traits):
    """(n_genomes x n_traits) out-of-fold score matrix for one arm/level, NaN where no prediction."""
    f = Path(oof_dir) / ("nulls.parquet" if arm.startswith("null_") else f"{arm}.parquet")
    d = pd.read_parquet(f)
    d = d[(d.arm == arm) & (d.level == level)]
    gid = {g: i for i, g in enumerate(lab.genome_id)}
    P = np.full((len(lab), len(traits)), np.nan)
    for j, t in enumerate(traits):
        dt = d[d.trait == t]
        P[dt.genome_id.map(gid).to_numpy(), j] = dt.prob
    return P


def paired_macro_delta(lab, traits, Y, Pa, Pb, level, n_boot=1000, seed=20260917):
    """Macro-AUPRC(a) - macro-AUPRC(b) with a paired cluster bootstrap over the level's groups."""
    obs = ~np.isnan(Y)
    groups = lab[f"group_{level}"].to_numpy()
    ug = np.unique(groups)
    members = {g: np.where(groups == g)[0] for g in ug}
    rng = np.random.default_rng(seed)

    def macro(P, idx):
        return np.nanmean([safe_ap(Y[idx, j][obs[idx, j]], P[idx, j][obs[idx, j]]) for j in range(len(traits))])

    allidx = np.arange(len(lab))
    point = macro(Pa, allidx) - macro(Pb, allidx)
    d = []
    for _ in range(n_boot):
        idx = np.concatenate([members[g] for g in rng.choice(ug, len(ug))])
        d.append(macro(Pa, idx) - macro(Pb, idx))
    return point, float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


class FittedFusion:
    """Released fusion model: equal-weight average of component models, each on its own feature columns."""

    def __init__(self, parts, cfg):
        self.parts, self.cfg = parts, cfg  # parts: list of FittedModel

    def predict_proba(self, Xdf):
        return np.mean([p.predict_proba(Xdf) for p in self.parts], axis=0)


def fit_calibrator(y, oof):
    """Platt scaling on inner out-of-fold predictions (no test data): logistic regression on the logit."""
    from sklearn.linear_model import LogisticRegression
    p = np.clip(np.asarray(oof, dtype=float), 1e-6, 1 - 1e-6)
    if len(set(y)) < 2:
        return None
    return LogisticRegression(C=1e6, solver="lbfgs").fit(np.log(p / (1 - p)).reshape(-1, 1), y)


def apply_calibrator(cal, p):
    if cal is None:
        return np.asarray(p, dtype=float)
    p = np.clip(np.asarray(p, dtype=float), 1e-6, 1 - 1e-6)
    return cal.predict_proba(np.log(p / (1 - p)).reshape(-1, 1))[:, 1]
