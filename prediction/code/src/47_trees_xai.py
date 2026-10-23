"""Interpretable tree models and XAI over KO profiles (docs/preregistration_trees_xai.md).

usage:
  python src/47_trees_xai.py cv      --task {da,pw,ox} --model {dt,rf,et,l1,figs,nn} [--jobs N]
                                    (nn runs on the Apple GPU via MPS; use --jobs 1 so one process owns the GPU)
  python src/47_trees_xai.py nulls   --task ox          (taxonomic / prevalence nulls for the new OX target)
  python src/47_trees_xai.py metrics --task {da,pw,ox}
  python src/47_trees_xai.py xai     --task {da,pw,ox} --model {dt,rf,et,l1}
Outputs: results/trees_xai/<task>/oof/<model>.parquet, metrics_*.csv, xai_*.csv
"""
import argparse
import json
import os
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from pathway_common import RES, TRAIT_MODULES, ko_matrix, load_pathways
from redox_common import (LEVELS, PROC, RULE_GROUPS, load_features, load_labelled, nested_cv, point_metrics,
                          prevalence_null, safe_ap, taxonomic_nn_null)

OUT = RES / "trees_xai"
# Pre-registration BT: TX_FEATURES=bacformer_large_mean runs the same learners on another feature block (da task);
# arms get the suffix _BF and are written to results/trees_xai/<task>_bf/oof/ so the KO-panel outputs stay untouched.
FEATURES = os.environ.get("TX_FEATURES")
SUFFIX = "_BF" if FEATURES else ""


def oof_dir(task):
    return OUT / (task + ("_bf" if FEATURES else "")) / "oof"


SEEDS = (0, 100, 200)
TASKS = ("da", "pw", "ox")
CONFIGS = {
    "dt": [{"model": "dt", "max_depth": d} for d in (2, 3, 4, 6)],
    "dt3": [{"model": "dt", "max_depth": d} for d in (2, 3)],   # Addendum L: T2 clause 2 ("a depth <= 3 tree")
    "rf": [{"model": "rf", "msl": m} for m in (1, 5)],
    "et": [{"model": "et", "msl": m} for m in (1, 5)],
    "xgb": [{"model": "xgb", "max_depth": d} for d in (2, 4)],   # Pre-registration BT: the RedoxKO XGBoost grid
    "l1": [{"model": "l1", "C": c} for c in (0.03, 0.1, 0.3, 1.0)],
    "figs": [{"model": "figs", "max_rules": r} for r in (5, 10, 20)],   # Addendum A, A3
    "nn": [{"model": "nn", "wd": w} for w in (1e-4, 1e-2)],              # Addendum C (Apple GPU)
    # Addendum S: a real hyper-parameter search for the MLP (width, dropout, learning rate, epochs, weight decay)
    "nn_wide": [{"model": "nn", "wd": w, "hidden": h, "dropout": d, "lr": lr, "epochs": e}
                for h in (64, 256, 512) for d in (0.1, 0.5) for w in (1e-4, 1e-2) for lr, e in ((1e-3, 40), (3e-4, 120))],
}
OXYGEN_USE = {"aerobic": 1, "obligate aerobic": 1, "facultative": 1, "microaerophilic": 1,
              "anaerobic": 0, "obligate anaerobic": 0}
TREE_JOBS = 2  # threads per forest; trait-level parallelism comes from --jobs


# ----------------------------------------------------------------------------------------------------------- data
def load_task(task):
    """lab (with folds, groups, taxonomy), traits, Y, feature DataFrame (rows aligned to lab)."""
    if task == "da":
        lab, traits, Y = load_labelled()
        X = load_features(FEATURES or "energy_panel", lab.genome_id)
    else:
        lab, traits, Y = load_pathways()
        if task == "ox":
            mad = pd.read_csv(PROC.parent / "raw/madin/condensed_species_NCBI.csv", low_memory=False) \
                    .drop_duplicates("species_tax_id").set_index("species_tax_id")
            ox = lab.species_tax_id.map(mad.metabolism).astype("string").str.lower().str.strip().map(OXYGEN_USE)
            keep = ox.notna().to_numpy()
            lab = lab[keep].reset_index(drop=True)
            traits = ["uses_O2"]
            Y = ox[keep].to_numpy(dtype=float).reshape(-1, 1)
        X = ko_matrix(lab.genome_id)
        prev = X.mean()
        X = X.loc[:, (prev >= 0.01) & (prev <= 0.99)]  # label-free prevalence filter (pre-registered)
    X = X.loc[:, X.std(axis=0) > 0]
    X.index = lab.genome_id.to_numpy()
    return lab, traits, Y, X


# ---------------------------------------------------------------------------------------------------------- models
def torch_device():
    """Apple-silicon GPU (MPS) when available; NN_DEVICE=cpu forces the CPU. See docs/deviations.md #22 on MPS
    reproducibility: dropout masks come from the device's own generator, so an MPS run is a different seed."""
    import os
    import torch
    return "cpu" if os.environ.get("NN_DEVICE") == "cpu" or not torch.backends.mps.is_available() else "mps"


class NNClassifier:
    """Addendum C: standardised KO input -> 128 ReLU (dropout 0.3) -> 1 logit, AdamW, balanced BCE, 40 epochs."""

    def __init__(self, wd, seed, epochs=40, hidden=128, dropout=0.3, lr=1e-3, batch=64):
        self.wd, self.seed, self.epochs, self.hidden, self.dropout, self.lr, self.batch = \
            wd, seed, epochs, hidden, dropout, lr, batch

    def fit(self, X, y):
        import torch
        from torch import nn
        torch.manual_seed(self.seed)
        self.dev = torch_device()
        self.mu, sd = X.mean(0), X.std(0)
        # zero-variance columns (a gene absent from - or present in - every training genome, common under
        # genus-blocked folds) get scale 1, as in sklearn's StandardScaler. Dividing by sd + 1e-6 turned such a
        # gene into a +-1e6 input on any test genome that differed, corrupting predictions and attributions.
        self.sd = np.where(sd > 1e-8, sd, 1.0)
        Xt = torch.tensor((X - self.mu) / self.sd, dtype=torch.float32, device=self.dev)
        yt = torch.tensor(y, dtype=torch.float32, device=self.dev)
        n_pos = float(y.sum())
        pw = torch.tensor(min((len(y) - n_pos) / max(n_pos, 1.0), 50.0), device=self.dev)
        self.net = nn.Sequential(nn.Linear(X.shape[1], self.hidden), nn.ReLU(), nn.Dropout(self.dropout),
                                 nn.Linear(self.hidden, 1)).to(self.dev)
        opt = torch.optim.AdamW(self.net.parameters(), lr=self.lr, weight_decay=self.wd)
        lossf = nn.BCEWithLogitsLoss(pos_weight=pw)
        g = torch.Generator().manual_seed(self.seed)          # batch order drawn on CPU: identical on any device
        self.net.train()
        for _ in range(self.epochs):
            for idx in torch.randperm(len(y), generator=g).to(self.dev).split(self.batch):
                opt.zero_grad()
                lossf(self.net(Xt[idx]).squeeze(1), yt[idx]).backward()
                opt.step()
        self.net.eval()
        return self

    def _x(self, X):
        import torch
        return torch.tensor((X - self.mu) / self.sd, dtype=torch.float32, device=self.dev)

    def predict_proba(self, X):
        import torch
        with torch.no_grad():
            p = torch.sigmoid(self.net(self._x(X)).squeeze(1)).cpu().numpy()
        return np.stack([1 - p, p], axis=1)

    def expected_gradients(self, X, Xref, n_samples=32, chunk=128, seed=0):
        """Expected Gradients (Addendum F): E over baselines x' ~ training genomes and alpha ~ U(0,1) of
        (x - x') * grad f(x' + alpha (x - x')). Reference = the training distribution, as for TreeSHAP."""
        import torch
        g = torch.Generator().manual_seed(seed)
        x = self._x(X)
        ref = self._x(Xref)
        out = torch.zeros_like(x)
        for i in range(0, len(x), chunk):
            xi = x[i:i + chunk]
            acc = torch.zeros_like(xi)
            for _ in range(n_samples):
                bi = torch.randint(len(ref), (len(xi),), generator=g).to(self.dev)
                a = torch.rand(len(xi), 1, generator=g).to(self.dev)
                x0 = ref[bi]
                z = (x0 + a * (xi - x0)).detach().requires_grad_(True)
                self.net(z).sum().backward()
                acc += (xi - x0) * z.grad
            out[i:i + chunk] = acc / n_samples
        return out.detach().cpu().numpy()

    def integrated_gradients(self, X, steps=32, chunk=256):
        """IG of the logit w.r.t. each KO, baseline = the all-KOs-absent genome (raw 0, standardised)."""
        import torch
        x = self._x(X)
        x0 = self._x(np.zeros((1, X.shape[1]), dtype=np.float32))
        out = torch.zeros_like(x)
        alphas = torch.linspace(1.0 / steps, 1.0, steps, device=self.dev)
        for i in range(0, len(x), chunk):
            xi = x[i:i + chunk]
            acc = torch.zeros_like(xi)
            for a in alphas:
                z = (x0 + a * (xi - x0)).detach().requires_grad_(True)
                self.net(z).sum().backward()
                acc += z.grad
            out[i:i + chunk] = (xi - x0) * acc / steps
        return out.detach().cpu().numpy()


def make(cfg, seed):
    from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.tree import DecisionTreeClassifier
    m = cfg["model"]
    if m == "dt":
        return DecisionTreeClassifier(max_depth=cfg["max_depth"], min_samples_leaf=3, class_weight="balanced",
                                      random_state=seed)
    if m in ("rf", "et"):
        cls = RandomForestClassifier if m == "rf" else ExtraTreesClassifier
        return cls(n_estimators=500, max_features="sqrt", min_samples_leaf=cfg["msl"],
                   class_weight="balanced_subsample", n_jobs=TREE_JOBS, random_state=seed)
    if m == "nn":
        return NNClassifier(wd=cfg["wd"], seed=seed, hidden=cfg.get("hidden", 128), dropout=cfg.get("dropout", 0.3),
                            lr=cfg.get("lr", 1e-3), epochs=cfg.get("epochs", 40))
    if m == "figs":
        from imodels import FIGSClassifier
        return FIGSClassifier(max_rules=cfg["max_rules"], random_state=seed)
    if m == "l1":
        return make_pipeline(StandardScaler(), LogisticRegression(penalty="l1", solver="liblinear", C=cfg["C"],
                                                                  class_weight="balanced", max_iter=2000))
    raise ValueError(cfg)


def fit_all(cfg, Xtr, ytr):
    """Fitted estimator(s): forests are averaged over 3 seeds, DT and L1 are deterministic."""
    if cfg["model"] == "xgb":
        from redox_common import fit_models
        return fit_models(cfg, Xtr, ytr, n_jobs=TREE_JOBS)
    seeds = SEEDS if cfg["model"] in ("rf", "et", "nn") else (0,)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return [make(cfg, s).fit(Xtr, ytr) for s in seeds]


def fit_predict(cfg, Xtr, ytr, Xte, seed=0):
    return np.mean([m.predict_proba(Xte)[:, 1] for m in fit_all(cfg, Xtr, ytr)], axis=0)


# -------------------------------------------------------------------------------------------------------------- cv
def cmd_cv(task, model, jobs):
    lab, traits, Y, Xdf = load_task(task)
    X = Xdf.to_numpy(dtype=np.float32)
    print(f"[{task}/{model}] {len(lab)} genomes, {len(traits)} traits, {X.shape[1]} features", flush=True)
    out = Parallel(n_jobs=jobs)(
        delayed(nested_cv)(model.upper() + SUFFIX, X, lab, [traits[j]], Y[:, [j]], CONFIGS[model], fit_predict, lvl)
        for lvl in LEVELS for j in range(len(traits)))
    d = pd.concat(out, ignore_index=True)
    oof_dir(task).mkdir(parents=True, exist_ok=True)
    d.to_parquet(oof_dir(task) / f"{model.upper()}{SUFFIX}.parquet", index=False)
    print(f"[{task}/{model}] written {len(d)} OOF rows", flush=True)


def cmd_nulls(task):
    lab, traits, Y, _ = load_task(task)
    d = pd.concat([f(lab, traits, Y, lvl) for lvl in LEVELS for f in (prevalence_null, taxonomic_nn_null)],
                  ignore_index=True)
    (OUT / task / "oof").mkdir(parents=True, exist_ok=True)
    d.to_parquet(OUT / task / "oof" / "nulls.parquet", index=False)
    print(f"[{task}] nulls written ({len(d)} rows)")


# --------------------------------------------------------------------------------------------------------- metrics
COMPARATORS = {"da": {"B2": RES / "oof/B2.parquet", "B2XL": RES / "oof/B2XL.parquet"},
               "pw": {"PB1": RES / "pathways/oof/PB1.parquet", "PB2": RES / "pathways/oof/PB2.parquet"},
               "ox": {}}
# pre-registered contrasts of the foundation-model addenda (pathway FM and fine-tuning addenda)
FM_CONTRASTS = [("PE", "PB2"), ("PE_BFL", "PB2"), ("PE_MLP", "PB2"), ("PE", "PE_BFL"),
                ("FTP_BFL2", "PB2"), ("FTP_BFL2", "FTP_BFL0"), ("FTP_FUSE", "PB2"),
                # Addenda O, Q, R: AutoML, tabular foundation model, multi-task transfer
                ("AM", "B2"), ("AM", "RF"), ("AM", "ET"), ("AMG", "B2"), ("AMG", "RF"), ("AMG", "ET"), ("AMG", "AM"), ("TPFN", "B2"), ("TPFN", "RF"), ("TPFN", "B2XL"),
                ("MT", "ST"), ("MT", "B2"), ("ST", "B2"),
                # Addenda S, T: wide MLP grid and self-supervised pretraining
                ("SSL", "ST"), ("SSL", "RF"), ("SSL", "B2"), ("NN_WIDE", "NN"), ("NN_WIDE", "RF"), ("NN_WIDE", "B2")]
NULL_FILE = {"da": RES / "oof/nulls.parquet", "pw": RES / "pathways/oof/nulls.parquet"}


def cmd_metrics(task, n_boot=1000, seed=20260918):
    import importlib
    m8 = importlib.import_module("08_metrics")
    lab, traits, Y, _ = load_task(task)
    frames = [pd.read_parquet(f) for f in sorted((OUT / task / "oof").glob("*.parquet"))]
    for name, f in COMPARATORS[task].items():
        if f.exists():
            frames.append(pd.read_parquet(f).assign(arm=name))
    if task in NULL_FILE and not (OUT / task / "oof/nulls.parquet").exists():
        frames.append(pd.read_parquet(NULL_FILE[task]))
    oof = pd.concat(frames, ignore_index=True)
    oof = oof[oof.trait.isin(traits) & oof.genome_id.isin(set(lab.genome_id))]
    gid = {g: i for i, g in enumerate(lab.genome_id)}
    rng = np.random.default_rng(seed)
    per_trait, macro, paired = [], [], []
    for level in LEVELS:
        groups = lab[f"group_{level}"].to_numpy()
        ug = np.unique(groups)
        mem = {g: np.where(groups == g)[0] for g in ug}
        boots = [np.concatenate([mem[g] for g in rng.choice(ug, len(ug), replace=True)]) for _ in range(n_boot)]
        P, T = {}, {}
        for arm in oof.arm.unique():
            d = oof[(oof.arm == arm) & (oof.level == level)]
            pr, th = np.full(Y.shape, np.nan), np.full(Y.shape, np.nan)
            for j, t in enumerate(traits):
                dt = d[d.trait == t]
                idx = dt.genome_id.map(gid).to_numpy()
                pr[idx, j], th[idx, j] = dt.prob, dt.threshold
            if np.isnan(pr).any():
                print(f"  skip {arm}/{level}: incomplete ({int(np.isnan(pr).sum())} missing)")
                continue
            P[arm], T[arm] = pr, th
        B = {a: np.stack([[safe_ap(Y[b, j], P[a][b, j]) for b in boots] for j in range(len(traits))], 1)
             for a in P}
        for arm in P:
            rows = []
            for j, t in enumerate(traits):
                y, p, th = Y[:, j], P[arm][:, j], T[arm][:, j]
                pm = point_metrics(y, p, th)
                has_thr = not np.all(np.isnan(th))
                from sklearn.metrics import matthews_corrcoef
                bm = [matthews_corrcoef(Y[b, j], P[arm][b, j] >= T[arm][b, j]) for b in boots] if has_thr else []
                dnull = B[arm][:, j] - B["null_taxNN"][:, j]
                r = dict(task=task, arm=arm, split_level=level, trait=t, n_pos=int(y.sum()), prevalence=y.mean(),
                         auprc=pm["auprc"], auroc=pm["auroc"], mcc=pm["mcc"], f1=pm["f1"],
                         threshold_median=float(np.nanmedian(th)) if has_thr else np.nan)
                r["auprc_ci_low"], r["auprc_ci_high"] = m8.ci(B[arm][:, j])
                r["mcc_ci_low"], r["mcc_ci_high"] = m8.ci(bm) if has_thr else (np.nan, np.nan)
                r["null_taxNN_auprc"] = safe_ap(y, P["null_taxNN"][:, j])
                r["delta_vs_null"] = pm["auprc"] - r["null_taxNN_auprc"]
                r["delta_vs_null_ci_low"], r["delta_vs_null_ci_high"] = m8.ci(dnull)
                r["p_vs_null"] = m8.boot_p(dnull) if arm != "null_taxNN" else np.nan
                rows.append(r)
            for r, qq in zip(rows, m8.bh([r["p_vs_null"] for r in rows])):
                r["q_vs_null"] = qq
            per_trait += rows
            lo, hi = m8.ci(np.nanmean(B[arm], axis=1))
            macro.append(dict(task=task, arm=arm, split_level=level, n_traits=len(traits),
                              macro_auprc=np.mean([r["auprc"] for r in rows]), macro_auprc_ci_low=lo,
                              macro_auprc_ci_high=hi, macro_mcc=np.nanmean([r["mcc"] for r in rows]),
                              macro_delta_vs_null=np.mean([r["delta_vs_null"] for r in rows])))
        ref = [a for a in COMPARATORS[task] if a in P]
        pairs = [(a, b) for a in ("DT", "RF", "ET", "L1", "FIGS", "NN") if a in P
                 for b in ref + [x for x in ("RF", "DT") if x in P and x != a]]
        pairs += [p for p in FM_CONTRASTS if p[0] in P and p[1] in P]   # foundation-model addenda
        for a, b in pairs:
            dm = np.nanmean(B[a], axis=1) - np.nanmean(B[b], axis=1)
            lo, hi = m8.ci(dm)
            paired.append(dict(task=task, split_level=level, contrast=f"{a} - {b}",
                               delta_macro_auprc=float(np.nanmean(dm)), ci_low=lo, ci_high=hi,
                               tied=bool(lo <= 0 <= hi)))
    pt = pd.DataFrame(per_trait)
    # the pre-registered "predictable" gate, at family level
    fam = pt[pt.split_level == "family"].set_index(["arm", "trait"])
    pt["predictable"] = [bool(fam.loc[(a, t), "mcc_ci_low"] > 0.2 and fam.loc[(a, t), "q_vs_null"] < 0.05
                              and fam.loc[(a, t), "delta_vs_null"] > 0) if (a, t) in fam.index else False
                         for a, t in zip(pt.arm, pt.trait)]
    pt.to_csv(OUT / task / "metrics_per_trait.csv", index=False, float_format="%.4f")
    mc = pd.DataFrame(macro)
    npred = pt[pt.split_level == "family"].groupby("arm").predictable.sum()
    mc["n_predictable"] = mc.arm.map(npred)
    mc.to_csv(OUT / task / "metrics_macro.csv", index=False, float_format="%.4f")
    pd.DataFrame(paired).to_csv(OUT / task / "paired.csv", index=False, float_format="%.4f")
    print(mc.sort_values(["split_level", "macro_auprc"], ascending=[True, False]).round(3).to_string(index=False))
    if paired:
        print(pd.DataFrame(paired).round(3).to_string(index=False))


# ------------------------------------------------------------------------------------------------------------- xai
def known_markers(task, trait, reference_groups):
    """KOs that literature/curation already ties to this trait (the validation set for XAI)."""
    known = set()
    for mods in (TRAIT_MODULES.get(trait, []),):
        if mods:
            kmods = pd.read_csv(PROC / "kegg_modules_energy.tsv", sep="\t").set_index("module")
            for m in mods:
                if m in kmods.index:
                    known |= set(__import__("re").findall(r"K\d{5}", kmods.at[m, "definition"]))
    for kos in RULE_GROUPS.get(trait, {}).values():
        known |= set(kos)
    for g, kos in reference_groups.items():
        if REFERENCE_TRAIT.get(g) in (trait,) or (trait in REFERENCE_PATHWAY.get(g, ())):
            known |= set(kos)
    if trait == "uses_O2":
        known |= {"K02274", "K02275", "K02276", "K00404", "K00405", "K00406", "K00425", "K00426",
                  "K15408", "K02297", "K02298", "K02299", "K02300"}  # terminal oxidases (aa3, cbb3, bd, bo3)
    return known


REFERENCE_TRAIT = {"Sulfur oxidizing (SOB)": "donor_sulfur_oxidation", "Ammonia oxidizing (AOB/AOA)":
                   "donor_ammonia_oxidation", "Nitrite oxidizing (NOB)": "donor_nitrite_oxidation",
                   "Hydrogen oxidizing": "donor_hydrogen_oxidation", "CO oxidizing (carboxydotroph)":
                   "donor_CO_oxidation", "Aerobic respiration": "acceptor_oxygen",
                   "Nitrate-reducing (denitrifiers)": "acceptor_nitrate", "Nitrate-reducing (DNRA)": "acceptor_nitrate",
                   "Sulfate-reducing (SRB)": "acceptor_sulfate", "Sulfur-reducing": "acceptor_sulfur",
                   "Methanogens (hydrogenotrophic)": "acceptor_CO2_methanogen_acetogen",
                   "Acetogens (homoacetogens)": "acceptor_CO2_methanogen_acetogen"}
REFERENCE_PATHWAY = {"Sulfur oxidizing (SOB)": ("sulfur_oxidation",), "Ammonia oxidizing (AOB/AOA)": ("nitrification",),
                     "Hydrogen oxidizing": ("hydrogen_oxidation",), "Nitrate-reducing (denitrifiers)":
                     ("denitrification", "nitrate_reduction"), "Nitrate-reducing (DNRA)": ("dnra", "nitrate_reduction"),
                     "Sulfate-reducing (SRB)": ("sulfate_reduction",), "Sulfur-reducing": ("sulfur_reduction",),
                     "Methanogens (hydrogenotrophic)": ("methanogenesis",), "Acetogens (homoacetogens)":
                     ("acetogenesis",), "Methane oxidizing (methanotroph)": ("methanotrophy", "methylotrophy"),
                     "Aerobic respiration": ("uses_O2", "aerobic_chemoheterotrophy")}


def reference_groups():
    import importlib.util
    spec = importlib.util.spec_from_file_location("jc", RES.parent / "src/40_jaccard_correlations.py")
    jc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(jc)
    return jc.load_html(RES.parent / "jaccard_similarity.html")[1]


def permutation_on(est, Xte, yte, candidates, seed, repeats=10):
    """Drop in held-out AUPRC when one KO column is shuffled, for the candidate columns only (Addendum E)."""
    rng = np.random.default_rng(seed)
    base = safe_ap(yte, est.predict_proba(Xte)[:, 1])
    out = np.zeros(Xte.shape[1])
    for c in candidates:
        Xp = Xte.copy()
        drops = []
        for _ in range(repeats):
            Xp[:, c] = rng.permutation(Xte[:, c])
            drops.append(base - safe_ap(yte, est.predict_proba(Xp)[:, 1]))
        out[c] = float(np.mean(drops))
    return out


def explain_fold(model, cfg, Xtr, ytr, Xte, yte, cols, seed, known=()):
    """Held-out importance for one outer fold (two-stage for RF/ET, Addendum E); rules (DT); IG (NN)."""
    ms = fit_all(cfg, Xtr, ytr)
    est = ms[0]
    res = {}
    if model in ("rf", "et"):
        res["impurity_importance"] = np.mean([m.feature_importances_ for m in ms], axis=0)
        screen = res["impurity_importance"]
        if model == "rf":
            import shap
            rng = np.random.default_rng(seed)
            pos, neg = np.where(yte == 1)[0], np.where(yte == 0)[0]
            sub = np.concatenate([pos, rng.choice(neg, min(200, len(neg)), replace=False)])
            sv = shap.TreeExplainer(est).shap_values(Xte[sub], check_additivity=False)
            sv = sv[..., 1] if np.ndim(sv) == 3 else (sv[1] if isinstance(sv, list) else sv)
            res["mean_abs_shap"] = np.abs(sv).mean(axis=0)
            screen = res["mean_abs_shap"]
            present = Xte[sub] > 0
            with np.errstate(invalid="ignore"):
                res["shap_direction"] = np.array([sv[present[:, k], k].mean() - sv[~present[:, k], k].mean()
                                                  if present[:, k].any() and (~present[:, k]).any() else 0.0
                                                  for k in range(len(cols))])
        cand = set(np.argsort(-screen)[:50]) | {i for i, c in enumerate(cols) if c in set(known)}
        res["perm_importance"] = (permutation_on(est, Xte, yte, sorted(cand), seed) if yte.sum() > 0
                                  else np.zeros(len(cols)))
        res["permutation_tested"] = np.isin(np.arange(len(cols)), sorted(cand)).astype(float)
    elif model == "dt":
        res["impurity_importance"] = est.feature_importances_
        res["rules"] = tree_rules(est, cols, Xte, yte)
    elif model == "figs":
        res["impurity_importance"] = np.asarray(est.feature_importances_, dtype=float)
    elif model == "nn":
        # NN_ATTRIB=ig_zero -> pre-registered Integrated Gradients, all-absent baseline (Addendum C);
        # default -> Expected Gradients against the training fold (Addendum F, post hoc). Averaged over the 3 seeds.
        import os
        if os.environ.get("NN_ATTRIB") == "ig_zero":
            ig = np.mean([m.integrated_gradients(Xte) for m in ms], axis=0)
        else:
            ig = np.mean([m.expected_gradients(Xte, Xtr, seed=seed + i) for i, m in enumerate(ms)], axis=0)
        res["mean_abs_ig"] = np.abs(ig).mean(axis=0)
        res["impurity_importance"] = res["mean_abs_ig"]
        present = Xte > 0
        res["ig_direction"] = np.array([ig[present[:, k], k].mean() if present[:, k].any() else 0.0
                                        for k in range(len(cols))])
    elif model == "l1":
        coef = est[-1].coef_.ravel()
        res["coef"] = coef
        res["impurity_importance"] = np.abs(coef)
    return res


def tree_rules(est, cols, Xte, yte):
    """Root-to-leaf paths that predict the positive class, scored on the held-out fold."""
    t = est.tree_
    rules = []

    def walk(node, conds):
        if t.children_left[node] == -1:
            v = t.value[node][0]
            if v[1] > v[0]:
                m = np.ones(len(Xte), bool)
                for f_, present in conds:
                    m &= (Xte[:, f_] > 0.5) == present
                rules.append(dict(rule=" AND ".join(f"{cols[f_]} {'present' if p else 'absent'}" for f_, p in conds)
                                  or "(all)", depth=len(conds), heldout_support=int(m.sum()),
                                  heldout_precision=float(yte[m].mean()) if m.any() else np.nan,
                                  heldout_recall=float((m & (yte == 1)).sum() / max(yte.sum(), 1))))
            return
        f_ = t.feature[node]
        walk(t.children_left[node], conds + [(f_, False)])   # x <= 0.5 -> absent
        walk(t.children_right[node], conds + [(f_, True)])
    walk(0, [])
    return rules


def cmd_xai(task, model, top_k=10):
    from scipy.stats import hypergeom
    import importlib
    m8 = importlib.import_module("08_metrics")
    lab, traits, Y, Xdf = load_task(task)
    X = Xdf.to_numpy(dtype=np.float32)
    cols = list(Xdf.columns)
    oof = pd.read_parquet(OUT / task / "oof" / f"{model.upper()}.parquet")
    oof = oof[oof.level == "genus"]
    fold = lab["fold_genus"].to_numpy()
    cg = reference_groups()
    imp_rows, rule_rows, enr_rows = [], [], []

    def one(j, t):
        dt = oof[oof.trait == t]
        rows_i, rows_r = [], []
        for k in sorted(dt.fold.unique()):
            cfg = json.loads(dt[dt.fold == k].config.iloc[0])
            if not isinstance(cfg, dict) or "model" not in cfg:
                continue
            tr, te = np.where(fold != k)[0], np.where(fold == k)[0]
            r = explain_fold(model, cfg, X[tr], Y[tr, j], X[te], Y[te, j], cols, seed=int(k),
                             known=known_markers(task, t, cg))
            score = r.get("perm_importance", r.get("impurity_importance"))
            order = np.argsort(-score)
            for rank, c in enumerate(order[:50]):
                row = dict(task=task, model=model.upper(), trait=t, fold=int(k), KO=cols[c], rank=rank + 1,
                           score=float(score[c]), prevalence_in_train=float(X[tr, c].mean()))
                for key in ("impurity_importance", "mean_abs_shap", "shap_direction", "coef", "mean_abs_ig",
                            "ig_direction", "permutation_tested"):
                    if key in r:
                        row[key] = float(r[key][c])
                rows_i.append(row)
            for rr in r.get("rules", []):
                rows_r.append(dict(task=task, trait=t, fold=int(k), **rr))
        return rows_i, rows_r

    out = Parallel(n_jobs=1 if model == "nn" else 4)(delayed(one)(j, t) for j, t in enumerate(traits))
    for ri, rr in out:
        imp_rows += ri
        rule_rows += rr
    import os
    sfx = "_zero_baseline" if model == "nn" and os.environ.get("NN_ATTRIB") == "ig_zero" else ""
    imp = pd.DataFrame(imp_rows)
    n_folds = imp.groupby("trait").fold.nunique()
    top = imp[imp["rank"] <= top_k]
    stab = (top.groupby(["trait", "KO"]).agg(n_folds_top10=("fold", "nunique"), mean_rank=("rank", "mean"),
                                             mean_score=("score", "mean"))
            .reset_index())
    stab["stability"] = stab.n_folds_top10 / stab.trait.map(n_folds)
    for key in ("shap_direction", "coef", "mean_abs_shap", "ig_direction", "mean_abs_ig"):
        if key in top.columns:
            stab = stab.merge(top.groupby(["trait", "KO"])[key].mean().rename(key).reset_index(), on=["trait", "KO"])
    stable = stab[stab.stability >= 0.6].copy()
    background = set(cols)
    for t in traits:
        known = known_markers(task, t, cg) & background
        s = stable[stable.trait == t]
        hits = set(s.KO) & known
        stable.loc[stable.trait == t, "known_marker"] = stable.loc[stable.trait == t, "KO"].isin(known)
        N, K, n, x = len(background), len(known), len(s), len(hits)
        p = hypergeom.sf(x - 1, N, K, n) if n and K else np.nan
        enr_rows.append(dict(task=task, model=model.upper(), trait=t, n_known_markers_in_features=K,
                             n_stable_top=n, n_stable_top_known=x,
                             expected_by_chance=n * K / N if N else np.nan, p_hypergeom=p,
                             stable_top_kos=", ".join(s.sort_values("mean_rank").KO)))
    enr = pd.DataFrame(enr_rows)
    enr["q_bh"] = m8.bh(enr.p_hypergeom.to_numpy())
    (OUT / task).mkdir(parents=True, exist_ok=True)
    imp.to_csv(OUT / task / f"xai_importance_{model}{sfx}.csv", index=False, float_format="%.5f")
    stable.sort_values(["trait", "mean_rank"]).to_csv(OUT / task / f"xai_stable_{model}{sfx}.csv", index=False,
                                                      float_format="%.4f")
    enr.to_csv(OUT / task / f"xai_enrichment_{model}{sfx}.csv", index=False, float_format="%.4g")
    if rule_rows:
        pd.DataFrame(rule_rows).to_csv(OUT / task / "xai_rules_dt.csv", index=False, float_format="%.4f")
    print(enr[["trait", "n_stable_top", "n_stable_top_known", "expected_by_chance", "q_bh", "stable_top_kos"]]
          .round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["cv", "nulls", "metrics", "xai"])
    ap.add_argument("--task", choices=TASKS, required=True)
    ap.add_argument("--model", choices=list(CONFIGS))
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    (OUT / a.task).mkdir(parents=True, exist_ok=True)
    if a.cmd == "cv":
        cmd_cv(a.task, a.model, a.jobs)
    elif a.cmd == "nulls":
        cmd_nulls(a.task)
    elif a.cmd == "metrics":
        cmd_metrics(a.task)
    else:
        cmd_xai(a.task, a.model)


if __name__ == "__main__":
    main()
