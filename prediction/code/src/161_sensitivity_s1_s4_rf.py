"""Sensitivity analyses S1 and S4 (docs/deviations.md; src/run_sensitivity.sh, src/run_sensitivity_metrics.sh) for the
final model RedoxKO (random forest on the 852-KO panel), as they were run for B2:
  S1  label_mode = "block_nan": a genome's all-negative donor or acceptor block is treated as unobserved (NaN); the
      forest is re-run in nested CV on the frozen folds with those labels (per category, unobserved genomes dropped,
      as src/06_run_arms.py --label-mode block_nan does), configurations = 47_trees_xai CONFIGS["rf"].
  S4  evaluation restricted to genomes with quality_ok in data/processed/genome_qc.tsv (CheckM2 completeness >= 90%,
      contamination <= 5%); no refit: the forest's primary out-of-fold scores (results/trees_xai/da/oof/RF.parquet) are
      re-scored on that subset, as src/08_metrics.py --subset quality_ok did for B2.
Metrics: macro AUPRC with 95% cluster-bootstrap intervals exactly as src/08_metrics.py (1000 resamples of genus / family
groups, seed 20260917, shared by all arms), recomputed for B2 and the nulls from their existing OOF files (and checked
against results/sensitivity/frozen__block_nan and frozen__main__quality_ok), plus paired RF - B2 and RF - null_taxNN.
usage: python src/161_sensitivity_s1_s4_rf.py [s1|s4|both] [--jobs N] [--smoke]
Output: results/sensitivity_rf/frozen__block_nan/{oof/RF.parquet, metrics_macro.csv, metrics_per_trait.csv,
        paired_comparisons.csv}; results/sensitivity_rf/frozen__main__quality_ok/{metrics_macro.csv, ...}
"""
import argparse
import importlib
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import LEVELS, PROC, RES, load_features, load_labelled, nested_cv, safe_ap  # noqa: E402

tx = importlib.import_module("47_trees_xai")
m8 = importlib.import_module("08_metrics")
tx.TREE_JOBS = 1
CFG = tx.CONFIGS["rf"]


def fp(cfg, Xtr, ytr, Xte, seed=0):
    return tx.fit_predict(cfg, Xtr, ytr, Xte)


def run_s1(jobs, out_dir, traits_limit=None):
    lab, traits, Y = load_labelled("frozen", "block_nan")
    X = load_features("energy_panel", lab.genome_id)
    X = X.loc[:, X.std(axis=0) > 0].to_numpy(dtype=np.float32)
    js = range(len(traits)) if traits_limit is None else range(traits_limit)

    def one(level, j):
        m = ~np.isnan(Y[:, j])
        return nested_cv("RF", X[m], lab.loc[m].reset_index(drop=True), [traits[j]], Y[m][:, [j]], CFG, fp, level,
                         log=lambda s: None)
    oof = pd.concat(Parallel(n_jobs=jobs)(delayed(one)(lv, j) for lv in LEVELS for j in js), ignore_index=True)
    (out_dir / "oof").mkdir(parents=True, exist_ok=True)
    oof.to_parquet(out_dir / "oof" / "RF.parquet", index=False)
    print("S1 forest OOF written", len(oof), flush=True)


def metrics(label_mode, subset, oof_files, out_dir, check_dir):
    """src/08_metrics.py macro block for the given arms (same bootstrap draws)."""
    lab, traits, Y = load_labelled("frozen", label_mode)
    if subset:
        qc = pd.read_csv(PROC / "genome_qc.tsv", sep="\t").set_index("genome_id")[subset].astype(bool)
        keep = lab.genome_id.map(qc).fillna(False).to_numpy(dtype=bool)
        Y = Y.copy()
        Y[~keep] = np.nan
        print(f"subset {subset}: evaluating {keep.sum()} / {len(keep)} genomes")
    obs = ~np.isnan(Y)
    oof = pd.concat([pd.read_parquet(f) for f in oof_files], ignore_index=True)
    arms = ["RF", "B2", "null_prevalence", "null_taxNN"]
    gid = {g: i for i, g in enumerate(lab.genome_id)}
    rng = np.random.default_rng(m8.SEED)
    macro, per_trait, paired = [], [], []
    for level in LEVELS:
        groups = lab[f"group_{level}"].to_numpy()
        ug = np.unique(groups)
        members = {g: np.where(groups == g)[0] for g in ug}
        boots = [np.concatenate([members[g] for g in rng.choice(ug, len(ug), replace=True)]) for _ in range(m8.N_BOOT)]
        P = {}
        for arm in arms:
            d = oof[(oof.arm == arm) & (oof.level == level)]
            probs, ychk = np.full(Y.shape, np.nan), np.full(Y.shape, np.nan)
            for j, t in enumerate(traits):
                dt = d[d.trait == t]
                idx = dt.genome_id.map(gid).to_numpy()
                probs[idx, j], ychk[idx, j] = dt.prob, dt.y
            assert not np.isnan(probs[obs]).any(), f"{arm}/{level}: unscored labels"
            assert np.array_equal(ychk[obs], Y[obs]), f"{arm}/{level}: OOF labels differ from label matrix"
            P[arm] = probs
        B = {a: np.stack([[safe_ap(Y[b, j][obs[b, j]], P[a][b, j][obs[b, j]]) for b in boots]
                          for j in range(len(traits))], axis=1) for a in P}
        pt = {a: np.array([safe_ap(Y[obs[:, j], j], P[a][obs[:, j], j]) for j in range(len(traits))]) for a in P}
        for a in P:
            lo, hi = m8.ci(np.nanmean(B[a], axis=1))
            macro.append(dict(arm=a, split_level=level, n_traits=len(traits), macro_auprc=float(np.mean(pt[a])),
                              macro_auprc_ci_low=lo, macro_auprc_ci_high=hi,
                              macro_delta_vs_null=float(np.mean(pt[a] - pt["null_taxNN"]))))
            for j, t in enumerate(traits):
                lo, hi = m8.ci(B[a][:, j])
                per_trait.append(dict(arm=a, split_level=level, trait=t, n_pos=int(np.nansum(Y[:, j])),
                                      n_observed=int(obs[:, j].sum()), auprc=pt[a][j], auprc_ci_low=lo, auprc_ci_high=hi))
        for a, b in (("RF", "B2"), ("RF", "null_taxNN"), ("B2", "null_taxNN")):
            dm = np.nanmean(B[a], axis=1) - np.nanmean(B[b], axis=1)
            lo, hi = m8.ci(dm)
            paired.append(dict(split_level=level, contrast=f"{a} - {b}", delta_macro_auprc=float(np.mean(pt[a] - pt[b])),
                               ci_low=lo, ci_high=hi, p_boot=m8.boot_p(dm), tied=bool(lo <= 0 <= hi)))
    out_dir.mkdir(parents=True, exist_ok=True)
    mac = pd.DataFrame(macro)
    mac.to_csv(out_dir / "metrics_macro.csv", index=False, float_format="%.4f")
    pd.DataFrame(per_trait).to_csv(out_dir / "metrics_per_trait.csv", index=False, float_format="%.4f")
    pr = pd.DataFrame(paired)
    pr.to_csv(out_dir / "paired_comparisons.csv", index=False, float_format="%.4f")
    if check_dir is not None:
        old = pd.read_csv(check_dir / "metrics_macro.csv").set_index(["arm", "split_level"])
        for _, r in mac[mac.arm != "RF"].iterrows():
            o = old.loc[(r.arm, r.split_level)]
            assert abs(o.macro_auprc - r.macro_auprc) < 6e-5 and abs(o.macro_auprc_ci_low - r.macro_auprc_ci_low) < 6e-5, \
                (r.arm, r.split_level, o.macro_auprc, r.macro_auprc)
        print(f"B2 and null values reproduce {check_dir.relative_to(ROOT)}/metrics_macro.csv")
    print(mac.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(pr.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", nargs="?", default="both", choices=["s1", "s4", "both"])
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--smoke", action="store_true", help="tiny forests on 2 categories; output to $SMOKE_DIR")
    args = ap.parse_args()
    base = RES / "sensitivity_rf"
    if args.smoke:
        _make = tx.make
        tx.make = lambda cfg, seed: _make(cfg, seed).set_params(n_estimators=10)
        sd = Path(os.environ.get("SMOKE_DIR", "/tmp/rf_smoke")) / "sensitivity_rf/frozen__block_nan"
        run_s1(args.jobs, sd, traits_limit=2)
        return
    if args.what in ("s1", "both"):
        s1 = base / "frozen__block_nan"
        if not (s1 / "oof/RF.parquet").exists():
            run_s1(args.jobs, s1)
        old = RES / "sensitivity/frozen__block_nan"
        metrics("block_nan", None, [s1 / "oof/RF.parquet", old / "oof/B2.parquet", old / "oof/nulls.parquet"], s1, old)
    if args.what in ("s4", "both"):
        metrics("main", "quality_ok", [RES / "trees_xai/da/oof/RF.parquet", RES / "oof/B2.parquet", RES / "oof/nulls.parquet"],
                base / "frozen__main__quality_ok", RES / "sensitivity/frozen__main__quality_ok")


if __name__ == "__main__":
    main()
