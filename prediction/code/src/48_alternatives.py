"""Alternative routes to gene-trait relationships (docs/preregistration_trees_xai.md, Addendum A).

usage:
  python src/48_alternatives.py cmh       --task {da,pw,ox}   A1: per-KO Cochran-Mantel-Haenszel, phylum strata
  python src/48_alternatives.py boruta    --task {da,pw,ox}   A2: Boruta all-relevant selection on training folds
  python src/48_alternatives.py consensus --task {da,pw,ox}   robust = RF-stable + CMH-significant (+ Boruta)
Outputs: results/trees_xai/<task>/alt_cmh.csv, alt_boruta.csv, consensus.csv
"""
import argparse
import importlib
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import chi2

tx = importlib.import_module("47_trees_xai")
m8 = importlib.import_module("08_metrics")
OUT = tx.OUT


def genus_level(lab, X, Y):
    """One data point per split genus: KO present / trait present if >= half its genomes have it."""
    g = lab.group_genus.to_numpy()
    Xg = pd.DataFrame(X, index=g).groupby(level=0).mean().ge(0.5).astype(np.int8)
    Yg = pd.DataFrame(Y, index=g).groupby(level=0).mean().ge(0.5).astype(np.int8)
    ph = pd.Series(lab.phylum.fillna("?").to_numpy(), index=g).groupby(level=0) \
           .agg(lambda v: v.value_counts().index[0]).reindex(Xg.index)
    return Xg, Yg, ph


def cmh_all(Xg, y, strata):
    """Vectorised CMH test of every KO column against one binary trait, stratified by phylum.

    Per stratum s and KO k (a = KO+ trait+, b = KO+ trait-, c = KO- trait+, d = KO- trait-, n = genera):
      E[a] = (a+b)(a+c)/n,  Var[a] = (a+b)(c+d)(a+c)(b+d) / (n^2 (n-1))
      CMH = (|sum(a - E[a])| - 0.5)^2 / sum Var[a]   (continuity-corrected, 1 df)
      Mantel-Haenszel common odds ratio = sum(a d / n) / sum(b c / n)
    Strata with n < 2, or where the trait or the KO does not vary, contribute nothing.
    """
    K = Xg.shape[1]
    s_dev, s_var, s_num, s_den = np.zeros(K), np.zeros(K), np.zeros(K), np.zeros(K)
    for s in np.unique(strata):
        m = strata == s
        n = m.sum()
        if n < 2:
            continue
        Xs, ys = Xg[m], y[m]
        t = ys.sum()
        if t == 0 or t == n:
            continue
        r1 = Xs.sum(axis=0)                        # KO+ genera
        a = (Xs * ys[:, None]).sum(axis=0)         # KO+ trait+
        b, c = r1 - a, t - a
        d = n - r1 - c
        E = r1 * t / n
        V = r1 * (n - r1) * t * (n - t) / (n ** 2 * (n - 1))
        s_dev += a - E
        s_var += V
        s_num += a * d / n
        s_den += b * c / n
    with np.errstate(divide="ignore", invalid="ignore"):
        stat = np.where(s_var > 0, (np.maximum(np.abs(s_dev) - 0.5, 0)) ** 2 / s_var, np.nan)
        orr = np.where(s_den > 0, s_num / s_den, np.inf)
    p = chi2.sf(stat, 1)
    return stat, p, orr, np.sign(s_dev)


def cmd_cmh(task, adjust_size=False):
    """A1 (phylum strata) or A1b (phylum x genome gene-count quintile strata; Addendum B)."""
    lab, traits, Y, Xdf = tx.load_task(task)
    Xg, Yg, ph = genus_level(lab, Xdf.to_numpy(), Y)
    strata = ph.to_numpy()
    if adjust_size:
        # genome breadth confounds presence/absence: aerobes carry ~40 % more KOs even within a phylum, so almost
        # any KO is "enriched" in them. Strata = phylum x quintile of the genus's mean KO count.
        n_ko = pd.Series(Xdf.sum(axis=1).to_numpy(), index=lab.group_genus.to_numpy()).groupby(level=0).mean() \
                 .reindex(Xg.index)
        q5 = pd.qcut(n_ko, 5, labels=False, duplicates="drop").astype(int).to_numpy()
        strata = np.array([f"{a}|q{b}" for a, b in zip(strata, q5)])
    rows = []
    for j, t in enumerate(traits):
        stat, p, orr, sign = cmh_all(Xg.to_numpy(), Yg[j].to_numpy(), strata)
        rows.append(pd.DataFrame({"task": task, "trait": t, "KO": Xdf.columns, "cmh_stat": stat, "p": p,
                                  "mh_odds_ratio": orr, "direction": sign,
                                  "n_genera_with_KO": Xg.sum(axis=0).to_numpy()}))
    d = pd.concat(rows, ignore_index=True)
    d["q_bh"] = np.nan
    for t in traits:  # BH within trait, over all KOs tested
        m = d.trait == t
        d.loc[m, "q_bh"] = m8.bh(d.loc[m, "p"].to_numpy())
    d = d.sort_values(["trait", "p"])
    name = "alt_cmh_size_adjusted.csv" if adjust_size else "alt_cmh.csv"
    d.to_csv(OUT / task / name, index=False, float_format="%.4g")
    sig = d[d.q_bh < 0.05]
    print(f"[{task}] CMH{' (size-adjusted)' if adjust_size else ''}: {len(Xg)} genera, {len(np.unique(strata))} "
          f"strata; significant KO-trait pairs "
          f"(q<0.05): {len(sig)}  (+ {int((sig.direction > 0).sum())} / - {int((sig.direction < 0).sum())})")
    print(sig.groupby("trait").head(3)[["trait", "KO", "mh_odds_ratio", "direction", "q_bh"]]
          .round(4).to_string(index=False))


def boruta_fold(Xtr, ytr, seed):
    from boruta import BorutaPy
    from sklearn.ensemble import RandomForestClassifier
    rf = RandomForestClassifier(n_jobs=2, class_weight="balanced_subsample", max_depth=7, random_state=seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = BorutaPy(rf, n_estimators="auto", max_iter=100, alpha=0.05, random_state=seed, verbose=0)
        b.fit(Xtr, ytr)
    return b.support_, b.support_weak_


def cmd_boruta(task, jobs=3):
    lab, traits, Y, Xdf = tx.load_task(task)
    X = Xdf.to_numpy(dtype=np.float32)
    fold = lab["fold_genus"].to_numpy()
    tasks = [(j, k) for j in range(len(traits)) for k in sorted(np.unique(fold))]

    def run(j, k):
        tr = fold != k
        if Y[tr, j].sum() < 5:
            return j, k, np.zeros(X.shape[1], bool), np.zeros(X.shape[1], bool)
        conf, weak = boruta_fold(X[tr], Y[tr, j].astype(int), seed=int(k))
        return j, k, conf, weak

    res = Parallel(n_jobs=jobs)(delayed(run)(j, k) for j, k in tasks)
    rows = []
    for j, t in enumerate(traits):
        confs = np.array([c for jj, k, c, w in res if jj == j])
        weaks = np.array([w for jj, k, c, w in res if jj == j])
        n_conf, n_weak = confs.sum(axis=0), weaks.sum(axis=0)
        for c in np.where((n_conf + n_weak) > 0)[0]:
            rows.append(dict(task=task, trait=t, KO=Xdf.columns[c], n_folds_confirmed=int(n_conf[c]),
                             n_folds_tentative=int(n_weak[c]), n_folds=len(confs),
                             boruta_confirmed=bool(n_conf[c] >= 3)))
    d = pd.DataFrame(rows).sort_values(["trait", "n_folds_confirmed"], ascending=[True, False])
    d.to_csv(OUT / task / "alt_boruta.csv", index=False)
    print(f"[{task}] Boruta: KOs confirmed in >= 3 of 5 training folds per trait:")
    print(d[d.boruta_confirmed].groupby("trait").KO.apply(lambda s: ", ".join(s.head(12))).to_string())


def cmd_consensus(task):
    stable = pd.read_csv(OUT / task / "xai_stable_rf.csv")
    # primary CMH criterion = size-adjusted (Addendum B); phylum-only A1 kept as a sensitivity column
    cmh = pd.read_csv(OUT / task / "alt_cmh_size_adjusted.csv")
    cmh_phy = pd.read_csv(OUT / task / "alt_cmh.csv")[["trait", "KO", "q_bh"]].rename(columns={"q_bh": "q_bh_phylum_only"})
    bor_f = OUT / task / "alt_boruta.csv"
    bor = pd.read_csv(bor_f) if bor_f.exists() else None
    d = stable.merge(cmh[["trait", "KO", "q_bh", "mh_odds_ratio", "direction"]], on=["trait", "KO"], how="left")
    d = d.merge(cmh_phy, on=["trait", "KO"], how="left")
    d["cmh_significant"] = d.q_bh < 0.05
    d["sign_agrees"] = np.sign(d.shap_direction) == d.direction
    if bor is not None:
        d = d.merge(bor[["trait", "KO", "n_folds_confirmed", "boruta_confirmed"]], on=["trait", "KO"], how="left")
        d["boruta_confirmed"] = d.boruta_confirmed.fillna(False).astype(bool)
    else:
        d["boruta_confirmed"] = np.nan
    crit = d.cmh_significant & d.sign_agrees & (d.boruta_confirmed if bor is not None else True)
    d["status"] = np.where(crit, "robust", np.where(d.cmh_significant & d.sign_agrees, "supported (2 of 3)",
                                                    "model-specific"))
    # Addendum C: agreement across model classes (tree ensemble vs neural network)
    import os
    sfx = "_zero_baseline" if os.environ.get("NN_ATTRIB") == "ig_zero" else ""
    nn_f = OUT / task / f"xai_stable_nn{sfx}.csv"
    if nn_f.exists():
        nn = pd.read_csv(nn_f)[["trait", "KO", "ig_direction"]]
        d = d.merge(nn.assign(in_nn_stable=True), on=["trait", "KO"], how="left")
        d["in_nn_stable"] = d.in_nn_stable.fillna(False).astype(bool)
        d["robust_across_model_classes"] = d.in_nn_stable & (np.sign(d.ig_direction) == np.sign(d.shap_direction))
        rows = []
        for t in d.trait.unique():
            rf_set = set(d.loc[d.trait == t, "KO"])
            nn_set = set(nn.loc[nn.trait == t, "KO"])
            if nn_set:
                rows.append(dict(trait=t, n_rf_stable=len(rf_set), n_nn_stable=len(nn_set),
                                 n_shared=len(rf_set & nn_set), share_of_nn_in_rf=len(rf_set & nn_set) / len(nn_set)))
        ov = pd.DataFrame(rows)
        ov.to_csv(OUT / task / f"rf_nn_overlap{sfx}.csv", index=False, float_format="%.3f")
        if len(ov):
            print(f"C-P1: median share of NN-stable KOs also RF-stable = {ov.share_of_nn_in_rf.median():.2f} "
                  f"over {len(ov)} traits")
    d = d.sort_values(["trait", "status", "mean_rank"])
    d.to_csv(OUT / task / "consensus.csv", index=False, float_format="%.4g")
    n = d.groupby("trait").status.value_counts().unstack(fill_value=0)
    print(f"[{task}] consensus over RF stable top-10 KOs:")
    print(n.to_string())
    share = (d.cmh_significant & d.sign_agrees).mean()
    print(f"A-P1: share of RF-stable KOs that are also CMH-significant with the same sign: {share:.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["cmh", "boruta", "consensus"])
    ap.add_argument("--task", choices=tx.TASKS, required=True)
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--adjust-size", action="store_true", help="A1b: strata = phylum x genome KO-count quintile")
    a = ap.parse_args()
    (OUT / a.task).mkdir(parents=True, exist_ok=True)
    {"cmh": lambda: cmd_cmh(a.task, a.adjust_size), "boruta": lambda: cmd_boruta(a.task, a.jobs),
     "consensus": lambda: cmd_consensus(a.task)}[a.cmd]()


if __name__ == "__main__":
    main()
