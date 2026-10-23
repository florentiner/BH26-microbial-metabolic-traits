"""Rules versus the best model (docs/preregistration_trees_xai.md, Addendum D).

For each task whose best model is trained, three rule sets are compared with it on the same held-out genomes:
tree rules (the DT arm), textbook marker rules, and a depth-3 surrogate tree that imitates the best model.

usage: python src/49_rules_vs_best.py --task {da,ox,pw}
Outputs: results/trees_xai/<task>/rules_vs_best.csv, surrogate_rules.csv, surrogate_fidelity.csv
"""
import argparse
import importlib
import json

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.metrics import cohen_kappa_score, matthews_corrcoef
from sklearn.model_selection import GroupKFold
from sklearn.tree import DecisionTreeClassifier

tx = importlib.import_module("47_trees_xai")
from pathway_common import ko_matrix  # noqa: E402

OUT = tx.OUT
N_BOOT = 1000
SEED = 20260918
AA3_CBB3 = ["K02274", "K02275", "K02276", "K00404", "K00405", "K00406"]
# Marker rules fixed in Addendum D before this comparison was run: trait -> KOs (rule = any present)
MARKER_RULES = {
    "donor_ammonia_oxidation": ["K10945"], "acceptor_sulfate": ["K11181"],
    "acceptor_CO2_methanogen_acetogen": ["K00399", "K14138"], "acceptor_oxygen": AA3_CBB3, "uses_O2": AA3_CBB3,
    "acceptor_nitrate": ["K00370", "K02567"], "acceptor_sulfur": ["K17219", "K08352"],
    "donor_hydrogen_oxidation": ["K06281", "K05922"], "methanogenesis": ["K00399"],
    "sulfate_reduction": ["K11181"], "nitrogen_fixation": ["K02588"],
    "denitrification": ["K00368", "K15864", "K00376"], "acetogenesis": ["K14138"],
    "photoautotrophy": ["K08928", "K08929", "K02689", "K02703"],
    "photoheterotrophy": ["K08928", "K08929", "K02689", "K02703"], "methylotrophy": ["K14028"],
    "hydrogen_oxidation": ["K06281", "K05922"],
}
EXISTING = {"B2": "oof/B2.parquet", "B2XL": "oof/B2XL.parquet", "PB1": "pathways/oof/PB1.parquet",
            "PB2": "pathways/oof/PB2.parquet", "PB3": "pathways/oof/PB3.parquet", "PB4": "pathways/oof/PB4.parquet"}


def best_model(task):
    m = pd.read_csv(OUT / task / "metrics_macro.csv")
    m = m[(m.split_level == "genus") & ~m.arm.str.startswith("null")]
    return m.sort_values("macro_auprc", ascending=False).arm.iloc[0]


def oof_calls(task, arm, lab, traits, level="genus"):
    f = (tx.RES / EXISTING[arm]) if arm in EXISTING else (OUT / task / "oof" / f"{arm}.parquet")
    d = pd.read_parquet(f)
    d = d[(d.level == level) & d.trait.isin(traits)]
    gid = {g: i for i, g in enumerate(lab.genome_id)}
    P = np.full((len(lab), len(traits)), np.nan)
    C = np.full((len(lab), len(traits)), np.nan)
    for j, t in enumerate(traits):
        dt = d[d.trait == t]
        idx = dt.genome_id.map(gid).dropna().astype(int).to_numpy()
        keep = dt.genome_id.isin(gid).to_numpy()
        P[idx, j] = dt.prob.to_numpy()[keep]
        C[idx, j] = (dt.prob.to_numpy()[keep] >= dt.threshold.to_numpy()[keep]).astype(float)
    return P, C


def full_ko(task, lab):
    """KO presence for marker rules: the full profile where it exists (PW/OX), else the DA feature panel."""
    if task == "da":
        return tx.load_task("da")[3]
    X = ko_matrix(lab.genome_id)
    X.index = lab.genome_id.to_numpy()
    return X


def boot_delta_mcc(y, a, b, groups, rng):
    ug = np.unique(groups)
    mem = {g: np.where(groups == g)[0] for g in ug}
    out = []
    for _ in range(N_BOOT):
        i = np.concatenate([mem[g] for g in rng.choice(ug, len(ug), replace=True)])
        if 0 < y[i].sum() < len(i):
            out.append(matthews_corrcoef(y[i], a[i]) - matthews_corrcoef(y[i], b[i]))
    return (np.percentile(out, 2.5), np.percentile(out, 97.5)) if out else (np.nan, np.nan)


def score(y, c):
    tp = int(((c == 1) & (y == 1)).sum())
    fp = int(((c == 1) & (y == 0)).sum())
    fn = int(((c == 0) & (y == 1)).sum())
    return dict(precision=tp / (tp + fp) if tp + fp else np.nan, recall=tp / (tp + fn) if tp + fn else np.nan,
                mcc=matthews_corrcoef(y, c), n_called=int(c.sum()))


def surrogate(X, cols, z, y, groups):
    """Depth-3 tree imitating the model's calls z; fidelity measured on held-out genera."""
    pred = np.zeros(len(z))
    for tr, te in GroupKFold(n_splits=5).split(X, z, groups):
        if len(np.unique(z[tr])) < 2:
            pred[te] = z[tr][0]
            continue
        s = DecisionTreeClassifier(max_depth=3, min_samples_leaf=3, class_weight="balanced", random_state=0)
        pred[te] = s.fit(X[tr], z[tr]).predict(X[te])
    fid = dict(fidelity_accuracy=float((pred == z).mean()),
               fidelity_kappa=float(cohen_kappa_score(z, pred)) if len(np.unique(z)) > 1 else np.nan,
               fidelity_on_model_positives=float((pred[z == 1] == 1).mean()) if z.sum() else np.nan,
               surrogate_mcc_vs_labels=float(matthews_corrcoef(y, pred)))
    full = DecisionTreeClassifier(max_depth=3, min_samples_leaf=3, class_weight="balanced", random_state=0).fit(X, z)
    rules = []
    t = full.tree_

    def walk(node, conds):
        if t.children_left[node] == -1:
            v = t.value[node][0]
            if v[1] > v[0]:
                m = np.ones(len(X), bool)
                for f_, present in conds:
                    m &= (X[:, f_] > 0.5) == present
                rules.append(dict(rule=" AND ".join(f"{cols[f_]} {'present' if p else 'absent'}" for f_, p in conds),
                                  support=int(m.sum()), agrees_with_model=float(z[m].mean()) if m.any() else np.nan,
                                  precision_vs_labels=float(y[m].mean()) if m.any() else np.nan))
            return
        walk(t.children_left[node], conds + [(t.feature[node], False)])
        walk(t.children_right[node], conds + [(t.feature[node], True)])
    walk(0, [])
    return fid, rules


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=tx.TASKS, required=True)
    a = ap.parse_args()
    task = a.task
    lab, traits, Y, Xdf = tx.load_task(task)
    best = best_model(task)
    print(f"[{task}] best model (genus-level macro AUPRC): {best}")
    rng = np.random.default_rng(SEED)
    groups = lab.group_genus.to_numpy()
    Pb, Cb = oof_calls(task, best, lab, traits)
    _, Cdt = oof_calls(task, "DT", lab, traits)
    KO = full_ko(task, lab)
    X, cols = Xdf.to_numpy(dtype=np.float32), list(Xdf.columns)
    predictable = set()
    mt = pd.read_csv(OUT / task / "metrics_per_trait.csv")
    predictable = set(mt[(mt.arm == best) & (mt.split_level == "family") & mt.predictable].trait)
    rows, srules, sfid = [], [], []
    for j, t in enumerate(traits):
        ok = ~np.isnan(Cb[:, j]) & ~np.isnan(Cdt[:, j])
        y, cb, cdt, g = Y[ok, j].astype(int), Cb[ok, j].astype(int), Cdt[ok, j].astype(int), groups[ok]
        rule_sets = {"tree rules (DT)": cdt}
        mk = [k for k in MARKER_RULES.get(t, []) if k in KO.columns]
        if mk:
            rule_sets[f"marker rule ({' or '.join(mk)})"] = KO.loc[lab.genome_id[ok], mk].max(axis=1).to_numpy() \
                .astype(int)
        sb = score(y, cb)
        for name, cr in rule_sets.items():
            sr = score(y, cr)
            lo, hi = boot_delta_mcc(y, cb, cr, g, rng)
            b_ = int(((cr == y) & (cb != y)).sum())   # rule right, model wrong
            c_ = int(((cr != y) & (cb == y)).sum())   # rule wrong, model right
            p_mc = binomtest(b_, b_ + c_, 0.5).pvalue if b_ + c_ else 1.0
            rows.append(dict(task=task, trait=t, best_model=best, predictable=t in predictable, rule_set=name,
                             best_precision=sb["precision"], best_recall=sb["recall"], best_mcc=sb["mcc"],
                             rule_precision=sr["precision"], rule_recall=sr["recall"], rule_mcc=sr["mcc"],
                             delta_mcc_best_minus_rule=sb["mcc"] - sr["mcc"], delta_ci_low=lo, delta_ci_high=hi,
                             kappa_rule_vs_model=cohen_kappa_score(cr, cb) if len(np.unique(cr)) > 1 else np.nan,
                             n_rule_right_model_wrong=b_, n_rule_wrong_model_right=c_, p_mcnemar=p_mc))
        fid, rules = surrogate(X[ok], cols, cb, y, g)
        sfid.append(dict(task=task, trait=t, best_model=best, predictable=t in predictable, **fid))
        for r in rules:
            srules.append(dict(task=task, trait=t, best_model=best, **r))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / task / "rules_vs_best.csv", index=False, float_format="%.4f")
    pd.DataFrame(sfid).to_csv(OUT / task / "surrogate_fidelity.csv", index=False, float_format="%.4f")
    pd.DataFrame(srules).to_csv(OUT / task / "surrogate_rules.csv", index=False, float_format="%.4f")
    show = ["trait", "rule_set", "best_mcc", "rule_mcc", "delta_mcc_best_minus_rule", "delta_ci_low", "delta_ci_high",
            "kappa_rule_vs_model", "p_mcnemar"]
    print(d[show].round(3).to_string(index=False))
    print(pd.DataFrame(sfid)[["trait", "fidelity_accuracy", "fidelity_kappa", "fidelity_on_model_positives",
                              "surrogate_mcc_vs_labels"]].round(3).to_string(index=False))
    json.dump({"task": task, "best_model": best}, open(OUT / task / "rules_vs_best.json", "w"))


if __name__ == "__main__":
    main()
