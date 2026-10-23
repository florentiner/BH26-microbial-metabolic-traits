"""Experiments U and S (docs/preregistration_U_usage_and_stacking.md): top-k screening precision (U1), decision-curve
net benefit (U2), selective prediction with conformal labels (U3), and a learned stack of RedoxKO and METABOLIC (S).
Output: results/additional/usage_stacking/results.json (+ net_benefit.csv)
"""
import json
import math
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import matthews_corrcoef, roc_auc_score
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from metabolic_calls import calls_from_rows  # noqa: E402
from redox_common import apply_calibrator, best_mcc_threshold, load_labelled  # noqa: E402

OUT = ROOT / "results/additional/usage_stacking"
OUT.mkdir(parents=True, exist_ok=True)
SEED, B = 20261017, 2000
b2 = joblib.load(ROOT / "results/models/B2/model.joblib")
logit = lambda p: np.log(np.clip(p, 1e-6, 1 - 1e-6) / (1 - np.clip(p, 1e-6, 1 - 1e-6)))


def ci(a):
    return [float(x) for x in np.nanpercentile(a, [2.5, 97.5])]


def boot(d, cluster, fn):
    g = {k: x for k, x in d.groupby(cluster)}
    keys = np.array(list(g))
    rng = np.random.default_rng(SEED)
    return np.array([fn(pd.concat([g[k] for k in rng.choice(keys, len(keys))])) for _ in range(B)])


# ---------------------------------------------------------------- data
p = pd.read_csv(ROOT / "results/external_validation_predictions.csv").query("arm == 'B2'")
lab_e = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t").set_index("genome_id")
Me = calls_from_rows(pd.read_parquet(ROOT / "data/additional/metabolic_hmm/external.parquet"))
p = p[p.genome_id.isin(Me.index)]
ext = pd.DataFrame({"category": p.trait.to_numpy(), "y": p.y.astype(int).to_numpy(), "s": p.score.to_numpy(),
                    "c_M": [int(Me.at[g, c]) for g, c in zip(p.genome_id, p.trait)],
                    "cl": lab_e.loc[p.genome_id, "genus"].to_numpy()})


def from_calls(f):
    x = pd.read_csv(ROOT / f)
    return pd.DataFrame({"category": x.category, "y": x.expected.astype(int), "s": x.score_B2,
                         "c_M": x.call_METABOLIC.astype(int), "cl": x.paper})


new, hold = from_calls("results/additional/x1_calls.csv"), from_calls("results/additional/hybrid/holdout_calls.csv")
for d in (ext, new, hold):
    d["c_R"] = (d.s >= d.category.map({t: b["threshold"] for t, b in b2.items()})).astype(int)
    d["p_cal"] = np.nan
    for t in d.category.unique():
        m = d.category == t
        d.loc[m, "p_cal"] = apply_calibrator(b2[t]["calibrator"], d.loc[m, "s"].to_numpy())
res = {}


# ---------------------------------------------------------------- U1 top-k precision
def prec_at_k(g, k):
    k = min(k, len(g))
    r = g.sort_values("s", ascending=False, kind="mergesort").head(k).y.mean()
    nm, tpm = int(g.c_M.sum()), int(g.y[g.c_M == 1].sum())
    if k <= nm:
        m = tpm / nm
    else:
        neg = g[g.c_M == 0]
        m = (tpm + (k - nm) * (neg.y.mean() if len(neg) else 0)) / k
    return r, m


def u1(d, kfun, min_n=1):
    rs, ms = [], []
    for t, g in d.groupby("category"):
        if len(g) < min_n or g.y.sum() == 0:
            continue
        r, m = prec_at_k(g, kfun(len(g)))
        rs.append(r)
        ms.append(m)
    return float(np.mean(rs)), float(np.mean(ms))


r, m = u1(ext, lambda n: 20)
bs = boot(ext, "cl", lambda s: np.subtract(*u1(s, lambda n: 20)))
res["U1_external_k20"] = {"RedoxKO": r, "METABOLIC": m, "diff": r - m, "ci": ci(bs), "success": bool(ci(bs)[0] > 0)}
kq = lambda n: math.ceil(n / 4)
r, m = u1(new, kq, 8)
bs = boot(new, "cl", lambda s: np.subtract(*u1(s, kq, 8)))
res["U1_newspecies_topquarter"] = {"RedoxKO": r, "METABOLIC": m, "diff": r - m, "ci": ci(bs)}

# ---------------------------------------------------------------- U2 decision curve
PT = np.round(np.arange(0.05, 0.91, 0.05), 2)


def nb(d):
    out = []
    for pt in PT:
        w = pt / (1 - pt)
        rr, mm = [], []
        for t, g in d.groupby("category"):
            n = len(g)
            cr = g.p_cal >= pt
            rr.append(((cr & (g.y == 1)).sum() - w * (cr & (g.y == 0)).sum()) / n)
            mm.append((((g.c_M == 1) & (g.y == 1)).sum() - w * ((g.c_M == 1) & (g.y == 0)).sum()) / n)
        out.append(np.mean(rr) - np.mean(mm))
    return np.array(out)


rows = []
for name, d in (("external", ext), ("new species", new)):
    pt_d = nb(d)
    bs = boot(d, "cl", nb)
    lo, hi = np.nanpercentile(bs, 2.5, axis=0), np.nanpercentile(bs, 97.5, axis=0)
    for i, pt in enumerate(PT):
        rows.append({"set": name, "p_t": pt, "nb_diff": pt_d[i], "ci_low": lo[i], "ci_high": hi[i],
                     "winner": "RedoxKO" if lo[i] > 0 else ("METABOLIC" if hi[i] < 0 else "tie")})
nbdf = pd.DataFrame(rows)
nbdf.to_csv(OUT / "net_benefit.csv", index=False)
res["U2"] = {s: {w: [float(x) for x in g[g.winner == w].p_t] for w in ("RedoxKO", "METABOLIC", "tie")}
             for s, g in nbdf.groupby("set")}

# ---------------------------------------------------------------- U3 selective prediction
q = json.loads((ROOT / "results/additional/x3_quantiles.json").read_text())


def single(d):
    d = d[d.category.isin(q)].copy()
    in1 = np.array([(1 - s) <= q[t]["1"] for t, s in zip(d.category, d.s)])
    in0 = np.array([s <= q[t]["0"] for t, s in zip(d.category, d.s)])
    d["lab"] = np.where(in1 & ~in0, 1, np.where(in0 & ~in1, 0, -1))
    return d


def u3(d):
    s = d[d.lab >= 0]
    return np.array([(s.lab == s.y).mean() - (s.c_M == s.y).mean(), len(s) / len(d)])


for name, d in (("new species", new), ("external", ext)):
    s = single(d)
    ss = s[s.lab >= 0]
    bs = boot(s, "cl", u3)
    res[f"U3_{name}"] = {"covered": float(len(ss) / len(s)), "acc_RedoxKO_confident": float((ss.lab == ss.y).mean()),
                         "acc_METABOLIC_same_tests": float((ss.c_M == ss.y).mean()),
                         "diff": float((ss.lab == ss.y).mean() - (ss.c_M == ss.y).mean()), "ci": ci(bs[:, 0]),
                         "success": bool(ci(bs[:, 0])[0] > 0)}

# ---------------------------------------------------------------- S stacking
lab, traits, Y = load_labelled()
oof = pd.read_parquet(ROOT / "results/oof/B2.parquet")
oof = oof[oof.level == "genus"].groupby(["trait", "genome_id"]).prob.mean()
Mt = calls_from_rows(pd.read_parquet(ROOT / "data/processed/genome_features_metabolic_hmm.parquet")).loc[lab.genome_id]
groups = lab.group_genus.to_numpy()
stack = {}
for t in Mt.columns:
    j = traits.index(t)
    Xs = np.c_[logit(oof.loc[t].reindex(lab.genome_id).to_numpy()), Mt[t].to_numpy()]
    y = Y[:, j]
    o = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(Xs, y, groups):
        o[te] = LogisticRegression(C=1.0, class_weight="balanced").fit(Xs[tr], y[tr]).predict_proba(Xs[te])[:, 1]
    lr = LogisticRegression(C=1.0, class_weight="balanced").fit(Xs, y)
    stack[t] = {"model": lr, "threshold": best_mcc_threshold(y, o), "coef": lr.coef_[0].tolist()}
joblib.dump(stack, OUT / "stack_models.joblib")


def apply_stack(d):
    d = d.copy()
    d["s_S"], d["c_S"] = d.s, d.c_R
    for t, st in stack.items():
        m = d.category == t
        if m.any():
            pr = st["model"].predict_proba(np.c_[logit(d.loc[m, "s"].to_numpy()), d.loc[m, "c_M"].to_numpy()])[:, 1]
            d.loc[m, "s_S"], d.loc[m, "c_S"] = pr, (pr >= st["threshold"]).astype(int)
    return d


def macro(d, f):
    v = [f(g) for _, g in d.groupby("category") if g.y.nunique() == 2]
    return float(np.mean(v))


auc = lambda col: (lambda g: roc_auc_score(g.y, g[col]))
mccf = lambda col: (lambda g: matthews_corrcoef(g.y, g[col]))
ns, hs, es = apply_stack(new), apply_stack(hold), apply_stack(ext)
pooled = lambda s, a, b: matthews_corrcoef(s.y, s[a]) - matthews_corrcoef(s.y, s[b])
res["S_coefficients"] = {t: st["coef"] for t, st in stack.items()}
for name, d in (("HS1_primary", ns), ("HS2_withheld", hs)):
    bs = boot(d, "cl", lambda s: [pooled(s, "c_S", "c_M"), pooled(s, "c_S", "c_R"),
                                  macro(s, auc("s_S")) - macro(s, auc("s"))])
    res[name] = {"pooled_mcc_stack": float(matthews_corrcoef(d.y, d.c_S)),
                 "pooled_mcc_METABOLIC": float(matthews_corrcoef(d.y, d.c_M)),
                 "pooled_mcc_RedoxKO": float(matthews_corrcoef(d.y, d.c_R)),
                 "d_vs_METABOLIC": pooled(d, "c_S", "c_M"), "ci_vs_METABOLIC": ci(bs[:, 0]),
                 "d_vs_RedoxKO": pooled(d, "c_S", "c_R"), "ci_vs_RedoxKO": ci(bs[:, 1]),
                 "HS4_d_macro_auroc_vs_RedoxKO": macro(d, auc("s_S")) - macro(d, auc("s")), "HS4_ci": ci(bs[:, 2]),
                 "sens_stack": float((d.c_S[d.y == 1] == 1).mean()), "spec_stack": float((d.c_S[d.y == 0] == 0).mean())}
res["HS1_primary"]["success"] = bool(res["HS1_primary"]["ci_vs_METABOLIC"][0] > 0)
res["HS2_withheld"]["same_sign"] = bool(np.sign(res["HS2_withheld"]["d_vs_METABOLIC"]) ==
                                        np.sign(res["HS1_primary"]["d_vs_METABOLIC"]))
bs = boot(es, "cl", lambda s: [macro(s, mccf("c_S")) - macro(s, mccf("c_M")), macro(s, mccf("c_S")) - macro(s, mccf("c_R")),
                               macro(s, auc("s_S")) - macro(s, auc("s"))])
res["HS3_external"] = {"macro_mcc_stack": macro(es, mccf("c_S")), "macro_mcc_METABOLIC": macro(es, mccf("c_M")),
                       "macro_mcc_RedoxKO": macro(es, mccf("c_R")),
                       "d_vs_METABOLIC": macro(es, mccf("c_S")) - macro(es, mccf("c_M")), "ci_vs_METABOLIC": ci(bs[:, 0]),
                       "d_vs_RedoxKO": macro(es, mccf("c_S")) - macro(es, mccf("c_R")), "ci_vs_RedoxKO": ci(bs[:, 1]),
                       "HS4_d_macro_auroc_vs_RedoxKO": macro(es, auc("s_S")) - macro(es, auc("s")), "HS4_ci": ci(bs[:, 2])}
(OUT / "results.json").write_text(json.dumps(res, indent=2, default=float))
print(json.dumps(res, indent=1, default=lambda x: round(float(x), 4)))
