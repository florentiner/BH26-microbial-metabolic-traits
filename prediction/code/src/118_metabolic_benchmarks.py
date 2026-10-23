"""Experiment M (docs/preregistration_M_metabolic_benchmarks.md), parts M1, M2 and M4: RedoxKO versus METABOLIC on the
external 562 genomes and per category on the 100 primary new species; speed from logged runs.
`validate` re-checks the METABOLIC call mapping against real METABOLIC output (143 genomes).
Outputs: results/additional/metabolic_bench/{m1_m2_m4.json, scorecard.csv}
"""
import importlib
import json
import re
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from metabolic_calls import calls_from_rows  # noqa: E402

m95 = importlib.import_module("95_metabolic_score")
OUT = ROOT / "results/additional/metabolic_bench"
OUT.mkdir(parents=True, exist_ok=True)
SEED, B = 20261017, 2000


def validate():
    tot, agree = 0, 0
    R = pd.read_parquet(ROOT / "data/additional/metabolic_hmm/newspecies.parquet")
    for W in [m95.W2, m95.W2.parents[2] / "holdout_out/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet2.tsv"]:
        w = m95.function_table(W)
        pres = [c for c in w.columns if c.endswith("Function presence")]
        w[pres] = w[pres].eq("Present")
        for col in pres:
            g = col.replace(" Function presence", "")
            C = calls_from_rows(R.loc[[g]]).iloc[0]
            for cat, pairs in m95.MAP.items():
                m = w.Function.str.strip().isin([f for _, f in pairs])
                tot += 1
                agree += int(C[cat] == int(w.loc[m, col].any()))
            ox = w.Category.str.startswith("Oxygen metabolism") & w.Function.str.contains("|".join(m95.OXY_KEYS))
            tot += 1
            agree += int(C["acceptor_oxygen"] == int(w.loc[ox, col].any()))
    print(f"METABOLIC call mapping vs real METABOLIC: {agree}/{tot}")
    return agree, tot


def mcc(y, c):
    return float(matthews_corrcoef(y, c)) if len(set(y)) == 2 else np.nan


def per_cat(d):
    rows = []
    for t, g in d.groupby("category"):
        if g.y.nunique() < 2:
            continue
        r = {"category": t, "n": len(g), "pos": int(g.y.sum())}
        for tool, s, c in (("RedoxKO", "s_R", "c_R"), ("METABOLIC", "c_M", "c_M")):
            y = g.y.to_numpy()
            r[f"auroc_{tool}"] = roc_auc_score(y, g[s])
            r[f"mcc_{tool}"] = mcc(y, g[c])
            r[f"sens_{tool}"] = float((g[c][y == 1] == 1).mean())
            r[f"spec_{tool}"] = float((g[c][y == 0] == 0).mean())
        rows.append(r)
    return pd.DataFrame(rows)


def stats(d, cats=None):
    s = d if cats is None else d[d.category.isin(cats)]
    pc = per_cat(s)
    y = s.y.to_numpy()
    out = {f"macro_auroc_{t}": float(pc[f"auroc_{t}"].mean()) for t in ("RedoxKO", "METABOLIC")}
    out |= {f"macro_mcc_{t}": float(pc[f"mcc_{t}"].mean()) for t in ("RedoxKO", "METABOLIC")}
    for t, c in (("RedoxKO", "c_R"), ("METABOLIC", "c_M")):
        out[f"sens_{t}"] = float((s[c][y == 1] == 1).mean())
        out[f"spec_{t}"] = float((s[c][y == 0] == 0).mean())
        out[f"prec_{t}"] = float((y[s[c] == 1] == 1).mean())
    return out, pc


def boot(d, cluster, fn, cats=None):
    groups = {k: g for k, g in d.groupby(cluster)}
    keys = np.array(list(groups))
    rng = np.random.default_rng(SEED)
    res = []
    for _ in range(B):
        s = pd.concat([groups[k] for k in rng.choice(keys, len(keys))])
        res.append(fn(s))
    return np.array(res)


def ci(a):
    return [float(x) for x in np.nanpercentile(a, [2.5, 97.5])]


def main():
    agree, tot = validate()
    rel = [t for t, b in joblib.load(ROOT / "results/models/B2/model.joblib").items() if b.get("deployed", False)]
    released = sorted(set(pd.read_csv(ROOT / "results/additional/x1_calls.csv").query("deployed").category))
    thr = {t: b["threshold"] for t, b in joblib.load(ROOT / "results/models/B2/model.joblib").items()}
    # ---------------- external
    p = pd.read_csv(ROOT / "results/external_validation_predictions.csv").query("arm == 'B2'")
    lab = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t").set_index("genome_id")
    M = calls_from_rows(pd.read_parquet(ROOT / "data/additional/metabolic_hmm/external.parquet"))
    p = p[p.genome_id.isin(M.index)].copy()
    d = pd.DataFrame({"genome": p.genome_id, "category": p.trait, "y": p.y.astype(int), "s_R": p.score,
                      "genus": lab.loc[p.genome_id, "genus"].to_numpy()})
    d["c_R"] = (d.s_R >= d.category.map(thr)).astype(int)
    d["c_M"] = [int(M.at[g, c]) for g, c in zip(d.genome, d.category)]
    res = {"mapping_validation": f"{agree}/{tot}", "external": {"n_pairs": int(len(d)), "n_genomes": int(d.genome.nunique()),
                                                                 "categories": sorted(d.category.unique())}}
    full, pc_ext = stats(d)
    rel_ext = [c for c in released if c in set(d.category)]
    relst, _ = stats(d, rel_ext)
    fa = lambda s: stats(s)[0]["macro_auroc_RedoxKO"] - stats(s)[0]["macro_auroc_METABOLIC"]

    def both(s, cats=None):
        o, _ = stats(s, cats)
        return [o["macro_auroc_RedoxKO"] - o["macro_auroc_METABOLIC"], o["macro_mcc_RedoxKO"] - o["macro_mcc_METABOLIC"]]
    bs = boot(d, "genus", both)
    bsr = boot(d, "genus", lambda s: both(s, rel_ext))
    res["external"] |= {"all": full, "M1_d_macro_auroc": full["macro_auroc_RedoxKO"] - full["macro_auroc_METABOLIC"],
                        "M1_ci": ci(bs[:, 0]), "M1_success": bool(ci(bs[:, 0])[0] > 0),
                        "M1b_d_macro_mcc": full["macro_mcc_RedoxKO"] - full["macro_mcc_METABOLIC"], "M1b_ci": ci(bs[:, 1]),
                        "released_categories": rel_ext, "released": relst,
                        "M1b_released_d_macro_mcc": relst["macro_mcc_RedoxKO"] - relst["macro_mcc_METABOLIC"],
                        "M1b_released_ci": ci(bsr[:, 1])}
    # ---------------- M2 scorecard: per-category differences with intervals
    x1 = pd.read_csv(ROOT / "results/additional/x1_calls.csv")
    n = pd.DataFrame({"genome": x1.accession, "category": x1.category, "y": x1.expected.astype(int), "s_R": x1.score_B2,
                      "c_R": x1.call_B2.astype(int), "c_M": x1.call_METABOLIC.astype(int), "paper": x1.paper})
    rows = []
    for name, dd, cl in (("external", d, "genus"), ("new species", n, "paper")):
        _, pc = stats(dd)
        for t in pc.category:
            g = dd[dd.category == t]

            def f(s):
                if s.y.nunique() < 2:
                    return [np.nan] * 4
                pr = per_cat(s).iloc[0]
                return [pr[f"{m}_RedoxKO"] - pr[f"{m}_METABOLIC"] for m in ("auroc", "mcc", "sens", "spec")]
            bb = boot(g, cl, f)
            pr = pc[pc.category == t].iloc[0].to_dict()
            for k, m in enumerate(("auroc", "mcc", "sens", "spec")):
                lo, hi = ci(bb[:, k])
                p2 = 2 * min(np.nanmean(bb[:, k] <= 0), np.nanmean(bb[:, k] >= 0))
                rows.append({"set": name, "category": t, "n": pr["n"], "pos": pr["pos"], "metric": m,
                             "RedoxKO": pr[f"{m}_RedoxKO"], "METABOLIC": pr[f"{m}_METABOLIC"],
                             "diff": pr[f"{m}_RedoxKO"] - pr[f"{m}_METABOLIC"], "ci_low": lo, "ci_high": hi,
                             "p_boot": min(1.0, p2),
                             "winner": "RedoxKO" if lo > 0 else ("METABOLIC" if hi < 0 else "tie")})
    sc = pd.DataFrame(rows)
    from statsmodels.stats.multitest import multipletests
    sc["q_bh"] = multipletests(sc.p_boot.fillna(1), method="fdr_bh")[1]
    sc.to_csv(OUT / "scorecard.csv", index=False)
    res["M2_wins"] = sc.groupby(["set", "metric"]).winner.value_counts().unstack(fill_value=0).to_dict(orient="index")
    res["M2_wins"] = {f"{a}|{b}": v for (a, b), v in res["M2_wins"].items()}
    # ---------------- M4 speed
    def wall(log):
        m = re.search(r"total running time: (\d+):(\d+):(\d+)", (ROOT / log).read_text())
        return int(m.group(1)) * 3600 + int(m.group(2)) * 60 + int(m.group(3))
    res["M4_speed"] = {"METABOLIC_primary100_6threads_s_per_genome": wall("data/additional/metabolic/full_run.log") / 100,
                       "METABOLIC_withheld43_8threads_s_per_genome": wall("data/additional/metabolic/holdout_run.log") / 43,
                       "RedoxKO_panel_scan_s_per_genome": 3.3, "RedoxKO_full_prediction_from_accession_s": 14.6}
    (OUT / "m1_m2_m4.json").write_text(json.dumps(res, indent=2, default=float))
    print(json.dumps(res, indent=1, default=lambda x: round(float(x), 4)))


if __name__ == "__main__":
    validate() if sys.argv[1:] == ["validate"] else main()
