"""New-species numbers of the paper recomputed for the final model RedoxKO (random forest, results/models/RF).

Scores every growth experiment of the primary 100 species (results/additional/x1_calls.csv, with the other tools' calls),
the 43 withheld species (results/additional/hybrid/holdout_calls.csv) and all scored species (scored178_B2.csv) from
the 852-KO panels, with the final model's thresholds and released categories. Recomputes the old model (score_B2,
now RedoxKO-LX) from the same files as a check against the published tables.
usage: python src/144_rf_new_species_numbers.py
Output: results/additional/rf_paper/new_species.json and *_rf.csv copies of the input tables with score_RF/call_RF
"""
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
V, RES = ROOT / "validation_new_organisms", ROOT / "results"
OUT = RES / "additional/rf_paper"
SEED, B = 20260921, 2000


def rf_scores(accs):
    bundle = joblib.load(RES / "models/RF/model.joblib")
    meta = json.loads((RES / "models/RF/metadata.json").read_text())
    K = pd.concat([pd.read_parquet(V / "ko_panel_validation.parquet"), pd.read_parquet(V / "ko_panel_auto.parquet")])
    K = (K[~K.index.duplicated()] > 0).astype(np.float32)
    missing = sorted(set(accs) - set(K.index))
    assert not missing, f"no KO panel for {missing[:5]}"
    K = K.loc[sorted(set(accs))]
    S = {t: pd.Series(np.asarray(b["model"].predict_proba(K.reindex(columns=b["model"].columns).fillna(0))), index=K.index)
         for t, b in bundle.items()}
    thr = {t: float(b["threshold"]) for t, b in bundle.items()}
    return S, thr, meta["reliable_traits_by_level"]


def add_rf(d, S, thr, rel):
    d = d.copy()
    d["score_RF"] = [S[c][a] for a, c in zip(d.accession, d.category)]
    d["call_RF"] = [int(s >= thr[c]) for s, c in zip(d.score_RF, d.category)]
    if "family_in_training" in d:
        d["deployed_RF"] = [c in rel["genus" if f else "family"] for c, f in zip(d.category, d.family_in_training)]
    return d


def summary(d, score, call):
    au = [roc_auc_score(g.expected, g[score]) for _, g in d.groupby("category") if 0 < g.expected.sum() < len(g)]
    y, c = d.expected.to_numpy(), d[call].to_numpy()
    return dict(auroc=float(np.mean(au)), accuracy=float((y == c).mean()), mcc=float(matthews_corrcoef(y, c)),
                sensitivity=float(c[y == 1].mean()), specificity=float(1 - c[y == 0].mean()),
                ppv=float(y[c == 1].mean()), npv=float(1 - y[c == 0].mean()), n=int(len(d)))


def boot(d, fn, n=B, seed=SEED):
    rng = np.random.default_rng(seed)
    papers = d.paper.unique()
    idx = {p: np.where(d.paper.to_numpy() == p)[0] for p in papers}
    out = []
    for _ in range(n):
        x = d.iloc[np.concatenate([idx[p] for p in rng.choice(papers, len(papers))])]
        try:
            out.append(fn(x))
        except ValueError:
            pass
    return np.nanpercentile(out, [2.5, 97.5]).tolist()


def mcc(y, c):
    return float(matthews_corrcoef(y, c))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    x1 = pd.read_csv(RES / "additional/x1_calls.csv")
    ho = pd.read_csv(RES / "additional/hybrid/holdout_calls.csv")
    sc = pd.read_csv(V / "scored178_B2.csv")
    S, thr, rel = rf_scores(set(x1.accession) | set(ho.accession) | set(sc.accession))
    x1, ho, sc = add_rf(x1, S, thr, rel), add_rf(ho, S, thr, rel), add_rf(sc, S, thr, rel)
    x1.to_csv(OUT / "x1_calls_rf.csv", index=False)
    ho.to_csv(OUT / "holdout_calls_rf.csv", index=False)
    sc.to_csv(OUT / "scored178_rf.csv", index=False)
    r = {"released_categories": rel}
    # Table 3 (primary 100 species, 255 experiments)
    for name, s, c in (("RedoxKO", "score_RF", "call_RF"), ("RedoxKO-LX", "score_B2", "call_B2")):
        r[f"t3_{name}"] = summary(x1, s, c)
        r[f"t3_{name}"]["auroc_ci"] = boot(x1, lambda x: summary(x, s, c)["auroc"])
    # text: calls in released categories only
    dep = x1[x1.deployed_RF]
    r["released_only_RF"] = dict(calls=int(len(dep)), accuracy=float((dep.expected == dep.call_RF).mean()))
    iron = x1[x1.category == "acceptor_iron"]
    r["iron_RF"] = dict(tests=int(len(iron)), errors=int((iron.expected != iron.call_RF).sum()),
                        released="acceptor_iron" in rel["family"])
    # all 143 primary-analysis species, and without small categories
    p143 = sc[sc.blind & ~sc.category.isin(["donor_other", "acceptor_other"])]
    r["all143"] = {n: summary(p143, s, c)["auroc"] for n, s, c in (("RedoxKO", "score_RF", "call_RF"), ("RedoxKO-LX", "score", "call"))}
    big = x1.groupby("category").filter(lambda g: len(g) >= 10)
    r["t3_no_small_categories_auroc"] = {"RedoxKO": summary(big, "score_RF", "call_RF")["auroc"]}
    # Table 4 and tools text
    y = x1.expected
    r["t4"] = {"RedoxKO_mcc": mcc(y, x1.call_RF), "RedoxKO_auroc": r["t3_RedoxKO"]["auroc"],
               "METABOLIC_minus_RedoxKO": mcc(y, x1.call_METABOLIC) - mcc(y, x1.call_RF),
               "METABOLIC_minus_RedoxKO_ci": boot(x1, lambda x: mcc(x.expected, x.call_METABOLIC) - mcc(x.expected, x.call_RF))}
    relc = x1[x1.category.isin(rel["family"])]
    r["t4_released_categories"] = {"categories": rel["family"], "tests": int(len(relc)),
                                   "RedoxKO": mcc(relc.expected, relc.call_RF), "METABOLIC": mcc(relc.expected, relc.call_METABOLIC),
                                   "DRAM": mcc(relc.expected, relc.call_DRAM.fillna(0)),
                                   "METABOLIC_minus_RedoxKO_ci": boot(relc, lambda x: mcc(x.expected, x.call_METABOLIC) - mcc(x.expected, x.call_RF)),
                                   "DRAM_minus_RedoxKO_ci": boot(relc, lambda x: mcc(x.expected, x.call_DRAM.fillna(0)) - mcc(x.expected, x.call_RF))}
    md = x1[x1.call_MicrobeDecoder.notna() & x1.score_MicrobeDecoder.notna()]
    r["t4_microbe_decoder"] = {"categories": sorted(md.category.unique()), "tests": int(len(md)),
                               "MD_mcc": mcc(md.expected, md.call_MicrobeDecoder), "RedoxKO_mcc": mcc(md.expected, md.call_RF),
                               "MD_auroc": summary(md, "score_MicrobeDecoder", "call_MicrobeDecoder")["auroc"],
                               "RedoxKO_auroc": summary(md, "score_RF", "call_RF")["auroc"]}
    ox = x1[(x1.category == "acceptor_oxygen") & x1.score_BacDiveAI.notna()]
    r["t4_oxygen"] = {"tests": int(len(ox)), "BacDive_auroc": roc_auc_score(ox.expected, ox.score_BacDiveAI),
                      "RedoxKO_auroc": roc_auc_score(ox.expected, ox.score_RF),
                      "RedoxKO_mcc": mcc(ox.expected, ox.call_RF)}
    # within METABOLIC's positive calls (primary and withheld)
    for name, d in (("primary", x1), ("withheld", ho)):
        m = d[d.call_METABOLIC == 1]
        r[f"within_metabolic_positive_{name}"] = {
            "n": int(len(m)), "RedoxKO": roc_auc_score(m.expected, m.score_RF),
            "RedoxKO_ci": boot(m, lambda x: roc_auc_score(x.expected, x.score_RF), n=3000, seed=20261017)}
    # withheld 43 species
    r["withheld"] = {"species": int(ho.accession.nunique()), "tests": int(len(ho)),
                     **{k: v for k, v in summary(ho, "score_RF", "call_RF").items() if k in ("auroc", "mcc")},
                     "METABOLIC_mcc": mcc(ho.expected, ho.call_METABOLIC),
                     "METABOLIC_minus_RedoxKO_ci": boot(ho, lambda x: mcc(x.expected, x.call_METABOLIC) - mcc(x.expected, x.call_RF))}
    # per-category balanced accuracy (Fig. 3a)
    pc = {}
    for c, g in x1.groupby("category"):
        def bacc(col):
            yy, cc = g.expected.to_numpy(), g[col].to_numpy()
            parts = [cc[yy == 1].mean() if (yy == 1).any() else np.nan, 1 - cc[yy == 0].mean() if (yy == 0).any() else np.nan]
            return float(np.nanmean(parts))
        pc[c] = dict(tests=int(len(g)), RedoxKO=bacc("call_RF"), METABOLIC=bacc("call_METABOLIC"))
    r["per_category_balanced_accuracy"] = pc
    (OUT / "new_species.json").write_text(json.dumps(r, indent=1, default=float))
    print(json.dumps({k: r[k] for k in ("t3_RedoxKO", "t3_RedoxKO-LX", "released_only_RF", "iron_RF", "all143", "t4",
                                        "t4_released_categories", "t4_microbe_decoder", "t4_oxygen",
                                        "within_metabolic_positive_primary", "within_metabolic_positive_withheld",
                                        "withheld")}, indent=1, default=lambda v: round(float(v), 3)))


if __name__ == "__main__":
    main()
