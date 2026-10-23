"""Where, if anywhere, is RedoxKO better than the random forest on the same 852-KO panel?

Compared on the properties the paper relies on, with the same data as the main results:
  1. calibration of the probabilities (Brier score, expected calibration error with 10 bins), genus-level CV
  2. rare categories (< 10% positives) and per-category wins, genus- and family-level CV, external, new species
  3. complementarity with METABOLIC: ranking of METABOLIC's positive calls on new species (AUROC within them)
  4. stability: external macro AUPRC with three different random-seed sets for the forest
  5. model size and prediction time for 100 genomes
usage: python src/140_rf_vs_redoxko.py
Output: results/additional/rf_vs_redoxko.json
"""
import importlib
import io
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import RES, best_mcc_threshold, fit_models, inner_folds, load_features, load_labelled, oof_matrix, safe_ap  # noqa

tx = importlib.import_module("47_trees_xai")
tx.TREE_JOBS = 8
OUT = RES / "additional/rf_vs_redoxko.json"


def ece(y, p, bins=10):
    e = np.linspace(0, 1, bins + 1)
    idx = np.clip(np.digitize(p, e) - 1, 0, bins - 1)
    return float(sum(abs(y[idx == b].mean() - p[idx == b].mean()) * (idx == b).mean() for b in range(bins) if (idx == b).any()))


def main():
    res = {}
    lab, traits, Y = load_labelled()
    prev = np.nanmean(Y, axis=0)
    rare = [j for j in range(len(traits)) if prev[j] < 0.10]
    obs = ~np.isnan(Y)
    # 1-2: cross-validation
    for level in ("genus", "family"):
        Pb = oof_matrix(RES / "oof", "B2", level, lab, traits)
        Pr = oof_matrix(RES / "trees_xai/da/oof", "RF", level, lab, traits)
        r = {}
        for name, P in (("RedoxKO", Pb), ("RF", Pr)):
            ap = np.array([safe_ap(Y[obs[:, j], j], P[obs[:, j], j]) for j in range(len(traits))])
            y, p = Y[obs], P[obs]
            r[name] = dict(macro_auprc=float(np.nanmean(ap)), rare_auprc=float(np.nanmean(ap[rare])),
                           brier=float(np.mean((p - y) ** 2)), ece=ece(y, p), per_cat=ap.tolist())
        r["categories_RedoxKO_better"] = int(sum(a > b for a, b in zip(r["RedoxKO"]["per_cat"], r["RF"]["per_cat"])))
        res[f"cv_{level}"] = r
    # 2: external
    ext_b2 = pd.read_csv(RES / "external_validation_predictions.csv")
    ext_b2 = ext_b2[ext_b2.arm == "B2"]
    ext_rf = pd.read_csv(RES / "external_validation_missing_rows_predictions.csv")
    ext_rf = ext_rf[ext_rf.arm == "Random forest"]
    m = ext_b2.merge(ext_rf[["genome_id", "trait", "score"]], on=["genome_id", "trait"], suffixes=("_b2", "_rf"))
    pc = {}
    for t, g in m.groupby("trait"):
        pc[t] = dict(positives=int(g.y.sum()), share=float(g.y.mean()), RedoxKO=safe_ap(g.y.to_numpy(), g.score_b2.to_numpy()),
                     RF=safe_ap(g.y.to_numpy(), g.score_rf.to_numpy()),
                     brier_RedoxKO=float(np.mean((g.score_b2 - g.y) ** 2)), brier_RF=float(np.mean((g.score_rf - g.y) ** 2)))
    res["external_per_category"] = pc
    res["external_categories_RedoxKO_better"] = int(sum(v["RedoxKO"] > v["RF"] for v in pc.values()))
    res["external_brier"] = dict(RedoxKO=float(np.mean((m.score_b2 - m.y) ** 2)), RF=float(np.mean((m.score_rf - m.y) ** 2)))
    res["external_ece"] = dict(RedoxKO=ece(m.y.to_numpy(), m.score_b2.to_numpy()), RF=ece(m.y.to_numpy(), m.score_rf.to_numpy()))
    # 2-3: new species
    d = pd.read_csv(RES / "additional/rf_new_species_scores.csv")
    pcn = {}
    for t, g in d.groupby("category"):
        if 0 < g.expected.sum() < len(g):
            pcn[t] = dict(tests=len(g), RedoxKO=roc_auc_score(g.expected, g.score), RF=roc_auc_score(g.expected, g.score_rf))
    res["new_species_per_category_auroc"] = pcn
    res["new_species_categories_RedoxKO_better"] = int(sum(v["RedoxKO"] > v["RF"] for v in pcn.values()))
    res["new_species_sens_spec"] = {n: dict(sensitivity=float((d[c][d.expected == 1] == 1).mean()),
                                            specificity=float((d[c][d.expected == 0] == 0).mean()))
                                    for n, c in (("RedoxKO", "call"), ("RF", "call_rf"))}
    met = pd.read_csv(ROOT / "validation_new_organisms/novel100_model_scores.csv") \
        if (ROOT / "validation_new_organisms/novel100_model_scores.csv").exists() else None
    tools = RES / "additional/metabolic_bench"
    mc = None
    for f in [RES / "additional/x1_calls.csv"]:
        if f.exists():
            mc = pd.read_csv(f)
            break
    res["metabolic_calls_file_found"] = mc is not None
    if mc is not None:
        col = [c for c in mc.columns if c.lower() == "call_metabolic"][0]
        x = d.merge(mc[["accession", "category", col]], on=["accession", "category"])
        x = x[x[col] == 1]
        res["within_metabolic_positive_auroc"] = dict(n=len(x), RedoxKO=roc_auc_score(x.expected, x.score),
                                                      RF=roc_auc_score(x.expected, x.score_rf))
    # 4-5: seed stability, size and speed on the external panel
    X = load_features("energy_panel", lab.genome_id)
    keep = X.columns[X.std(axis=0) > 0]
    Xi = X[keep].to_numpy()
    j = traits.index("acceptor_sulfate")
    y = Y[:, j]
    t0 = time.time()
    lr = fit_models({"model": "logreg", "C": 0.1}, Xi, y)
    size_lr = len(pickle.dumps(lr))
    rf = tx.fit_all({"model": "rf", "msl": 1}, Xi, y)
    size_rf = len(pickle.dumps(rf))
    Xv = np.repeat(Xi[:1], 100, axis=0)
    t1 = time.time(); [m_.predict_proba(Xv) for m_ in lr]; t_lr = time.time() - t1
    t1 = time.time(); [m_.predict_proba(Xv) for m_ in rf]; t_rf = time.time() - t1
    res["one_category_model"] = dict(bytes_logreg=size_lr, bytes_rf=size_rf, seconds_100_genomes_logreg=t_lr,
                                     seconds_100_genomes_rf=t_rf)
    ext = m.copy()
    seeds_sets = [(0, 100, 200), (1, 101, 201), (2, 102, 202)]
    ext_lab = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t").set_index("genome_id")
    Xe = pd.read_csv(ROOT / "data/external/ko_matrix_panel.tsv", sep="\t", index_col=0).reindex(columns=keep).fillna(0)
    stab = []
    for ss in seeds_sets:
        tx.SEEDS = ss
        aps = []
        for t in ["acceptor_sulfate", "acceptor_sulfur", "acceptor_nitrate", "donor_sulfur_oxidation", "acceptor_oxygen"]:
            jj = traits.index(t)
            mm = ~np.isnan(Y[:, jj])
            s = tx.fit_predict({"model": "rf", "msl": 1}, Xi[mm], Y[mm, jj], Xe.to_numpy())
            yy = ext_lab.reindex(Xe.index)[t].to_numpy(dtype=float)
            ok = ~np.isnan(yy)
            aps.append(safe_ap(yy[ok], s[ok]))
        stab.append(float(np.mean(aps)))
    res["rf_seed_sets_external_auprc_5cats"] = stab
    OUT.write_text(json.dumps(res, indent=1, default=float))
    for k in ("cv_genus", "cv_family"):
        print(k, {n: {kk: round(v, 3) for kk, v in res[k][n].items() if kk != "per_cat"} for n in ("RedoxKO", "RF")},
              "RedoxKO better in", res[k]["categories_RedoxKO_better"], "of 14")
    print("external brier", res["external_brier"], "ece", res["external_ece"], "RedoxKO better in",
          res["external_categories_RedoxKO_better"], "of", len(pc))
    print("new species RedoxKO better in", res["new_species_categories_RedoxKO_better"], "of", len(pcn), res["new_species_sens_spec"])
    print("within METABOLIC positives", res.get("within_metabolic_positive_auroc"), "file found", res["metabolic_calls_file_found"])
    print("size/speed", res["one_category_model"])
    print("RF seed sets", stab)


if __name__ == "__main__":
    main()
