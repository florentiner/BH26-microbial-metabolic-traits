"""RF version of 24_label_concordance.py (final model RedoxKO = random forest).
Cross-database label concordance = an empirical ceiling for this task (audit round 3, M4).

Two independent curations of the same organisms: the provided electron donor/acceptor database (our labels) and
Madin et al. 2020. For training genomes whose species is also in Madin:
  1. per-trait agreement, Cohen's kappa, and the disagreement rate among species called positive by either source;
  2. the released model's out-of-fold predictions (held-out genus) scored against BOTH label sets, compared
     like-for-like with one curation predicting the other.
Audit round 4 fixes: (a) only species with an actual Madin `pathways` record are used (a missing record is not a
curated negative; 53/232 were previously scored as negatives); (b) the model is binarised at its own MCC threshold
before being compared with a binary curation, because average precision of a two-valued score is not comparable with
that of a ranked score (the continuous value is reported too); (c) Madin's per-species `data_source` is reported, so
the FAPROTAX-derived share is visible; (d) acceptor_oxygen is scored as well.
Interpretation: a model trained on curation X is not bounded by how well curation Y predicts X. The comparable
quantities are *model vs the independent curation* and *curation vs curation* - both measure how well an outside
source reproduces a given curation. Both label sets are closed-world, so disagreement mixes genuine differences
(definitions, strain variation) with annotation incompleteness.
Output: results/label_concordance_madin.csv, results/qc/label_concordance_summary.json
"""
import importlib
import json

import numpy as np
import pandas as pd
from sklearn.metrics import cohen_kappa_score

from redox_common import RES, ROOT, load_labelled, oof_matrix, safe_ap

prep = importlib.import_module("18_external_prepare")


def main():
    lab, traits, Y = load_labelled()
    d = pd.read_csv(ROOT / "data/raw/madin/condensed_species_NCBI.csv", low_memory=False)
    d = d[d.pathways.notna()].copy()  # a missing pathways record is not a curated negative
    paths = d.pathways.fillna("").str.split(",").apply(lambda xs: [x.strip() for x in xs])
    for t, keys in prep.MAP.items():
        d[t] = paths.apply(lambda ps: float(any(p == k or p.startswith(k + "_") for p in ps for k in keys)))
    met = d.metabolism.fillna("")
    d["acceptor_oxygen"] = np.where(met.isin(prep.AEROBIC), 1.0, np.where(met.isin(prep.ANAEROBIC), 0.0, np.nan))
    d["key"] = d.species.str.lower().str.strip()
    d = d.drop_duplicates("key").set_index("key")

    lab = lab.copy()
    lab["key"] = lab.species.str.lower().str.strip()
    m = lab.key.isin(d.index).to_numpy()
    print(f"training genomes with a Madin species record: {m.sum()} / {len(lab)}")
    P = oof_matrix(RES / "trees_xai/da/oof", "RF", "genus", lab, traits)

    thr = {}
    bundle_file = RES / "models/RF/model.joblib"
    if bundle_file.exists():
        import joblib
        thr = {t: b["threshold"] for t, b in joblib.load(bundle_file).items()}
    rows = []
    for t in list(prep.MAP) + ["acceptor_oxygen"]:
        if t == "nitrification":
            continue
        j = traits.index(t)
        other = lab.loc[m, "key"].map(d[t]).to_numpy(dtype=float)
        ours = Y[m, j]
        ok = ~np.isnan(other)
        a, b = ours[ok], other[ok]
        both_pos = ((a == 1) | (b == 1))
        p_cont = P[m][ok, j]
        p_bin = (p_cont >= thr.get(t, 0.5)).astype(float)
        fap = lab.loc[m, "key"].map(d.data_source.fillna("").str.contains("faprotax")).to_numpy()[ok]
        rows.append(dict(trait=t, n_species=int(ok.sum()), n_pos_ours=int(a.sum()), n_pos_madin=int(b.sum()),
                         frac_madin_species_with_faprotax_source=float(np.mean(fap)),
                         agreement=float((a == b).mean()),
                         kappa=float(cohen_kappa_score(a, b)) if len(set(a)) > 1 and len(set(b)) > 1 else np.nan,
                         disagreement_among_either_positive=float((a != b)[both_pos].mean()) if both_pos.any() else np.nan,
                         auprc_madin_predicting_ours=safe_ap(a, b), auprc_ours_predicting_madin=safe_ap(b, a),
                         model_auprc_vs_our_labels_continuous=safe_ap(a, p_cont),
                         model_auprc_vs_madin_labels_continuous=safe_ap(b, p_cont),
                         model_auprc_vs_our_labels_binarised=safe_ap(a, p_bin),
                         model_auprc_vs_madin_labels_binarised=safe_ap(b, p_bin),
                         kappa_model_vs_madin=float(cohen_kappa_score(b, p_bin)) if len(set(p_bin)) > 1 else np.nan))
    out = pd.DataFrame(rows)
    out.to_csv(RES / "additional/rf_paper/rf_label_concordance.csv", index=False, float_format="%.4f")
    out["evaluable"] = (out.n_pos_ours >= 10) & (out.n_pos_madin >= 10)
    ev = out[out.evaluable]
    ev_noh2 = ev[ev.trait != "donor_hydrogen_oxidation"]
    summary = {
        "n_training_genomes_with_madin_pathway_record": int(m.sum()),
        "traits_evaluable": ev.trait.tolist(),
        "mean_frac_madin_species_with_faprotax_source": float(ev.frac_madin_species_with_faprotax_source.mean()),
        "mean_agreement": float(ev.agreement.mean()), "mean_kappa": float(ev.kappa.mean()),
        "mean_kappa_excluding_H2_definitional_mismatch": float(ev_noh2.kappa.mean()),
        "mean_disagreement_among_either_positive": float(ev.disagreement_among_either_positive.mean()),
        "like_for_like_binary": {
            "madin_predicting_our_labels": float(ev.auprc_madin_predicting_ours.mean()),
            "our_labels_predicting_madin": float(ev.auprc_ours_predicting_madin.mean()),
            "model_vs_our_labels_binarised": float(ev.model_auprc_vs_our_labels_binarised.mean()),
            "model_vs_madin_labels_binarised": float(ev.model_auprc_vs_madin_labels_binarised.mean()),
            "mean_kappa_model_vs_madin": float(ev.kappa_model_vs_madin.mean()),
        },
        "continuous_scores": {
            "model_vs_our_labels": float(ev.model_auprc_vs_our_labels_continuous.mean()),
            "model_vs_madin_labels": float(ev.model_auprc_vs_madin_labels_continuous.mean()),
        },
        "interpretation": "Comparable quantities are model-vs-independent-curation and curation-vs-curation, both "
                          "binarised. A model trained on one curation is not bounded by how well another curation "
                          "predicts it; what these numbers bound is how well any outside source reproduces a given "
                          "curation. Closed-world on both sides.",
    }
    (RES / "qc").mkdir(parents=True, exist_ok=True)
    (RES / "additional/rf_paper/rf_label_concordance_summary.json").write_text(json.dumps(summary, indent=2))
    print(out.round(3).to_string())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
