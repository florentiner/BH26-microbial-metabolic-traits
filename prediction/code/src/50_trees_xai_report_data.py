"""Collect tree / XAI / correlation / rules results into one JSON for the report, and score every pre-registered
prediction of docs/preregistration_trees_xai.md. Missing pieces are marked "pending", never guessed.

usage: python src/50_trees_xai_report_data.py
Output: results/trees_xai/report_data.json, results/trees_xai/predictions_scored.csv
"""
import gzip
import json

import numpy as np
import pandas as pd

from pathway_common import RES

OUT = RES / "trees_xai"
TASK_NAME = {"da": "Electron donors & acceptors (358 genomes, 14 traits)",
             "ox": "Oxygen use (2 927 genomes, 1 trait)",
             "pw": "Energy-metabolism pathways (3 378 genomes, 20 traits)"}
LEARNED = ("RF", "ET", "NN", "L1", "DT", "FIGS", "B2", "B2XL", "PB1", "PB2", "PB3", "PB4")


def ko_names():
    with gzip.open(RES.parent / "data/raw/kofam/ko_list.gz", "rt") as f:
        d = pd.read_csv(f, sep="\t", usecols=["knum", "definition"])
    return dict(zip(d.knum, d.definition.str.slice(0, 70)))


def rd(path):
    return pd.read_csv(path) if path.exists() else None


def rec(df, n=None):
    if df is None:
        return []
    df = df.head(n) if n else df
    return json.loads(df.replace({np.inf: None, -np.inf: None}).where(pd.notna(df), None).to_json(orient="records"))


def task_block(task, names):
    d = OUT / task
    mm, pt, pa = rd(d / "metrics_macro.csv"), rd(d / "metrics_per_trait.csv"), rd(d / "paired.csv")
    if mm is None:
        return {"status": "pending"}
    g = mm[mm.split_level == "genus"]
    learned = g[g.arm.isin(LEARNED)].sort_values("macro_auprc", ascending=False)
    best = learned.arm.iloc[0] if len(learned) else None
    out = {"status": "done", "name": TASK_NAME[task], "best_model": best,
           "macro": rec(mm.sort_values(["split_level", "macro_auprc"], ascending=[True, False])),
           "paired": rec(pa)}
    if pt is not None and best:
        fam = pt[(pt.split_level == "family") & (pt.arm == best)].sort_values("auprc", ascending=False)
        gen = pt[(pt.split_level == "genus") & (pt.arm == best)].set_index("trait")
        dt = pt[(pt.split_level == "family") & (pt.arm == "DT")].set_index("trait")
        nul = pt[(pt.split_level == "family") & (pt.arm == "null_taxNN")].set_index("trait")
        rows = []
        for r in fam.itertuples():
            rows.append(dict(trait=r.trait, n_pos=r.n_pos, prevalence=r.prevalence, auprc_family=r.auprc,
                             auprc_genus=gen.at[r.trait, "auprc"] if r.trait in gen.index else None,
                             mcc=r.mcc, mcc_ci_low=r.mcc_ci_low, threshold=r.threshold_median,
                             q_vs_null=r.q_vs_null, predictable=bool(r.predictable),
                             dt_auprc=dt.at[r.trait, "auprc"] if r.trait in dt.index else None,
                             null_auprc=nul.at[r.trait, "auprc"] if r.trait in nul.index else None))
        out["per_trait_best"] = rows
        out["n_predictable_best"] = int(fam.predictable.sum())
    # XAI of the forest: stable genes with names, direction and validation flags
    st = rd(d / "xai_stable_rf.csv")
    if st is not None:
        st["name"] = st.KO.map(names)
        out["xai_rf_stable"] = rec(st.sort_values(["trait", "mean_rank"]))
    for m in ("rf", "nn", "dt", "l1", "et", "figs"):
        e = rd(d / f"xai_enrichment_{m}.csv")
        if e is not None:
            out[f"enrichment_{m}"] = rec(e)
    rules = rd(d / "xai_rules_dt.csv")
    if rules is not None:
        agg = (rules.groupby(["trait", "rule"]).agg(folds=("fold", "nunique"), support=("heldout_support", "sum"),
                                                     precision=("heldout_precision", "mean"),
                                                     recall=("heldout_recall", "mean"), depth=("depth", "first"))
               .reset_index())
        agg = agg[(agg.folds >= 3) & (agg.support >= 5)].sort_values(["trait", "precision"], ascending=[True, False])
        out["dt_rules_recurring"] = rec(agg.groupby("trait").head(3))
    for f in ("rules_vs_best", "surrogate_fidelity", "surrogate_rules", "consensus", "rf_nn_overlap"):
        x = rd(d / f"{f}.csv")
        if x is not None:
            if "rule" in x.columns and f == "surrogate_rules":
                x = x.sort_values(["trait", "support"], ascending=[True, False]).groupby("trait").head(3)
            if f == "consensus":
                x["name"] = x.KO.map(names)
            out[f] = rec(x)
    for f, key in (("alt_cmh.csv", "cmh_phylum"), ("alt_cmh_size_adjusted.csv", "cmh_size")):
        x = rd(d / f)
        if x is not None:
            s = x[x.q_bh < 0.05]
            out[key] = dict(n_tests=int(len(x)), n_significant=int(len(s)), n_positive=int((s.direction > 0).sum()),
                            n_negative=int((s.direction < 0).sum()))
    cs = rd(d / "alt_cmh_size_adjusted.csv")
    if cs is not None:
        top = cs[cs.q_bh < 0.05].sort_values("p").groupby("trait").head(5).copy()
        top["name"] = top.KO.map(names)
        out["cmh_top"] = rec(top)
    fr = rd(RES / "rulebook" / f"{task}_forest_rules.csv")
    if fr is not None:
        out["forest_rules"] = rec(fr[fr["rank"] <= 2])
    ovz = rd(d / "rf_nn_overlap_zero_baseline.csv")
    if ovz is not None:
        out["rf_nn_overlap_zero"] = rec(ovz)
    b = rd(d / "alt_boruta.csv")
    if b is not None:
        out["boruta"] = dict(n_confirmed=int(b.boruta_confirmed.sum()),
                             per_trait=b[b.boruta_confirmed].groupby("trait").size().to_dict())
    return out


SINGLE_MARKER = {"da": ["donor_ammonia_oxidation", "acceptor_sulfate"],          # T2 / D-P2, as pre-registered
                 "pw": ["methanogenesis", "sulfate_reduction"], "ox": []}
HC = {"K02274", "K02275", "K02276", "K00404", "K00405", "K00406"}
BD = {"K00425", "K00426"}


def score_predictions(blocks):
    """Every pre-registered clause for every task it covers, verdicts binary as specified (audit M2/M3): 'supported'
    or 'not supported'; 'pending' only while its inputs do not exist. Post-hoc rows carry preregistered=False and are
    kept out of the tally."""
    rows = []

    def add(pid, text, verdict, evidence, prereg=True):
        rows.append(dict(id=pid, prediction=text, verdict=verdict, evidence=evidence, preregistered=prereg))

    def predictable(task):
        return [r["trait"] for r in blocks.get(task, {}).get("per_trait_best", []) if r["predictable"]]

    for task in ("da", "ox", "pw"):
        b = blocks.get(task, {})
        done = b.get("status") == "done"
        m = pd.DataFrame(b.get("macro", []))
        g = m[m.split_level == "genus"].set_index("arm") if len(m) else pd.DataFrame()
        # T1 (DA and PW only)
        if task in ("da", "pw"):
            ref = {"da": ("B2", "B2XL"), "pw": ("PB1", "PB2")}[task]
            have = [r for r in ref if r in g.index]
            if "RF" in g.index and have:
                best_ref = max(have, key=lambda a: g.at[a, "macro_auprc"])
                dlt = g.at["RF", "macro_auprc"] - g.at[best_ref, "macro_auprc"]
                add(f"T1 ({task})", "RF within ±0.03 macro AUPRC of the best existing arm",
                    "supported" if abs(dlt) <= 0.03 else "not supported", f"RF − {best_ref} = {dlt:+.3f} (genus)")
            else:
                add(f"T1 ({task})", "RF within ±0.03 macro AUPRC of the best existing arm", "pending", "")
        # T2, both clauses, every task
        if "RF" in g.index and "DT" in g.index:
            dlt = g.at["DT", "macro_auprc"] - g.at["RF", "macro_auprc"]
            c1 = dlt <= -0.05
            ev = f"DT − RF = {dlt:+.3f} (genus)"
            c2, ev2 = True, []
            oof_f = OUT / task / "oof" / "DT3.parquet"     # depth-2/3 trees, same nested CV and cut-off rule
            if SINGLE_MARKER[task] and oof_f.exists():
                o = pd.read_parquet(oof_f)
                o = o[o.level == "genus"]
                for t in SINGLE_MARKER[task]:
                    ot = o[o.trait == t]
                    if not len(ot):
                        continue
                    call = ot.prob >= ot.threshold
                    prec = float((ot.y[call] == 1).mean()) if call.any() else float("nan")
                    c2 &= bool(prec >= 0.7)
                    ev2.append(f"{t}: precision {prec:.2f} ({int(call.sum())} calls)")
                ev += "; second clause (depth ≤ 3 tree, held-out precision ≥ 0.7): " + "; ".join(ev2)
            elif SINGLE_MARKER[task]:
                ev += "; second clause not yet evaluated (depth ≤ 3 tree not run)"
                c2 = None
            add(f"T2 ({task})", "DT loses ≥ 0.05 to RF, but a depth ≤ 3 tree reaches precision ≥ 0.7 for "
                "single-marker traits", "pending" if c2 is None and c1 else
                ("supported" if c1 and c2 else "not supported"), ev)
        else:
            add(f"T2 ({task})", "DT loses ≥ 0.05 to RF; depth ≤ 3 tree precision ≥ 0.7 for single-marker traits",
                "pending", "")
        # A-P2, every task
        if all(a in g.index for a in ("RF", "FIGS", "DT")):
            a1 = g.at["FIGS", "macro_auprc"] - g.at["DT", "macro_auprc"]
            a2 = g.at["FIGS", "macro_auprc"] - g.at["RF", "macro_auprc"]
            add(f"A-P2 ({task})", "FIGS beats DT by ≥ 0.03 and is within 0.05 of RF",
                "supported" if a1 >= 0.03 and a2 >= -0.05 else "not supported",
                f"FIGS − DT = {a1:+.3f}; FIGS − RF = {a2:+.3f}")
        else:
            add(f"A-P2 ({task})", "FIGS beats DT by ≥ 0.03 and is within 0.05 of RF", "pending", "")
        pred = predictable(task)
        # T3
        e = pd.DataFrame(b.get("enrichment_rf", []))
        if len(e) and pred:
            ee = e[e.trait.isin(pred)]
            share = float((ee.q_bh < 0.05).mean()) if len(ee) else np.nan
            add(f"T3 ({task})", "stable RF top-10 enriched for known markers in ≥ 70 % of predictable traits",
                "supported" if share >= 0.7 else "not supported",
                f"{int((ee.q_bh < 0.05).sum())} of {len(ee)} predictable traits enriched (q < 0.05)")
        else:
            add(f"T3 ({task})", "stable RF top-10 enriched for known markers in ≥ 70 % of predictable traits",
                "pending", "")
        # A-P1, as pre-registered: stable RF top-10 of predictable traits, A1 (unadjusted CMH) significance
        st_f, cmh_f = OUT / task / "xai_stable_rf.csv", OUT / task / "alt_cmh.csv"
        if st_f.exists() and cmh_f.exists() and pred:
            st, cm = pd.read_csv(st_f), pd.read_csv(cmh_f)
            st = st[st.trait.isin(pred)]
            sig = set(zip(cm[cm.q_bh < 0.05].trait, cm[cm.q_bh < 0.05].KO))
            hit = [(t, k) in sig for t, k in zip(st.trait, st.KO)]
            share = float(np.mean(hit)) if hit else np.nan
            add(f"A-P1 ({task})", "≥ 50 % of RF stable top-10 genes of predictable traits are A1 (CMH) significant",
                "supported" if share >= 0.5 else "not supported", f"{sum(hit)} of {len(hit)} ({share:.0%})")
        else:
            add(f"A-P1 ({task})", "≥ 50 % of RF stable top-10 genes of predictable traits are A1 significant",
                "pending", "")
        # D-P1..3
        rv = pd.DataFrame(b.get("rules_vs_best", []))
        sf = pd.DataFrame(b.get("surrogate_fidelity", []))
        if len(rv):
            dtv = rv[(rv.rule_set == "tree rules (DT)") & rv.predictable]
            if len(dtv):
                share = float((dtv.delta_ci_low > 0).mean())
                add(f"D-P1 ({task})", "best model beats tree rules on MCC for ≥ 70 % of predictable traits",
                    "supported" if share >= 0.7 else "not supported",
                    f"significantly better on {int((dtv.delta_ci_low > 0).sum())} of {len(dtv)} "
                    f"(point estimate better on {int((dtv.delta_mcc_best_minus_rule > 0).sum())})")
            mk = rv[rv.rule_set.str.startswith("marker") & rv.trait.isin(SINGLE_MARKER[task])]
            if len(mk):
                ok = mk.delta_mcc_best_minus_rule <= 0.10
                add(f"D-P2 ({task})", "marker rule within 0.10 MCC of the best model for single-marker traits",
                    "supported" if ok.all() else "not supported",
                    f"{int(ok.sum())} of {len(mk)} within 0.10: " +
                    "; ".join(f"{r.trait}: Δ {r.delta_mcc_best_minus_rule:+.2f}" for r in mk.itertuples()))
        elif task != "ox" or not done:
            add(f"D-P1/D-P2 ({task})", "rules versus the best model", "pending", "")
        if len(sf):
            s2 = sf[sf.predictable]
            if len(s2):
                ok = s2.fidelity_accuracy >= 0.90
                add(f"D-P3 ({task})", "depth-3 surrogate fidelity ≥ 0.90 for predictable traits",
                    "supported" if ok.all() else "not supported",
                    f"{int(ok.sum())} of {len(s2)} traits ≥ 0.90 (range {s2.fidelity_accuracy.min():.2f}–"
                    f"{s2.fidelity_accuracy.max():.2f})")
        elif not len(rv):
            add(f"D-P3 ({task})", "depth-3 surrogate fidelity ≥ 0.90 for predictable traits", "pending", "")
        # C-P1 (NN tasks), pre-registered on predictable traits with the all-absent IG baseline
        zb = OUT / task / "rf_nn_overlap_zero_baseline.csv"
        if zb.exists():
            z = pd.read_csv(zb)
            z = z[z.trait.isin(pred)] if pred else z
            med = float(z.share_of_nn_in_rf.median())
            note = ""
            sz, iz = OUT / task / "xai_stable_nn_zero_baseline.csv", OUT / task / "xai_importance_nn_zero_baseline.csv"
            if sz.exists() and iz.exists():   # measured per task, not assumed
                prev = pd.read_csv(iz).groupby("KO").prevalence_in_train.mean()
                mp = float(prev.reindex(pd.read_csv(sz).KO).median())
                note = (f" — its stable genes are near-universal (median prevalence {mp:.2f}): the all-absent baseline "
                        f"favours housekeeping genes regardless of trait")
            add(f"C-P1 ({task})", "≥ 40 % of NN-stable genes also RF-stable, predictable traits (IG, all-absent "
                "baseline)", "supported" if med >= 0.4 else "not supported",
                f"median {med:.0%} over {len(z)} predictable traits{note}")
        elif task == "pw":
            add("C-P1 (pw)", "≥ 40 % of NN-stable genes also RF-stable (NN on PW only 'if time allows')", "not run", "",
                prereg=False)
        ov = pd.DataFrame(b.get("rf_nn_overlap", []))
        if len(ov):
            ov = ov[ov.trait.isin(pred)] if pred else ov
            med = float(ov.share_of_nn_in_rf.median())
            add(f"C-P1* ({task})", "same, POST HOC with Expected Gradients (training-genome baselines)",
                "supported" if med >= 0.4 else "not supported",
                f"median {med:.0%} over {len(ov)} predictable traits (post hoc, not in the tally)", prereg=False)
    # T4 (ox), as pre-registered: top-3 by mean |SHAP|
    imp_f = OUT / "ox" / "xai_importance_rf.csv"
    if imp_f.exists():
        imp = pd.read_csv(imp_f)
        top3 = list(imp.groupby("KO").mean_abs_shap.mean().sort_values(ascending=False).index[:3])
        add("T4 (ox)", "top-3 SHAP genes include an aa3/cbb3 subunit and not cytochrome bd",
            "supported" if (set(top3) & HC) and not (BD & set(top3)) else "not supported",
            "top-3 by mean |SHAP|: " + ", ".join(top3))
    else:
        add("T4 (ox)", "top-3 SHAP genes include an aa3/cbb3 subunit and not cytochrome bd", "pending", "")
    cm = pd.DataFrame(blocks.get("ox", {}).get("cmh_top", []))
    a1 = OUT / "ox" / "alt_cmh.csv"
    if a1.exists():                                   # A-P3 is on A1, the unadjusted CMH
        c = pd.read_csv(a1)
        pos = c[(c.direction > 0) & (c.q_bh < 0.05)].sort_values("p")
        top = pos.KO.iloc[0] if len(pos) else None
        bd_p = pos[pos.KO.isin(BD)].p.min() if len(pos) else np.nan
        hc_p = pos[pos.KO.isin(HC)].p.min() if len(pos) else np.nan
        ok = top in HC and (np.isnan(bd_p) or bd_p > hc_p)
        add("A-P3 (ox)", "strongest positive A1 (CMH) gene is an aa3/cbb3 subunit; bd weaker",
            "supported" if ok else "not supported", f"strongest positive: {top}")
    # T5
    st = pd.DataFrame(blocks.get("da", {}).get("xai_rf_stable", []))
    if len(st):
        hao = st[(st.trait == "donor_ammonia_oxidation") & (st.KO == "K10535")]
        add("T5 (da)", "hao (K10535) not among the stable RF top-10 for ammonia oxidation",
            "supported" if hao.empty else "not supported",
            "hao absent from the stable set" if hao.empty else f"hao stability {hao.stability.iloc[0]:.1f}")
    return pd.DataFrame(rows)


def main():
    names = ko_names()
    blocks = {t: task_block(t, names) for t in ("da", "ox", "pw")}
    preds = score_predictions(blocks)
    preds.to_csv(OUT / "predictions_scored.csv", index=False)
    (OUT / "report_data.json").write_text(json.dumps({"tasks": blocks, "predictions": rec(preds)}, indent=1,
                                                     default=float))
    print({t: b.get("status") for t, b in blocks.items()})
    print(preds.to_string(index=False))


if __name__ == "__main__":
    main()
