"""Collect the rulebook into one JSON for the docx and score predictions G-P1 / G-P2 (Addendum G).

usage: python src/52_rulebook_report_data.py
Output: results/rulebook/rulebook_data.json
"""
import json

import numpy as np
import pandas as pd

from pathway_common import RES

OUT = RES / "rulebook"
TN = {"da": "Electron donors & acceptors", "ox": "Oxygen use", "pw": "Energy-metabolism pathways"}
# known explanations for specific external failures, verified in the data (see the report text)
FAILURE_NOTES = {
    ("donor_ammonia_oxidation", "K10945"): "3 of the 4 lineage-novel external genomes carrying amoB are methanotrophs "
    "(Methylomarinovum, Methylomarinum, Methyloglobulus): their methane monooxygenase has the same KEGG ID as ammonia "
    "monooxygenase.",
    ("donor_sulfur_oxidation", "K27925"): "only 2 external genomes; one is Acidiferrobacter, whose type species "
    "A. thiooxydans oxidises sulfur - probably a missing label rather than a rule error.",
}


def rec(df):
    return json.loads(df.replace({np.inf: None, -np.inf: None}).where(pd.notna(df), None).to_json(orient="records"))


def main():
    import importlib
    ext_x, _ = importlib.import_module("51_rulebook").load_external()
    data = {"tasks": {}, "predictions": [], "n_external_novel": int(len(ext_x))}
    for task in ("da", "ox", "pw"):
        f = OUT / f"{task}_rules.csv"
        if not f.exists():
            data["tasks"][task] = {"status": "pending", "name": TN[task]}
            continue
        r = pd.read_csv(f)
        notes = []
        for row in r.itertuples():
            n = ""
            for (t, ko), txt in FAILURE_NOTES.items():
                if row.trait == t and ko in row.rule:
                    n = txt
            notes.append(n)
        r["failure_note"] = notes
        r["probability"] = r.external_precision.where(r.external_precision.notna(), r.heldout_precision) \
            if "external_precision" in r.columns else r.heldout_precision
        r["probability_source"] = np.where(r.get("external_precision", pd.Series(np.nan, index=r.index)).notna(),
                                           "external genomes", "held-out genera")
        block = {"status": "done", "name": TN[task], "rules": rec(r)}
        p = OUT / f"{task}_procedure.csv"
        if p.exists():
            block["procedure"] = rec(pd.read_csv(p))
        tb = OUT / f"{task}_transfer_by_level.csv"
        if tb.exists():
            block["transfer_by_level"] = rec(pd.read_csv(tb))
        # rules the random forest itself uses (Addendum I)
        fr = OUT / f"{task}_forest_rules.csv"
        if fr.exists():
            block["forest_rules"] = rec(pd.read_csv(fr))
            # audit round 2 N2: is the forest's most frequent first question a gene of the trait's own pathway?
            from gene_mechanisms import MECH, PATHWAY
            f1 = pd.read_csv(fr)
            f1 = f1[f1["rank"] == 1].drop_duplicates("trait")
            block["forest_top"] = [dict(trait=r.trait, gene=r.genes, name=MECH.get(r.genes, (r.genes,))[0],
                                        on_pathway=bool(r.genes in PATHWAY.get(r.trait, set())),
                                        tree_frequency=float(r.tree_frequency)) for r in f1.itertuples()]
        # what the forest looks at (XAI) and the recurring single-tree rules
        import gzip
        with gzip.open(RES.parent / "data/raw/kofam/ko_list.gz", "rt") as fh:
            names = dict(pd.read_csv(fh, sep="\t", usecols=["knum", "definition"]).values)
        st = RES / "trees_xai" / task / "xai_stable_rf.csv"
        if st.exists():
            x = pd.read_csv(st)
            x["name"] = x.KO.map(lambda k: str(names.get(k, ""))[:60])
            block["forest_xai"] = rec(x.sort_values(["trait", "mean_rank"]))
        dtr = RES / "trees_xai" / task / "xai_rules_dt.csv"
        if dtr.exists():
            x = pd.read_csv(dtr)
            agg = (x.groupby(["trait", "rule"]).agg(folds=("fold", "nunique"), support=("heldout_support", "sum"),
                                                    precision=("heldout_precision", "mean"), depth=("depth", "first"))
                   .reset_index())
            agg = agg[(agg.folds >= 3) & (agg.support >= 5)].sort_values(["trait", "precision"],
                                                                          ascending=[True, False])
            block["dt_rules"] = rec(agg.groupby("trait").head(3))
        ti = OUT / f"{task}_trusted_interpretation.csv"
        if ti.exists():
            block["trusted"] = rec(pd.read_csv(ti))
        rvb = RES / "trees_xai" / task / "rules_vs_best.csv"
        if rvb.exists():
            block["rules_vs_best"] = rec(pd.read_csv(rvb))
        mt_f = RES / "trees_xai" / task / "metrics_per_trait.csv"
        if mt_f.exists():
            m = pd.read_csv(mt_f)
            best = m[m.split_level == "genus"].groupby("arm").auprc.mean()
            best = best[~best.index.str.startswith("null")].idxmax()
            block["best_model"] = best
            block["best_model_family"] = rec(m[(m.arm == best) & (m.split_level == "family")][
                ["trait", "auprc", "mcc", "mcc_ci_low", "predictable", "prevalence"]])
        data["tasks"][task] = block
        if task == "da":
            mt = pd.read_csv(RES / "trees_xai/da/metrics_per_trait.csv")
            pred = set(mt[(mt.arm == "RF") & (mt.split_level == "family") & mt.predictable].trait)
            pos = r[(r.cls == "has trait") & r.trait.isin(pred)]
            # G-P1 as pre-registered: "high probability" = held-out (or external) precision >= 0.80, Wilson LB >= 0.60;
            # the rule-level flag uses the external value where one exists, the held-out value otherwise
            ok = pos[pos.high_probability].trait.nunique()
            ok_h = pos[(pos.heldout_precision >= 0.8) & (pos.heldout_ci_low >= 0.6)].trait.nunique()
            ok_e = pos[(pos.external_precision >= 0.8) & (pos.external_ci_low >= 0.6)].trait.nunique()
            data["predictions"].append(dict(
                id="G-P1 (da)", prediction="≥ 1 high-probability positive rule for ≥ 5 of the 7 predictable traits",
                verdict="supported" if ok_h >= 5 else "not supported", preregistered=True,
                evidence=f"{ok_h} of {len(pred)} traits on held-out genera (the reading G-P2 implies: 'those rules' "
                         f"are then checked externally); on the external genomes alone {ok_e} of {len(pred)} — "
                         f"the verdict depends on the reading"))
            # G-P2 as pre-registered (audit C2): "those rules" = the G-P1 positive rules of predictable traits,
            # selected on their HELD-OUT precision (never on the external outcome); every one must transfer
            sel = pos[(pos.heldout_precision >= 0.8) & (pos.heldout_ci_low >= 0.6) & pos.external_precision.notna()]
            diff = (sel.external_precision - sel.heldout_precision).abs()
            within = diff <= 0.10
            data["predictions"].append(dict(
                id="G-P2 (da)", prediction="external precision of those rules within 0.10 of held-out precision",
                verdict="supported" if len(sel) and within.all() else "not supported", preregistered=True,
                evidence=f"{int(within.sum())} of {len(sel)} rules within 0.10 (median |difference| "
                         f"{float(diff.median()):.2f}; rules selected on held-out precision, external genomes "
                         f"restricted to novel lineages)"))
    # I-P1: the most frequent forest path shares a gene with the rulebook rules of the same trait
    for task in ("da", "ox", "pw"):
        fr = OUT / f"{task}_forest_rules.csv"
        mt = RES / "trees_xai" / task / "metrics_per_trait.csv"
        if fr.exists() and mt.exists():
            d, m = pd.read_csv(fr), pd.read_csv(mt)
            best = m[m.split_level == "genus"].groupby("arm").auprc.mean()
            best = best[~best.index.str.startswith("null")].idxmax()
            pred = set(m[(m.arm == best) & (m.split_level == "family") & m.predictable].trait)
            first = d[(d["rank"] == 1) & d.trait.isin(pred)]
            if len(first):
                n = int(first.shares_gene_with_rulebook.sum())
                data["predictions"].append(dict(
                    id=f"I-P1 ({task})", prediction="the forest's most frequent path shares a gene with the rulebook",
                    verdict="supported" if n == len(first) else "not supported", preregistered=True,
                    evidence=f"{n} of {len(first)} predictable traits"))
    (OUT / "rulebook_data.json").write_text(json.dumps(data, indent=1, default=float))
    print({t: b["status"] for t, b in data["tasks"].items()})
    for p in data["predictions"]:
        print(p)


if __name__ == "__main__":
    main()
