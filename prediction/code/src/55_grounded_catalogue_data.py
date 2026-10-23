"""Data for the grounded rule catalogue: every rule judged on two axes and given one verdict.

Axis 1 - trust (statistics): precision on the genomes it fires on (external genomes where available, otherwise
held-out genera), its 95 % interval, how many genomes it fires on, and the share of trait carriers it catches.
Axis 2 - biological grounding (computed in src/54_rule_interpretation.py): does every gene connect to the trait
through known biology?
Verdict:
  Use                          trusted and grounded
  Use with care                trusted, partly grounded (part of the rule reflects lineage)
  Do not rely on it            trusted but not grounded (accurate here, no biological reason: likely lineage)
  Sound, not reliable enough   grounded but not trusted
  Reject                       neither

usage: python src/55_grounded_catalogue_data.py
Output: results/rulebook/grounded_catalogue.json, grounded_catalogue.csv
"""
import importlib
import json
import re

import numpy as np
import pandas as pd

from pathway_common import RES

ri = importlib.import_module("54_rule_interpretation")
OUT = RES / "rulebook"
TN = {"da": "Electron donors & acceptors", "ox": "Oxygen use", "pw": "Energy-metabolism pathways"}
ORDER = ["Use", "Use with care", "Sound, not reliable enough", "Do not rely on it", "Reject"]


def verdict(trusted, grounding):
    if trusted:
        return {"grounded": "Use", "partly grounded": "Use with care"}.get(grounding, "Do not rely on it")
    return "Sound, not reliable enough" if grounding == "grounded" else "Reject"


def rec(df):
    return json.loads(df.replace({np.inf: None, -np.inf: None}).where(pd.notna(df), None).to_json(orient="records"))


def main():
    out, allrows = {"tasks": {}}, []
    out["n_external_novel"] = int(len(importlib.import_module("51_rulebook").load_external()[0]))
    for task in ("da", "ox", "pw"):
        f = OUT / f"{task}_all_interpretation.csv"
        if not f.exists():
            out["tasks"][task] = {"status": "pending", "name": TN[task]}
            continue
        d = pd.read_csv(f)
        _, trs, Yt, Xt = ri.tx.load_task(task)
        ri.register_task(Xt, trs, Yt)          # forest-rule grounding below uses carrier shares (audit round 2 N3)
        d["verdict"] = [verdict(t, g) for t, g in zip(d.trusted, d.grounding)]
        d["verdict_rank"] = d.verdict.map({v: i for i, v in enumerate(ORDER)})
        d = d.sort_values(["trait", "verdict_rank", "probability"], ascending=[True, True, False])
        allrows.append(d.assign(source="rulebook"))
        block = {"status": "done", "name": TN[task], "rules": rec(d),
                 "counts": d.verdict.value_counts().reindex(ORDER, fill_value=0).to_dict()}
        # the forest's own single-gene questions, judged the same way
        fr_f = OUT / f"{task}_forest_rules.csv"
        if fr_f.exists():
            fr = pd.read_csv(fr_f)
            rows = []
            for r in fr.itertuples():
                conds = [(k, s == "present") for k, s in re.findall(r"\((K\d{5})\) (present|absent)", r.rule)]
                cls = "has trait" if r.direction == "raises" else "lacks trait"
                gv, _, _, why = ri.grounding(conds, cls, r.trait)
                use_ext = pd.notna(getattr(r, "p_external", np.nan))    # probability, base rate and test: same genomes
                p = r.p_external if use_ext else r.p_heldout
                base = r.external_base_rate if use_ext else r.base_rate
                inf = r.informative_external if use_ext else r.informative_heldout
                # a forest rule is "trusted" when it is informative and decisive: >= 80 % for raising rules,
                # <= 20 % (i.e. >= 80 % correct at ruling out) for lowering rules
                decisive = (p >= 0.8) if cls == "has trait" else (p <= 0.2)
                n_sup = r.external_support if use_ext else r.heldout_support
                trusted = bool(inf) and bool(decisive) and n_sup >= 8     # same minimum as rulebook rules (round 2 n4)
                rows.append(dict(trait=r.trait, rule=r.rule, cls=cls, p_trait=p, base_rate=base, trusted=trusted,
                                 p_source="external genomes" if use_ext else "held-out genera", support=n_sup,
                                 grounding=gv, grounding_reasons=why, tree_frequency=r.tree_frequency,
                                 stability=r.stability, verdict=verdict(trusted, gv)))
            fdf = pd.DataFrame(rows)
            fdf["verdict_rank"] = fdf.verdict.map({v: i for i, v in enumerate(ORDER)})
            block["forest"] = rec(fdf.sort_values(["trait", "verdict_rank"]))
        out["tasks"][task] = block
    if allrows:
        a = pd.concat(allrows)
        a.to_csv(OUT / "grounded_catalogue.csv", index=False, float_format="%.4f")
        out["overall"] = a.verdict.value_counts().reindex(ORDER, fill_value=0).to_dict()
        # does grounding predict transfer? (donor/acceptor rules with an external test)
        da = pd.read_csv(OUT / "da_rules.csv") if (OUT / "da_rules.csv").exists() else None
        if da is not None and "external_precision" in da.columns:
            g = pd.read_csv(OUT / "da_all_interpretation.csv")[["rule", "grounding"]]
            m = da.merge(g, on="rule")
            m = m[m.external_precision.notna() & m.heldout_precision.notna() & (m.cls == "has trait")]
            tr = m.groupby("grounding").agg(n=("rule", "size"), heldout=("heldout_precision", "mean"),
                                            external=("external_precision", "mean")).reset_index()
            out["transfer_by_grounding"] = rec(tr)
    (OUT / "grounded_catalogue.json").write_text(json.dumps(out, indent=1, default=float))
    print({t: b.get("counts") for t, b in out["tasks"].items()})
    print("overall:", out.get("overall"))
    print("transfer by grounding:", out.get("transfer_by_grounding"))


if __name__ == "__main__":
    main()
