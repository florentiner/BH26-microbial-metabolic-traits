"""Data for the universal-rules document (Addendum J): rules that hold across phyla, with interpretation.

usage: python src/57_universal_data.py
Output: results/rulebook/universal_data.json
"""
import importlib
import json

import numpy as np
import pandas as pd

from pathway_common import RES

OUT = RES / "rulebook"
TN = {"da": "Electron donors & acceptors", "ox": "Oxygen use", "pw": "Energy-metabolism pathways"}
# why a rule built on these genes can hold across lineages - textbook physiology, used only as interpretation
WHY_GENERAL = {
    "methanogenesis": "Methyl-CoM reductase performs the final, methane-forming step and is found in every methanogen "
                      "lineage; outside methanogens it occurs in anaerobic methane oxidisers (running it in reverse) and "
                      "in archaea that oxidise short alkanes with related enzymes, so its presence indicates "
                      "methanogenesis in most, not all, carriers.",
    "sulfate_reduction": "Dissimilatory sulfite reductase and APS reductase are the conserved core of sulfate "
                         "respiration in every lineage that does it — bacterial and archaeal — because the chemistry "
                         "of reducing sulfate to sulfide has one solution.",
    "acceptor_sulfate": "Dissimilatory sulfite reductase and APS reductase are the conserved core of sulfate "
                        "respiration in every lineage that does it.",
    "nitrogen_fixation": "Nitrogenase is the only enzyme family that fixes N₂; the Mo-nitrogenase NifHDK is conserved across "
                         "diazotrophs from cyanobacteria to archaea (alternative V- and Fe-only nitrogenases exist "
                         "alongside it), so carrying complete NifHDK very likely means nitrogen fixation in any lineage.",
    "photoautotrophy": "Photosynthetic reaction centres and the Calvin cycle (RuBisCO, phosphoribulokinase) are "
                       "conserved wherever light-driven CO₂ fixation occurs.",
    "photoheterotrophy": "Type-II reaction centres (PufLM) are conserved across the anoxygenic phototrophs of several "
                         "phyla.",
    "denitrification": "The NO-forming nitrite reductases (NirK, NirS) and NO/N₂O reductases are committed steps of "
                       "denitrification and are shared across the phyla that denitrify.",
    "uses_O2": "Heme-copper terminal oxidases and the enzymes that build them exist to reduce O₂; aerobic metabolism "
               "has the same core in every phylum.",
    "acceptor_oxygen": "Heme-copper terminal oxidases and the bc1 complex are the conserved core of aerobic "
                       "respiration in every phylum.",
    "aerobic_chemoheterotrophy": "Heme-copper terminal oxidases are the conserved core of aerobic respiration in every "
                                 "phylum.",
    "donor_hydrogen_oxidation": "The [NiFe] uptake hydrogenases HyaB and HybC share one catalytic core across the "
                                "phyla that oxidise H₂.",
    "hydrogen_oxidation": "[NiFe] uptake hydrogenases share one catalytic core across the phyla that oxidise H₂.",
    "donor_sulfur_oxidation": "Sulfide:quinone oxidoreductase and the Sox system are core enzymes of sulfur oxidation in "
                              "both aerobic chemolithotrophs and anaerobic phototrophs.",
    "acceptor_sulfur": "Sulfur-compound reductases are shared across the phyla that respire sulfur.",
    "methylotrophy": "Methanol dehydrogenases (MxaF, XoxF) are the entry enzymes of methylotrophy across phyla.",
    "fermentation": "Pyruvate formate-lyase, lactate and alcohol dehydrogenases are the hallmark enzymes of "
                    "fermentation in every phylum that ferments.",
    "nitrate_reduction": "The respiratory nitrate reductase NarGHI is the conserved enzyme of nitrate respiration in "
                         "every phylum that respires nitrate; requiring the O₂-sensitive 2-oxoglutarate:ferredoxin "
                         "oxidoreductase (kor) to be absent keeps out strict anaerobes, so the rule picks the "
                         "facultative organisms that switch from O₂ to nitrate.",
}


def rec(df):
    return json.loads(df.replace({np.inf: None, -np.inf: None}).where(pd.notna(df), None).to_json(orient="records"))


LIMITED_BY_LABELS = {
    "nitrogen_fixation": "Nitrogenase is the only enzyme family that fixes N₂; a genome with the complete NifHDK complex "
                         "very likely fixes nitrogen. The labels record the trait for only a minority of carriers, "
                         "because nitrogen fixation is rarely tested and an untested trait counts as absent.",
    "denitrification": "NO reductase and N₂O reductase are committed steps of denitrification; organisms carrying both "
                       "denitrify at least partly. Denitrification is under-reported in trait databases.",
    "sulfate_reduction": "DsrAB with APS reductase is the conserved core of sulfate respiration, but the same genes "
                         "run in reverse in some sulfur oxidisers (reverse-Dsr) — so part of the gap is real biology, "
                         "part is missing labels.",
    "photoheterotrophy": "The PufLM reaction centre is the hallmark of anoxygenic phototrophy; almost every labelled "
                         "phototroph carries it, but many carriers are not labelled.",
    "photoautotrophy": "Photosystem II marks oxygenic photosynthesis (cyanobacteria); none of the carriers is labelled "
                       "with this trait in the source database.",
    "acetogenesis": "The Wood–Ljungdahl pathway is also used by sulfate reducers and methanogens, so carrying it does "
                    "not make an organism an acetogen — this gap is biology, not labels.",
    "methylotrophy": "Methanol dehydrogenase marks methylotrophs, but it is often not recorded as a trait.",
}


# the pre-registered (Addendum J) run as recorded in Addendum K, before any later change to the grounding logic:
# no positive universal rule; rule-out rules: 3 for oxygen use, 15 for the pathways, 0 for donors/acceptors
PREREG_RECORDED = {"da": 0, "ox": 3, "pw": 15}
J_P1_TRAITS = {"methanogenesis": "methanogenesis", "sulfate_reduction": "sulfate reduction",
               "nitrogen_fixation": "nitrogen fixation", "uses_O2": "oxygen use"}


def main():
    out = {"tasks": {}}
    out["n_external_novel"] = int(len(importlib.import_module("51_rulebook").load_external()[0]))
    for task in ("da", "ox", "pw"):
        f = OUT / f"universal_{task}_posthoc.csv"
        if not f.exists():
            out["tasks"][task] = {"status": "pending", "name": TN[task]}
            continue
        d = pd.read_csv(f)
        d = (d.sort_values("universal", ascending=False)                 # keep the universal variant (round 3)
             .drop_duplicates(["trait", "cls", "enzyme_signature"]))    # subunit variants of one enzyme rule
        lc = OUT / f"label_completeness_{task}.csv"
        pre = OUT / f"universal_{task}.csv"
        d["why_general"] = d.trait.map(WHY_GENERAL).fillna("")
        # a positive rule that requires the trait's own enzyme to be ABSENT is indirect: the generic "why" would
        # contradict it, so it gets an explicit note instead
        from gene_mechanisms import PATHWAY
        import re as _re
        indirect = []
        for r in d.itertuples():
            absent_path = [k for k, st in _re.findall(r"\((K\d{5})\) (present|absent)", r.rule)
                           if st == "absent" and k in PATHWAY.get(r.trait, set())]
            indirect.append(bool(r.cls == "has trait" and absent_path and not any(
                k in PATHWAY.get(r.trait, set()) for k, st in _re.findall(r"\((K\d{5})\) (present|absent)", r.rule)
                if st == "present")))
        d["indirect"] = indirect
        # rule-specific note: cytochrome bd is a poor marker alone, but with an aa3 oxidase it marks a branched chain
        both = d.rule.str.contains("K0042[56]") & d.rule.str.contains("K0227[456]") & (d.cls == "has trait")
        d.loc[both, "why_general"] = (d.loc[both, "why_general"] + " Cytochrome bd on its own is not a marker — many "
                                      "strict anaerobes keep it to scavenge traces of O₂ — but together with an aa3 "
                                      "oxidase it marks a branched respiratory chain: a low-affinity oxidase for "
                                      "plentiful O₂ and a high-affinity one for when O₂ runs low, typical of aerobes "
                                      "and facultative anaerobes. Pyruvate dehydrogenase adds the aerobic route into "
                                      "the citric-acid cycle.")
        d.loc[d.indirect, "why_general"] = ("Indirect rule: it predicts the trait from the absence of its canonical "
                                            "enzyme together with lifestyle markers. It picks organisms that carry out "
                                            "the physiology with other enzyme classes (here, anaerobes using "
                                            "hydrogenases other than HyaB/HybC). Statistically general, but not a "
                                            "direct gene-to-function link.")
        b = {"status": "done", "name": TN[task], "rules": rec(d.sort_values(["universal", "trait", "cls", "rank"],
                                                                             ascending=[False, True, True, True]))}
        if lc.exists():
            L = pd.read_csv(lc)
            L["why"] = L.trait.map(LIMITED_BY_LABELS).fillna("")
            b["label_completeness"] = rec(L)
        if pre.exists():   # the pre-registered run (Addendum J): its universal rules were all rule-out rules
            P = pd.read_csv(pre)
            b["preregistered_universal"] = rec(P[P.universal])
            b["preregistered_n_pos"] = int((P.universal & (P.cls == "has trait")).sum())
            b["preregistered_recorded_n_neg"] = PREREG_RECORDED[task]      # as recorded in Addendum K
            b["preregistered_rerun_n_neg"] = int((P.universal & (P.cls == "lacks trait")).sum())
        b["n_universal_pos"] = int(((d.universal) & (d.cls == "has trait")).sum())
        b["n_universal_neg"] = int(((d.universal) & (d.cls == "lacks trait")).sum())
        b["n_traits"] = int(d.trait.nunique())
        b["traits_with_universal"] = sorted(d[(d.universal) & (d.cls == "has trait")].trait.unique())
        out["tasks"][task] = b
    # J-P1 (pre-registered): >= 1 universal positive rule for each single-committed-marker trait
    have = {t for b in out["tasks"].values() if b.get("status") == "done" for t in b.get("traits_with_universal", [])}
    out["j_p1"] = dict(prediction="at least one universal rule for each of methanogenesis, sulfate reduction, nitrogen "
                                  "fixation and oxygen use",
                       preregistered_verdict="not supported (no positive universal rule at all; Addendum K)",
                       posthoc=f"{sum(t in have for t in J_P1_TRAITS)} of {len(J_P1_TRAITS)} traits have one post hoc: "
                               + ", ".join(f"{v} {'yes' if k in have else 'no'}" for k, v in J_P1_TRAITS.items()))
    (OUT / "universal_data.json").write_text(json.dumps(out, indent=1, default=float))
    print({t: (b.get("n_universal_pos"), b.get("n_universal_neg")) for t, b in out["tasks"].items()})


if __name__ == "__main__":
    main()
