"""Addendum A (X1 comparator): score METABOLIC-G v4 function presence against the 255 blind growth tests.
Reads data/additional/metabolic/full_out/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet2.tsv
(function presence per genome) and applies the pre-registered mapping. Writes results/additional/x1_metabolic_calls.csv.

usage: python src/95_metabolic_score.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
W2 = ROOT / "data/additional/metabolic/full_out/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet2.tsv"
OUT = ROOT / "results/additional"
MAP = {
    "donor_sulfur_oxidation": [("Sulfur cycling", "Sulfide oxidation"), ("Sulfur cycling", "Sulfur oxidation"),
                               ("Sulfur cycling", "Thiosulfate oxidation")],
    "donor_ammonia_oxidation": [("Nitrogen cycling", "Ammonia oxidation")],
    "donor_nitrite_oxidation": [("Nitrogen cycling", "Nitrite oxidation")],
    "donor_hydrogen_oxidation": [("Hydrogenases", "Ni-Fe Hydrogenase")],
    "donor_CO_oxidation": [("C1 metabolism", "Aerobic CO oxidation")],
    "donor_iron_oxidation": [("Iron cycling", "Iron oxidation")],
    "acceptor_nitrate": [("Nitrogen cycling", "Nitrate reduction")],
    "acceptor_sulfate": [("Sulfur cycling", "Sulfate reduction")],
    "acceptor_sulfur": [("Sulfur cycling", "Sulfur reduction"), ("Sulfur cycling", "Sulfite reduction")],
    "acceptor_CO2_methanogen_acetogen": [("Methane metabolism", "Methane production"),
                                         ("Carbon fixation", "Wood Ljungdahl pathway")],
    "acceptor_iron": [("Iron cycling", "Iron reduction"), ("Iron cycling", "Metal (Iron/Manganese) reduction")],
}
OXY_KEYS = ["aa3 type", "bo type", "caa3-type", "cbb3-type"]  # bd type excluded (pre-registered)


def function_table(w2_path):
    """Function presence per genome: worksheet 2 (function summary), plus functions that METABOLIC reports only per
    HMM in worksheet 1 (its iron-cycling HMMs), present when any HMM of that function is present (deviation D8)."""
    w = pd.read_csv(w2_path, sep="\t")
    pres = [c for c in w.columns if c.endswith("Function presence")]
    w1 = pd.read_csv(Path(str(w2_path).replace("worksheet2", "worksheet1")), sep="\t")
    hp = [c for c in w1.columns if c.endswith("Hmm presence")]
    extra = w1[~w1.Function.isin(w.Function)].copy()
    extra[hp] = extra[hp].eq("Present")
    agg = extra.groupby(["Category", "Function"])[hp].any().reset_index()
    agg.columns = ["Category", "Function"] + [c.replace("Hmm presence", "Function presence") for c in hp]
    agg[pres] = agg[pres].map(lambda v: "Present" if v else "Absent")
    return pd.concat([w[["Category", "Function"] + pres], agg[["Category", "Function"] + pres]], ignore_index=True)


def main():
    w = function_table(W2)
    pres = [c for c in w.columns if c.endswith("Function presence")]
    w[pres] = w[pres].eq("Present")
    genomes = {c: c.replace(" Function presence", "") for c in pres}
    rows = []
    for col, g in genomes.items():
        res = {}
        for cat, pairs in MAP.items():
            m = pd.Series(False, index=w.index)
            for cat_name, fun in pairs:
                m |= w.Function.str.strip() == fun  # category labels differ between sheets; names are unique
            if m.sum() == 0:
                raise SystemExit(f"mapping rows not found for {cat}: {pairs}")
            res[cat] = bool(w.loc[m, col].any())
        ox = w.Category.str.startswith("Oxygen metabolism") & w.Function.str.contains("|".join(OXY_KEYS))
        res["acceptor_oxygen"] = bool(w.loc[ox, col].any())
        rows.append({"accession": g, **res})
    P = pd.DataFrame(rows).set_index("accession")
    calls = pd.read_csv(OUT / "x1_calls.csv")
    calls["call_METABOLIC"] = [int(P.at[a, c]) if (a in P.index and c in P.columns) else 0
                               for a, c in zip(calls.accession, calls.category)]
    calls[["accession", "category", "call_METABOLIC"]].to_csv(OUT / "x1_metabolic_calls.csv", index=False)
    print(P.mean().round(2).to_string())


if __name__ == "__main__":
    main()
