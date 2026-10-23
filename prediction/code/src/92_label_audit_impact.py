"""Impact of the October 2026 audit of the phenotype tables on the donor/acceptor training labels.

Inputs: the pre-audit phenotype table (hack_res, 4 Oct) and the audited table (hack_res_fedotmas/final). Donor and
acceptor compounds are mapped to the 14 categories with a fixed keyword map; organisms are matched to the 358
training genomes by binomial. Reports (1) category calls that the audit flips for training species, and (2) how
each derived matrix agrees with data/processed/labels.tsv (the flags the models were trained on).

Note: the training flags came from the curated 586-genome table, not directly from these phenotype tables, so this
measures consistency, not a direct relabelling.

usage: python src/92_label_audit_impact.py
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/additional"
OUT.mkdir(parents=True, exist_ok=True)
BEFORE = Path("~/work/hack_res/phenotype_data_GRADED - phenotype_data_GRADED.csv.csv")
AFTER = Path("~/work/hack_res_fedotmas/final/phenotype_data_GRADED_VALIDATED.csv")
AFTER_OA = Path("~/work/hack_res_fedotmas/final/phenotype_data_GRADED_OPEN_ACCESS.csv")

DONOR = [("donor_hydrogen_oxidation", r"\bH2\b|hydrogen"),
         ("donor_ammonia_oxidation", r"NH3|NH4|ammoni"),
         ("donor_nitrite_oxidation", r"NO2-|nitrite"),
         ("donor_CO_oxidation", r"\bCO\b|carbon monoxide"),
         ("donor_iron_oxidation", r"Fe\(II\)|ferrous|pyrite"),
         ("donor_sulfur_oxidation", r"H2S|HS-|sulfide|S0|elemental sulfur|thiosulfate|S2O3|tetrathionate|S4O6|"
                                    r"trithionate|S3O6|sulfite|SO3|sulfur compounds")]
ACCEPTOR = [("acceptor_oxygen", r"\bO2\b|oxygen"),
            ("acceptor_nitrate", r"NO3|nitrate"),
            ("acceptor_sulfate", r"SO4|sulfate"),
            ("acceptor_sulfur", r"S0|elemental sulfur|thiosulfate|S2O3|sulfite|SO3|tetrathionate|polysulfide"),
            ("acceptor_CO2_methanogen_acetogen", r"\bCO2\b|carbon dioxide|methanogen|acetogen"),
            ("acceptor_iron", r"Fe\(III\)|ferric")]
CATS = [c for c, _ in DONOR] + [c for c, _ in ACCEPTOR]


def entries(cell):
    if not isinstance(cell, str) or not cell.strip() or cell.strip() in ("-", "nan"):
        return []
    return [e.strip() for e in cell.split(";") if e.strip()]


def categorise(cell, rules, other):
    cats = set()
    for e in entries(cell):
        e2 = re.sub(r"thiosulfate", "thiosulfate", e)
        hit = [c for c, rx in rules if re.search(rx, e2, flags=re.I)]
        # "sulfate" must not match thiosulfate
        if "acceptor_sulfate" in hit and re.search(r"thiosulfate", e2, re.I) and not re.search(r"\bSO4", e2):
            hit.remove("acceptor_sulfate")
        cats |= set(hit) if hit else {other}
    return cats


def matrix(path):
    p = pd.read_csv(path, low_memory=False)
    p["binomial"] = p.organism.astype(str).str.replace(r"\s+sp\. nov\.?", "", regex=True) \
        .str.replace(r"\(.*?\)", "", regex=True).str.split().str[:2].str.join(" ")
    rows = []
    for b, d in p.groupby("binomial"):
        s = set()
        for _, r in d.iterrows():
            s |= categorise(r.get("electron_donors"), DONOR, "donor_other")
            s |= categorise(r.get("electron_acceptors"), ACCEPTOR, "acceptor_other")
        rows.append({"binomial": b, **{c: int(c in s) for c in CATS + ["donor_other", "acceptor_other"]},
                     "any_info": int(bool(s))})
    return pd.DataFrame(rows).set_index("binomial")


def main():
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    lab = lab[lab.labelled == True].copy()  # noqa: E712
    lab["binomial"] = lab.species.fillna(lab.organism_name).astype(str).str.split().str[:2].str.join(" ")
    B, A, O = matrix(BEFORE), matrix(AFTER), matrix(AFTER_OA)
    common = sorted(set(lab.binomial) & set(B.index) & set(A.index))
    L = lab.drop_duplicates("binomial").set_index("binomial").loc[common, CATS].astype(int)
    rows, flips = [], []
    for c in CATS:
        b, a, o, y = B.loc[common, c], A.loc[common, c], O.loc[common, c], L[c]
        info = A.loc[common, "any_info"] == 1
        rows.append(dict(category=c, species_matched=len(common), species_with_info=int(info.sum()),
                         flipped_by_audit=int((b != a).sum()),
                         added_by_audit=int(((b == 0) & (a == 1)).sum()),
                         removed_by_audit=int(((b == 1) & (a == 0)).sum()),
                         train_label_pos=int(y.sum()),
                         agree_train_vs_before=float((b[info] == y[info]).mean()),
                         agree_train_vs_after=float((a[info] == y[info]).mean()),
                         agree_train_vs_open_access=float((o[info] == y[info]).mean()),
                         train_pos_but_after_neg=int(((y == 1) & (a == 0) & info).sum()),
                         train_neg_but_after_pos=int(((y == 0) & (a == 1) & info).sum())))
        for s in np.array(common)[(b != a).to_numpy()]:
            flips.append(dict(species=s, category=c, before=int(B.at[s, c]), after=int(A.at[s, c]),
                              training_label=int(L.at[s, c])))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "label_audit_impact_per_category.csv", index=False)
    F = pd.DataFrame(flips)
    F.to_csv(OUT / "label_audit_flips.csv", index=False)
    pd.set_option("display.width", 220)
    print(R.round(3).to_string(index=False))
    print("\nflips:", len(F), "| flips where the audit now disagrees with the training label:",
          int((F.after != F.training_label).sum()) if len(F) else 0)
    print(F.to_string(index=False) if len(F) else "")


if __name__ == "__main__" and len(__import__("sys").argv) == 1:
    main()


def rescore():
    """Re-score saved out-of-fold predictions with the audit's disagreeing cells applied to the labels."""
    from sklearn.metrics import average_precision_score
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    lab = lab[lab.labelled == True].copy()  # noqa: E712
    lab["binomial"] = lab.species.fillna(lab.organism_name).astype(str).str.split().str[:2].str.join(" ")
    F = pd.read_csv(OUT / "label_audit_flips.csv")
    F = F[F.after != F.training_label]
    change = {(g, r.category): r.after for _, r in F.iterrows() for g in lab.genome_id[lab.binomial == r.species]}
    rows = []
    for arm in ["B2", "B2XL", "B2XLE_large", "null_taxNN"]:
        f = ROOT / "results/oof" / ("nulls.parquet" if arm.startswith("null") else f"{arm}.parquet")
        d = pd.read_parquet(f)
        d = d[d.arm == arm]
        for level, dl in d.groupby("level"):
            old, new = [], []
            for t, dt in dl.groupby("trait"):
                y = dt.y.to_numpy()
                y2 = np.array([change.get((g, t), yy) for g, yy in zip(dt.genome_id, y)])
                old.append(average_precision_score(y, dt.prob))
                new.append(average_precision_score(y2, dt.prob))
            rows.append(dict(arm=arm, level=level, macro_auprc_original=np.mean(old),
                             macro_auprc_audited=np.mean(new), delta=np.mean(new) - np.mean(old)))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "label_audit_rescored_cv.csv", index=False)
    print("cells changed:", len(change))
    print(R.round(4).to_string(index=False))


if __name__ == "__main__" and len(__import__("sys").argv) > 1 and __import__("sys").argv[1] == "rescore":
    rescore()
