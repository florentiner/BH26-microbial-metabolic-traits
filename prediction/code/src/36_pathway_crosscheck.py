"""Cross-source validation of the pathway model against independently curated donor/acceptor labels.

The pathway labels come from Madin et al. 2020; the donor/acceptor labels come from the provided 586-genome
curation. They are different curations by different people from different literature. Where a pathway trait and a
donor/acceptor category name the same physiology, the *out-of-fold* pathway prediction (its genus was held out of
pathway training) can be scored against the donor/acceptor curation. Agreement there cannot come from label
leakage, since the two label sets were produced independently.

Genomes are matched by assembly accession first, then by species name (both sets are one genome per species).
Output: results/pathways/crosscheck_donor_acceptor.csv
usage: python src/36_pathway_crosscheck.py [--arm PB2]
"""
import argparse

import numpy as np
import pandas as pd

from pathway_common import RES, load_pathways
from redox_common import PROC, safe_ap

# pathway trait -> donor/acceptor category naming the same physiology (fixed before scoring)
PAIRS = {
    "methanogenesis": "acceptor_CO2_methanogen_acetogen",
    "acetogenesis": "acceptor_CO2_methanogen_acetogen",
    "sulfate_reduction": "acceptor_sulfate",
    "sulfur_reduction": "acceptor_sulfur",
    "nitrate_reduction": "acceptor_nitrate",
    "denitrification": "acceptor_nitrate",
    "dnra": "acceptor_nitrate",
    "iron_reduction": "acceptor_iron",
    "sulfur_oxidation": "donor_sulfur_oxidation",
    "hydrogen_oxidation": "donor_hydrogen_oxidation",
    "nitrification": "donor_ammonia_oxidation",
    "iron_oxidation": "donor_iron_oxidation",
    "aerobic_chemoheterotrophy": "acceptor_oxygen",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="PB2")
    args = ap.parse_args()
    lab, traits, Y = load_pathways()
    oof = pd.read_parquet(RES / f"pathways/oof/{args.arm}.parquet")

    da = pd.read_csv(PROC / "labels.tsv", sep="\t")
    da = da[da.labelled] if "labelled" in da.columns else da
    da_traits = [c for c in da.columns if c.startswith(("donor_", "acceptor_"))]
    by_acc = da.set_index("genome_id")
    # the donor/acceptor table carries organism_name (strain-level); the species name is its first two words
    sp = da.organism_name.fillna("").str.replace("Candidatus ", "", regex=False).str.split().str[:2].str.join(" ")
    by_species = da.assign(sp=sp).drop_duplicates("sp").set_index("sp")

    rows = []
    for level in ("genus", "family"):
        d = oof[oof.level == level]
        for pw, cat in PAIRS.items():
            if pw not in traits or cat not in da_traits:
                continue
            dt = d[d.trait == pw].merge(lab[["genome_id", "species"]], on="genome_id", how="left")
            y_ext, p, y_own = [], [], []
            for r in dt.itertuples():
                if r.genome_id in by_acc.index:
                    ext = by_acc.at[r.genome_id, cat]
                elif isinstance(r.species, str) and r.species.strip() in by_species.index:
                    ext = by_species.at[r.species.strip(), cat]
                else:
                    continue
                y_ext.append(float(ext))
                p.append(r.prob)
                y_own.append(r.y)
            if len(y_ext) < 20 or len(set(y_ext)) < 2:
                continue
            y_ext, p, y_own = np.array(y_ext), np.array(p), np.array(y_own)
            rows.append(dict(split_level=level, pathway_trait=pw, donor_acceptor_trait=cat, n_matched=len(y_ext),
                             n_pos_external=int(y_ext.sum()), prevalence_external=float(y_ext.mean()),
                             auprc_vs_external=safe_ap(y_ext, p), auprc_lift=safe_ap(y_ext, p) / y_ext.mean(),
                             auprc_vs_own_labels=safe_ap(y_own, p),
                             label_agreement=float((y_ext == y_own).mean())))
    out = pd.DataFrame(rows)
    (RES / "pathways").mkdir(parents=True, exist_ok=True)
    out.to_csv(RES / "pathways/crosscheck_donor_acceptor.csv", index=False, float_format="%.4f")
    print(out.round(3).to_string(index=False))
    if len(out):
        g = out[out.split_level == "genus"]
        print(f"\nmatched pairs: {len(g)}; mean AUPRC vs independent curation {g.auprc_vs_external.mean():.3f} "
              f"(lift {g.auprc_lift.mean():.2f}x); mean label agreement {g.label_agreement.mean():.3f}")


if __name__ == "__main__":
    main()
