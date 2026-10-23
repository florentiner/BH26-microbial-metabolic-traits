"""Label-quality audit and annotation concordance (audit round 1, item M1/M11).

1. Marker-KO concordance: provided marker table vs KOfam recomputation (per KO and overall).
2. Marker contradictions: labelled negatives that carry a complete, trait-specific marker complex, and
   biologically implausible label combinations (chemolithotrophic ammonia/nitrite oxidisers without O2 as acceptor).
3. High-confidence out-of-fold false positives of the primary KO arms (genus hold-out), with marker evidence.
These are *candidate* label errors for expert review; nothing here changes the labels used for training.

Output: results/qc/marker_ko_concordance.csv, results/qc/label_contradictions.csv,
        results/candidate_label_errors.csv, results/qc/label_audit_summary.json
"""
import json

import numpy as np
import pandas as pd

from redox_common import RES, RULE_GROUPS, load_features, load_labelled

SPECIFIC = {  # trait -> marker groups; complete complex in a labelled negative is flagged
    "acceptor_CO2_methanogen_acetogen": ["mcrABG"],
    "donor_ammonia_oxidation": ["amoABC"],
    "acceptor_sulfate": ["dsrAB"],
    "acceptor_oxygen": ["coxABC", "ccoNOQP"],
}
# Only mcrABG is treated as strong evidence. The rest are ambiguous: amoABC KOs = pmoABC (methanotrophs);
# dsrAB also works in reverse in sulfur oxidisers; heme-copper oxidases occur in facultative anaerobes/microaerophiles.
AMBIGUOUS = {"donor_ammonia_oxidation", "acceptor_sulfate", "acceptor_oxygen"}


def main():
    qc = RES / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    lab, traits, Y = load_labelled()
    ids = pd.read_csv(RES.parent / "ko_presence_586genomes.tsv", sep="\t", encoding="utf-8-sig").genome_id
    prov = load_features("markers_provided", ids)
    rec = load_features("markers_recomputed", ids)[prov.columns]
    conc = pd.DataFrame({"n_provided": prov.sum(), "n_recomputed": rec.sum(),
                         "provided1_recomputed0": ((prov == 1) & (rec == 0)).sum(),
                         "provided0_recomputed1": ((prov == 0) & (rec == 1)).sum()})
    conc.index.name = "KO"
    conc.sort_values("provided1_recomputed0", ascending=False).to_csv(qc / "marker_ko_concordance.csv")

    M = load_features("markers_recomputed", lab.genome_id)
    frac = {(t, g): M[k].mean(axis=1).to_numpy() for t, gs in RULE_GROUPS.items() for g, k in gs.items()}
    rows = []
    for t, groups in SPECIFIC.items():
        j = traits.index(t)
        for g in groups:
            hit = (frac[(t, g)] == 1.0) & (Y[:, j] == 0)
            for i in np.where(hit)[0]:
                rows.append(dict(genome_id=lab.genome_id[i], organism=lab.organism_name[i], genus=lab.genus[i],
                                 trait=t, label=0, evidence=f"complete {g}",
                                 ambiguous=t in AMBIGUOUS))
    ox, ja = traits.index("acceptor_oxygen"), [traits.index("donor_ammonia_oxidation"),
                                               traits.index("donor_nitrite_oxidation")]
    for i in np.where((Y[:, ja].max(axis=1) == 1) & (Y[:, ox] == 0))[0]:
        rows.append(dict(genome_id=lab.genome_id[i], organism=lab.organism_name[i], genus=lab.genus[i],
                         trait="acceptor_oxygen", label=0,
                         evidence="labelled ammonia/nitrite oxidiser (obligate aerobes) without O2 acceptor",
                         ambiguous=False))
    contra = pd.DataFrame(rows)
    contra.to_csv(qc / "label_contradictions.csv", index=False)

    oof = pd.concat([pd.read_parquet(RES / "oof" / f"{a}.parquet") for a in ("B1r", "B2")])
    oof = oof[(oof.level == "genus") & (oof.y == 0) & (oof.prob >= oof.threshold)]
    fp = (oof.groupby(["genome_id", "trait"]).agg(n_arms=("arm", "nunique"), mean_prob=("prob", "mean"))
          .reset_index())
    fp = fp[fp.n_arms == 2].merge(lab[["genome_id", "organism_name", "genus", "family"]], on="genome_id")
    gi = {g: i for i, g in enumerate(lab.genome_id)}
    fp["marker_groups_complete"] = [
        ",".join(g for g in RULE_GROUPS[t] if frac[(t, g)][gi[a]] == 1.0) for a, t in zip(fp.genome_id, fp.trait)]
    fp["block_all_zero"] = [
        bool(Y[gi[a], [j for j, x in enumerate(traits) if x.split("_")[0] == t.split("_")[0]]].sum() == 0)
        for a, t in zip(fp.genome_id, fp.trait)]
    fp["in_marker_contradiction_list"] = [
        ((contra.genome_id == a) & (contra.trait == t)).any() if len(contra) else False
        for a, t in zip(fp.genome_id, fp.trait)]
    fp = fp.sort_values("mean_prob", ascending=False)
    fp.to_csv(RES / "candidate_label_errors.csv", index=False)

    summary = {
        "marker_cells_total": int(prov.size), "marker_cells_agree": int((prov.values == rec.values).sum()),
        "provided1_recomputed0": int(conc.provided1_recomputed0.sum()),
        "provided0_recomputed1": int(conc.provided0_recomputed1.sum()),
        "labelled_genomes_zero_donor_block": int((Y[:, [j for j, t in enumerate(traits) if t.startswith("donor_")]]
                                                  .sum(1) == 0).sum()),
        "labelled_genomes_zero_acceptor_block": int((Y[:, [j for j, t in enumerate(traits)
                                                          if t.startswith("acceptor_")]].sum(1) == 0).sum()),
        "marker_contradictions_by_trait": contra.groupby(["trait", "ambiguous"]).size().rename("n").reset_index()
        .to_dict("records") if len(contra) else [],
        "oof_false_positives_both_B1r_B2_genus": int(len(fp)),
        "of_which_with_complete_marker_group": int((fp.marker_groups_complete != "").sum()),
        "of_which_block_all_zero": int(fp.block_all_zero.sum()),
    }
    (qc / "label_audit_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
