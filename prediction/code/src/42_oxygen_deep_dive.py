"""Why the aerobic-respiration markers track the curated O2 label but not the recorded oxygen requirement.

The same marker group scores MCC 0.39 against the curated `acceptor_oxygen` label and ~0.0 against Madin's oxygen
requirement (`metabolism`). Exactly one of three things must be true, and each is testable:
  1. the two phenotype sources disagree with each other  -> measured here as Cohen's kappa,
  2. the marker set is too broad (cytochrome bd is carried by anaerobes for O2 detoxification, not respiration)
     -> measured by scoring each oxidase KO separately against oxygen requirement,
  3. the association is real but diluted -> measured by restricting to genomes where both sources agree.

Output: results/jaccard/oxygen_deep_dive.csv, oxygen_label_agreement.csv
usage: python src/42_oxygen_deep_dive.py
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/jaccard"
SEED = 20260918

spec = importlib.util.spec_from_file_location("jc", ROOT / "src/40_jaccard_correlations.py")
jc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(jc)

# The 13 KOs of the "Aerobic respiration" group, split by what the enzyme actually implies. Cytochrome bd (cydAB)
# is a high-affinity oxidase that many strict anaerobes keep for O2 scavenging; aa3-type cytochrome c oxidase
# (coxABC) is the canonical aerobic respiratory chain; cbb3 (ccoNOQP) is the microaerobic one.
OXIDASE_CLASS = {
    "K02274": "aa3 cytochrome c oxidase (coxA)", "K02275": "aa3 cytochrome c oxidase (coxB)",
    "K02276": "aa3 cytochrome c oxidase (coxC)",
    "K00404": "cbb3 oxidase (ccoN)", "K00405": "cbb3 oxidase (ccoO)", "K00406": "cbb3 oxidase (ccoP)",
    "K00425": "cytochrome bd (cydA)", "K00426": "cytochrome bd (cydB)",
}


def kappa(a, b):
    a, b = np.asarray(a, dtype=int), np.asarray(b, dtype=int)
    po = float((a == b).mean())
    pe = float(a.mean() * b.mean() + (1 - a.mean()) * (1 - b.mean()))
    return (po - pe) / (1 - pe) if pe < 1 else np.nan, po


def main():
    rng = np.random.default_rng(SEED)
    X, groups, role, ranks, info = jc.load_html(ROOT / "jaccard_similarity.html")
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    lab = lab[lab.labelled] if "labelled" in lab.columns else lab
    lab = lab.set_index("genome_id")

    mad = pd.read_csv(ROOT / "data/raw/madin/condensed_species_NCBI.csv", low_memory=False)
    sp = info["species"].str.replace("Candidatus ", "", regex=False).str.split().str[:2].str.join(" ")
    ox = sp.map(mad.drop_duplicates("species").set_index("species").metabolism) \
           .map(lambda v: jc.OXYGEN.get(str(v).strip().lower(), np.nan))

    # ---- 1. do the two phenotype sources agree?
    both = [g for g in X.index if g in lab.index and pd.notna(ox.get(g))]
    curated = lab.loc[both, "acceptor_oxygen"].to_numpy(dtype=int)
    aerobe = (ox[both].isin(["aerobe", "facultative", "microaerophile"])).astype(int).to_numpy()
    strict_aerobe = (ox[both] == "aerobe").astype(int).to_numpy()
    rows = []
    for name, y in (("O2-capable (aerobe/facultative/microaerophile)", aerobe), ("strict aerobe only", strict_aerobe)):
        k, po = kappa(curated, y)
        rows.append(dict(comparison=f"curated acceptor_oxygen vs Madin {name}", n=len(both),
                         n_curated_pos=int(curated.sum()), n_madin_pos=int(y.sum()),
                         raw_agreement=po, cohens_kappa=k,
                         curated_pos_madin_neg=int(((curated == 1) & (y == 0)).sum()),
                         curated_neg_madin_pos=int(((curated == 0) & (y == 1)).sum())))
    agr = pd.DataFrame(rows)
    agr.to_csv(OUT / "oxygen_label_agreement.csv", index=False, float_format="%.4f")
    print("=== Do the two oxygen phenotype sources agree? ===")
    print(agr.round(3).to_string(index=False))

    # ---- 2. per-oxidase association with oxygen requirement, and 3. the agreeing subset
    genera_b = ranks["genus"].reindex(both).fillna("?").to_numpy()
    rows = []
    variants = {"Aerobic respiration (all 13 KOs)": groups["Aerobic respiration"],
                "aa3 cytochrome c oxidase only": ["K02274", "K02275", "K02276"],
                "cbb3 oxidase only": ["K00404", "K00405", "K00406"],
                "cytochrome bd only": ["K00425", "K00426"],
                "aa3 or cbb3 (true respiratory oxidases)": ["K02274", "K02275", "K02276",
                                                            "K00404", "K00405", "K00406"]}
    for vname, kos in variants.items():
        kos = [k for k in kos if k in X.columns]
        if not kos:
            continue
        x = X.loc[both, kos].max(axis=1).to_numpy()
        for pname, y in (("Madin O2-capable", aerobe), ("Madin strict aerobe", strict_aerobe),
                         ("curated acceptor_oxygen", curated)):
            r = jc.associate(x, y, genera_b, rng, marker_set=vname, n_kos=len(kos), phenotype=pname)
            rows.append(r)
    # the agreeing subset: genomes where both sources say the same thing
    agree_mask = curated == aerobe
    if agree_mask.sum() >= 50:
        idx = [g for g, m in zip(both, agree_mask) if m]
        x = X.loc[idx, groups["Aerobic respiration"]].max(axis=1).to_numpy()
        y = lab.loc[idx, "acceptor_oxygen"].to_numpy(dtype=int)
        rows.append(jc.associate(x, y, ranks["genus"].reindex(idx).fillna("?").to_numpy(), rng,
                                 marker_set="Aerobic respiration (all 13 KOs)",
                                 n_kos=len(groups["Aerobic respiration"]),
                                 phenotype=f"both sources agree (n={len(idx)})"))
    dd = pd.DataFrame(rows)
    dd["q_bh_permutation"] = jc.bh(dd.p_genus_permutation)
    cols = ["marker_set", "n_kos", "phenotype", "n", "n_marker_present", "n_phenotype_pos", "precision", "recall",
            "mcc_per_genome", "mcc_genus_level", "odds_ratio", "p_genus_permutation", "q_bh_permutation"]
    dd[cols].to_csv(OUT / "oxygen_deep_dive.csv", index=False, float_format="%.4f")
    print("\n=== Which oxidase actually carries oxygen information? ===")
    print(dd[cols].round(3).to_string(index=False))

    # prevalence of each oxidase among aerobes and anaerobes - the clearest picture of the dilution problem
    rows = []
    for ko, label in OXIDASE_CLASS.items():
        if ko not in X.columns:
            continue
        x = X.loc[both, ko].to_numpy()
        rows.append(dict(KO=ko, enzyme=label,
                         prevalence_in_O2_capable=float(x[aerobe == 1].mean()),
                         prevalence_in_anaerobes=float(x[aerobe == 0].mean()),
                         difference=float(x[aerobe == 1].mean() - x[aerobe == 0].mean())))
    pv = pd.DataFrame(rows).sort_values("difference", ascending=False)
    pv.to_csv(OUT / "oxidase_prevalence.csv", index=False, float_format="%.3f")
    print("\n=== Prevalence of each oxidase in O2-capable vs anaerobic organisms ===")
    print(pv.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
