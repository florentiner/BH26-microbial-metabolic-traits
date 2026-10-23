"""Gates G2 (label matrix) and G3 (frozen phylogeny-blocked splits).

Labels
  * 228 genomes flagged "Not in the electron donor / acceptor database" carry no trait
    information -> unlabelled (Y = NaN for every trait). They are never used for training
    or evaluation; they form the application set.
  * For the 358 database genomes a category flag of 0 is treated as a negative. This is
    a closed-world assumption (the database is taken as exhaustive for listed organisms)
    and is reported as a limitation.
  * Two "Not in database" columns are dropped; 14 traits remain (7 donor, 7 acceptor).
  * A trait is evaluable if it has >= 10 positives and >= 10 negatives.

Splits (5 outer folds, frozen, hashed)
  * genus level : groups = genus
  * family level: groups = family; a genome with missing family is its own genus-keyed
    group ("genus:<name>") so novel lineages are never lumped under NaN.
  Folds are assigned by a deterministic greedy multilabel-balanced group allocation.
  Every genome of a group lands in exactly one test fold.

Output: data/processed/labels.tsv, data/processed/splits.tsv, data/processed/splits.sha256,
        results/trait_counts.csv
"""
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
N_FOLDS = 5
SEED = 20260917
NOT_IN_DB = "Not in the electron donor / acceptor database"
SHORT = {
    "Sulfur oxidizing bacteria": "donor_sulfur_oxidation",
    "Ammonia oxidizing bacteria": "donor_ammonia_oxidation",
    "Nitrite oxidizing bacteria": "donor_nitrite_oxidation",
    "Hydrogen oxidizing bacteria": "donor_hydrogen_oxidation",
    "Carbon monoxide oxidizing bacteria (carboxydotrophs)": "donor_CO_oxidation",
    "Iron oxidizing bacteria": "donor_iron_oxidation",
    "Aerobes / aerobic respiration": "acceptor_oxygen",
    "Nitrate-reducing bacteria (denitrifiers / DNRA)": "acceptor_nitrate",
    "Sulfate-reducing bacteria (SRB)": "acceptor_sulfate",
    "Sulfur-reducing bacteria": "acceptor_sulfur",
    "Methanogens / Acetogens": "acceptor_CO2_methanogen_acetogen",
    "Iron-reducing bacteria (FeRB)": "acceptor_iron",
}


def trait_name(col):
    kind, cat = col.split(" categories: ")
    if cat == "Other":
        return f"{kind.lower()}_other"
    return SHORT[cat]


def greedy_group_folds(groups, Y, n_folds, seed):
    """Assign whole groups to folds, balancing genome counts and per-trait positive counts.

    Groups are processed from the one containing the rarest positives to the most common;
    each goes to the fold with the lowest normalised load. Ties are broken by a seeded RNG.
    """
    rng = np.random.default_rng(seed)
    uniq = pd.unique(groups)
    pos = Y.sum(0)
    weight = 1.0 / np.maximum(pos, 1)
    gidx = {g: np.where(groups == g)[0] for g in uniq}
    rarity = {g: (Y[gidx[g]] * weight).sum() for g in uniq}
    order = sorted(uniq, key=lambda g: (-rarity[g], -len(gidx[g]), rng.random()))
    load_n = np.zeros(n_folds)
    load_y = np.zeros((n_folds, Y.shape[1]))
    target_n = len(groups) / n_folds
    target_y = np.maximum(pos / n_folds, 1e-9)
    fold_of = {}
    for g in order:
        idx = gidx[g]
        cost = []
        for f in range(n_folds):
            ny = (load_y[f] + Y[idx].sum(0)) / target_y
            cost.append(((load_n[f] + len(idx)) / target_n) ** 2 + (ny ** 2 * (pos > 0)).sum() / max((pos > 0).sum(), 1))
        cost = np.array(cost) + rng.random(n_folds) * 1e-9
        f = int(np.argmin(cost))
        fold_of[g] = f
        load_n[f] += len(idx)
        load_y[f] += Y[idx].sum(0)
    return np.array([fold_of[g] for g in groups])


def main():
    df = pd.read_csv(ROOT / "ko_presence_586genomes.tsv", sep="\t", encoding="utf-8-sig")
    label_cols = [c for c in df.columns if " categories: " in c and not c.endswith(NOT_IN_DB)]
    unlabelled = df[f"Donor categories: {NOT_IN_DB}"].astype(bool)
    assert (unlabelled == df[f"Acceptor categories: {NOT_IN_DB}"].astype(bool)).all()
    assert (df.loc[unlabelled, label_cols].sum(axis=1) == 0).all()

    tax = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    labels = df[["genome_id", "organism_name"] + tax].copy()
    labels["labelled"] = ~unlabelled
    traits = []
    for c in label_cols:
        t = trait_name(c)
        traits.append(t)
        labels[t] = df[c].astype(float).where(~unlabelled)
    labels.to_csv(ROOT / "data/processed/labels.tsv", sep="\t", index=False)

    lab = labels[labels.labelled].reset_index(drop=True)
    Y = lab[traits].to_numpy()
    counts = pd.DataFrame({
        "trait": traits,
        "source_column": label_cols,
        "n_pos": Y.sum(0).astype(int),
        "n_neg": (Y == 0).sum(0),
        "n_unobserved": int(unlabelled.sum()),
    })
    counts["prevalence"] = counts.n_pos / (counts.n_pos + counts.n_neg)
    counts["evaluable"] = (counts.n_pos >= 10) & (counts.n_neg >= 10)
    counts["n_genera_with_pos"] = [lab.loc[lab[t] == 1, "genus"].nunique() for t in traits]
    counts["n_families_with_pos"] = [lab.loc[lab[t] == 1, "family"].nunique(dropna=False) for t in traits]
    (ROOT / "results").mkdir(exist_ok=True)
    counts.to_csv(ROOT / "results/trait_counts.csv", index=False)
    print(counts.to_string())

    genus_group = lab.genus.astype(str).to_numpy()
    fam_group = np.where(lab.family.isna(), "genus:" + lab.genus.astype(str), "family:" + lab.family.astype(str))
    splits = lab[["genome_id", "genus", "family"]].copy()
    splits["group_genus"] = genus_group
    splits["group_family"] = fam_group
    splits["fold_genus"] = greedy_group_folds(genus_group, Y, N_FOLDS, SEED)
    splits["fold_family"] = greedy_group_folds(fam_group, Y, N_FOLDS, SEED)
    for lvl in ("genus", "family"):
        per_group = splits.groupby(f"group_{lvl}")[f"fold_{lvl}"].nunique()
        assert (per_group == 1).all(), f"{lvl} group split across folds"
        print(lvl, "groups:", splits[f"group_{lvl}"].nunique())
        print(pd.DataFrame(Y, columns=traits).groupby(splits[f"fold_{lvl}"]).sum().T.to_string())
    out = ROOT / "data/processed/splits.tsv"
    splits.to_csv(out, sep="\t", index=False)
    digest = hashlib.sha256(out.read_bytes()).hexdigest()
    (ROOT / "data/processed/splits.sha256").write_text(f"{digest}  splits.tsv\n")
    print("splits sha256", digest)


if __name__ == "__main__":
    main()
