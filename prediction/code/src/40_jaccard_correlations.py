"""Correlation analysis for jaccard_similarity.html: do the marker-KO functional groups track real physiology?

The HTML carries, for 586 genomes, presence/absence of 84 marker KOs grouped into 15 functional groups
(7 electron-donor, 8 electron-acceptor), plus NCBI lineage. It carries no phenotype, so phenotypes are joined in:
  * curated electron donor / acceptor categories  (data/processed/labels.tsv, the project's own labels)
  * oxygen requirement and habitat                (Madin et al. 2020: `metabolism`, `isolation_source`)

For every (functional group, phenotype) pair the script reports the 2x2 table, MCC / phi, odds ratio with a
Haldane-Anscombe correction, Fisher's exact p, and BH q within each family of tests.

**Pseudo-replication is the main threat to these correlations.** 586 genomes are not 586 independent observations:
several genomes of one genus that share a marker and a label are close to one observation. Every
association is therefore reported three ways:
  1. per genome (naive, comparable to what a spreadsheet would give),
  2. per genus (one randomly drawn genome per genus, median over 200 draws) - removes within-genus replication,
  3. a genus-permutation test: labels are shuffled between whole genera, so the null keeps the phylogenetic
     structure and only breaks the group-phenotype link. This p-value is the honest one.

usage: python src/40_jaccard_correlations.py [--html jaccard_similarity.html]
Outputs (results/jaccard/):
  ko_presence.tsv, group_presence.tsv   the extracted matrices
  group_vs_trait.csv                    functional group vs curated donor/acceptor label
  group_vs_oxygen.csv                   functional group vs oxygen requirement
  group_vs_habitat.csv                  functional group vs habitat category
  ko_vs_trait.csv                       single-KO resolution (which KO carries a group's signal)
  summary.json
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/jaccard"
SEED = 20260918
N_PERM = 2000
N_GENUS_DRAWS = 200

# functional group -> the curated donor/acceptor category that names the same physiology (fixed before scoring)
GROUP_TRAIT = {
    "Sulfur oxidizing (SOB)": "donor_sulfur_oxidation",
    "Ammonia oxidizing (AOB/AOA)": "donor_ammonia_oxidation",
    "Nitrite oxidizing (NOB)": "donor_nitrite_oxidation",
    "Hydrogen oxidizing": "donor_hydrogen_oxidation",
    "CO oxidizing (carboxydotroph)": "donor_CO_oxidation",
    "Formate oxidizing": "donor_other",
    "Methane oxidizing (methanotroph)": "donor_other",
    "Aerobic respiration": "acceptor_oxygen",
    "Nitrate-reducing (denitrifiers)": "acceptor_nitrate",
    "Nitrate-reducing (DNRA)": "acceptor_nitrate",
    "Anammox": "acceptor_other",
    "Sulfate-reducing (SRB)": "acceptor_sulfate",
    "Sulfur-reducing": "acceptor_sulfur",
    "Methanogens (hydrogenotrophic)": "acceptor_CO2_methanogen_acetogen",
    "Acetogens (homoacetogens)": "acceptor_CO2_methanogen_acetogen",
}

# habitat category -> prefixes of Madin `isolation_source`, which is a controlled vocabulary of 72 tokens
# ("water_marine", "sediment_marine_hydrothermal", "host_animal_endotherm_intestinal", ...), not free text.
HABITAT = {
    "marine": ["water_marine", "sediment_marine"],
    "freshwater": ["water_fresh", "sediment_fresh", "water_groundwater"],
    "soil_sediment": ["soil", "sediment", "compost", "rock"],
    "host_associated": ["host"],
    "engineered": ["wastewater", "sludge", "bioreactor", "built_environment", "petroleum", "food", "mine"],
    "thermal_hypersaline": ["water_hotspring", "water_hypersaline", "sediment_hypersaline",
                            "sediment_marine_hydrothermal"],
    "anoxic_proxy": ["sediment", "sludge", "bioreactor", "host_animal_endotherm_intestinal",
                     "host_animal_endotherm_rumen"],
}
OXYGEN = {  # Madin `metabolism` -> coarse oxygen requirement
    "aerobic": "aerobe", "obligate aerobic": "aerobe", "microaerophilic": "microaerophile",
    "facultative": "facultative", "anaerobic": "anaerobe", "obligate anaerobic": "anaerobe",
    "strictly anaerobic": "anaerobe",
}


def load_html(path):
    html = Path(path).read_text(encoding="utf-8")
    d = json.loads(re.search(r'<script id="jaccard-data" type="application/json">(.*?)</script>', html, re.S).group(1))
    kos = d["koColumns"]
    X = pd.DataFrame(0, index=[g["id"] for g in d["genomes"]], columns=kos, dtype=np.int8)
    for g in d["genomes"]:
        X.iloc[X.index.get_loc(g["id"]), g["present"]] = 1
    X.index.name = "genome_id"
    groups = {g["functionalGroup"]: [kos[i] for i in g["kos"]] for g in d["groups"]}
    role = {g["functionalGroup"]: g["role"] for g in d["groups"]}
    ranks = pd.DataFrame(d["genomeRanks"], index=d["rankNames"]).T
    ranks.index.name = "genome_id"
    info = pd.DataFrame(d["genomeInfo"], index=d["genomeInfoLabels"]).T
    info.index.name = "genome_id"
    # genomeInfo / genomeRanks are stored in a different order from `genomes` (584 of 586 positions differ). Every
    # frame is reindexed to the marker matrix so that positional masks can never pair a genome's markers with
    # another genome's phenotype.
    ranks = ranks.reindex(X.index)
    info = info.reindex(X.index)
    assert info.index.equals(X.index) and ranks.index.equals(X.index)
    return X, groups, role, ranks, info


def mcc_phi(a, b):
    """MCC of two binary vectors (identical to the phi coefficient for 2x2)."""
    a, b = np.asarray(a, dtype=int), np.asarray(b, dtype=int)
    tp = int(((a == 1) & (b == 1)).sum())
    tn = int(((a == 0) & (b == 0)).sum())
    fp = int(((a == 1) & (b == 0)).sum())
    fn = int(((a == 0) & (b == 1)).sum())
    den = np.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    return (tp * tn - fp * fn) / den if den > 0 else np.nan, (tp, fp, fn, tn)


def genus_permutation_p(x, y, genera, rng, n=N_PERM):
    """p-value under a null that shuffles the phenotype between whole genera.

    Keeps the phylogenetic structure of both variables and breaks only their association, so a correlation that is
    merely 'both traits run in the same families' does not come out significant.
    """
    obs, _ = mcc_phi(x, y)
    if np.isnan(obs):
        return np.nan, obs
    gi = pd.Series(genera)
    gy = pd.Series(y).groupby(gi).mean()          # one phenotype value per genus
    order = gy.index.to_numpy()
    vals = gy.to_numpy()
    count = 0
    for _ in range(n):
        perm = dict(zip(order, rng.permutation(vals)))
        y_perm = (gi.map(perm).to_numpy() >= 0.5).astype(int)
        v, _ = mcc_phi(x, y_perm)
        if not np.isnan(v) and abs(v) >= abs(obs) - 1e-12:
            count += 1
    return (count + 1) / (n + 1), obs


def genus_level_mcc(x, y, genera, rng, draws=N_GENUS_DRAWS):
    """Median MCC over repeated draws of one genome per genus (removes within-genus replication)."""
    df = pd.DataFrame({"x": x, "y": y, "g": genera})
    out = []
    for _ in range(draws):
        s = df.groupby("g", sort=False).sample(1, random_state=int(rng.integers(1 << 31)))
        v, _ = mcc_phi(s.x.to_numpy(), s.y.to_numpy())
        out.append(v)
    out = np.array(out, dtype=float)
    return np.nanmedian(out), np.nanpercentile(out, 2.5), np.nanpercentile(out, 97.5)


def bh(p):
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if not ok.any():
        return q
    pv = p[ok]
    o = np.argsort(pv)
    r = pv[o] * len(pv) / np.arange(1, len(pv) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty_like(pv)
    out[o] = np.minimum(r, 1)
    q[ok] = out
    return q


def associate(x, y, genera, rng, **meta):
    phi, (tp, fp, fn, tn) = mcc_phi(x, y)
    try:
        p_fisher = fisher_exact([[tp, fp], [fn, tn]])[1]
    except ValueError:
        p_fisher = np.nan
    or_hald = ((tp + 0.5) * (tn + 0.5)) / ((fp + 0.5) * (fn + 0.5))
    p_perm, _ = genus_permutation_p(x, y, genera, rng)
    gm, glo, ghi = genus_level_mcc(x, y, genera, rng)
    return dict(**meta, n=len(x), n_marker_present=int(np.sum(x)), n_phenotype_pos=int(np.sum(y)),
                tp=tp, fp=fp, fn=fn, tn=tn,
                precision=tp / (tp + fp) if tp + fp else np.nan,
                recall=tp / (tp + fn) if tp + fn else np.nan,
                f1=2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else np.nan,
                mcc_per_genome=phi, odds_ratio=or_hald, p_fisher_naive=p_fisher,
                mcc_genus_level=gm, mcc_genus_ci_low=glo, mcc_genus_ci_high=ghi,
                p_genus_permutation=p_perm)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--html", default=str(ROOT / "jaccard_similarity.html"))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    X, groups, role, ranks, info = load_html(args.html)
    G = pd.DataFrame({name: X[kos].max(axis=1) for name, kos in groups.items()})
    X.to_csv(OUT / "ko_presence.tsv", sep="\t")
    G.to_csv(OUT / "group_presence.tsv", sep="\t")
    print(f"genomes {len(X)}, marker KOs {X.shape[1]}, functional groups {G.shape[1]}")
    print(f"KOs assigned to >1 group: {sum(1 for k in X.columns if sum(k in v for v in groups.values()) > 1)}")

    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    lab = lab[lab.labelled] if "labelled" in lab.columns else lab
    traits = [c for c in lab.columns if c.startswith(("donor_", "acceptor_"))]
    lab = lab.set_index("genome_id")

    # --- functional group vs the curated donor/acceptor category naming the same physiology
    rows = []
    shared = [g for g in X.index if g in lab.index]
    for name, trait in GROUP_TRAIT.items():
        if trait not in traits:
            continue
        idx = shared
        rows.append(associate(G.loc[idx, name].to_numpy(), lab.loc[idx, trait].to_numpy(dtype=int),
                              ranks["genus"].reindex(idx).fillna("?").to_numpy(), rng,
                              functional_group=name, role=role[name], phenotype=trait,
                              phenotype_source="curated donor/acceptor labels", n_kos_in_group=len(groups[name])))
    gt = pd.DataFrame(rows)
    gt["q_bh_permutation"] = bh(gt.p_genus_permutation)
    gt.sort_values("mcc_genus_level", ascending=False).to_csv(OUT / "group_vs_trait.csv", index=False,
                                                              float_format="%.4f")

    # --- single-KO resolution: which KO carries each group's signal
    rows = []
    for name, trait in GROUP_TRAIT.items():
        if trait not in traits:
            continue
        for ko in groups[name]:
            x = X.loc[shared, ko].to_numpy()
            if x.sum() == 0:
                continue
            rows.append(associate(x, lab.loc[shared, trait].to_numpy(dtype=int),
                                  ranks["genus"].reindex(shared).fillna("?").to_numpy(), rng,
                                  functional_group=name, KO=ko, phenotype=trait,
                                  phenotype_source="curated donor/acceptor labels", n_kos_in_group=len(groups[name])))
    kt = pd.DataFrame(rows)
    kt["q_bh_permutation"] = bh(kt.p_genus_permutation)
    kt.sort_values(["functional_group", "mcc_genus_level"], ascending=[True, False]).to_csv(
        OUT / "ko_vs_trait.csv", index=False, float_format="%.4f")

    # --- oxygen requirement and habitat from Madin, joined by species name
    mad = pd.read_csv(ROOT / "data/raw/madin/condensed_species_NCBI.csv", low_memory=False)
    sp = info["species"].str.replace("Candidatus ", "", regex=False).str.split().str[:2].str.join(" ")
    mad_by_sp = mad.drop_duplicates("species").set_index("species")
    ox = sp.map(mad_by_sp.metabolism).map(lambda v: OXYGEN.get(str(v).strip().lower(), np.nan))
    iso = sp.map(mad_by_sp.isolation_source).fillna("").str.lower()
    hab = pd.DataFrame({k: iso.apply(lambda s, ws=v: int(any(s.startswith(w) for w in ws)))
                        for k, v in HABITAT.items()})
    print(f"joined to Madin: oxygen {int(ox.notna().sum())}/{len(sp)}, "
          f"isolation source {int((iso != '').sum())}/{len(sp)}")

    rows = []
    for name in G.columns:
        for cls in ("aerobe", "facultative", "anaerobe", "microaerophile"):
            m = ox.notna().to_numpy()
            if m.sum() < 30 or (ox[m] == cls).sum() < 10:
                continue
            rows.append(associate(G.loc[m, name].to_numpy(), (ox[m] == cls).astype(int).to_numpy(),
                                  ranks["genus"].reindex(G.index[m]).fillna("?").to_numpy(), rng,
                                  functional_group=name, role=role[name], phenotype=f"oxygen:{cls}",
                                  phenotype_source="Madin metabolism", n_kos_in_group=len(groups[name])))
    go = pd.DataFrame(rows)
    go["q_bh_permutation"] = bh(go.p_genus_permutation)
    go.sort_values("mcc_genus_level", ascending=False).to_csv(OUT / "group_vs_oxygen.csv", index=False,
                                                              float_format="%.4f")

    rows = []
    for name in G.columns:
        for h in HABITAT:
            m = (iso != "").to_numpy()
            if m.sum() < 30 or hab.loc[m, h].sum() < 10:
                continue
            rows.append(associate(G.loc[m, name].to_numpy(), hab.loc[m, h].to_numpy(),
                                  ranks["genus"].reindex(G.index[m]).fillna("?").to_numpy(), rng,
                                  functional_group=name, role=role[name], phenotype=f"habitat:{h}",
                                  phenotype_source="Madin isolation_source", n_kos_in_group=len(groups[name])))
    gh = pd.DataFrame(rows)
    gh["q_bh_permutation"] = bh(gh.p_genus_permutation)
    gh.sort_values("mcc_genus_level", ascending=False).to_csv(OUT / "group_vs_habitat.csv", index=False,
                                                              float_format="%.4f")

    json.dump({"n_genomes": int(len(X)), "n_marker_kos": int(X.shape[1]), "n_groups": int(G.shape[1]),
               "n_with_curated_labels": len(shared), "n_with_oxygen": int(ox.notna().sum()),
               "n_with_isolation_source": int((iso != "").sum()),
               "n_permutations": N_PERM, "n_genus_draws": N_GENUS_DRAWS,
               "kos_in_multiple_groups": sorted(k for k in X.columns
                                                if sum(k in v for v in groups.values()) > 1)},
              open(OUT / "summary.json", "w"), indent=2)

    show = ["functional_group", "phenotype", "n_marker_present", "n_phenotype_pos", "precision", "recall",
            "mcc_per_genome", "mcc_genus_level", "p_genus_permutation", "q_bh_permutation"]
    print("\n=== functional group vs curated donor/acceptor label ===")
    print(gt.sort_values("mcc_genus_level", ascending=False)[show].round(3).to_string(index=False))
    print("\n=== functional group vs oxygen requirement (top 12) ===")
    print(go.sort_values("mcc_genus_level", ascending=False)[show].head(12).round(3).to_string(index=False))
    print("\n=== functional group vs habitat (top 12) ===")
    print(gh.sort_values("mcc_genus_level", ascending=False)[show].head(12).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
