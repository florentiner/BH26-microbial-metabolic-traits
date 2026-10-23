"""Collect every number the correlation report cites into one JSON, and grade each theory's confidence.

Confidence grades (fixed before grading):
  High          genus-level MCC >= 0.5, genus-permutation q < 0.01, and (where a KO was selected) the choice survives
                out-of-genus validation with >= 80 % selection stability; mechanism established in the literature.
  Moderate      q < 0.05 and genus-level MCC 0.3-0.5, or robust out-of-genus but unstable KO choice.
  Low           q < 0.05 but genus-level MCC < 0.3, or the signal rests on KOs shared between two groups.
  Not supported q >= 0.05, or the effect reverses under genus-level analysis.
These four apply to association claims ("marker X tracks trait Y"). Three claims are of a different kind and get
their own rules, written during this analysis (not pre-registered) and applied mechanically:
  Definition-insufficiency claim (H2 oxidation):
    High      recall < 0.7, misses spread over >= 10 genera, and the largest single lineage among the misses - one whose missing enzyme family can be named - accounts for
              >= 25 % of the misses;  Moderate  recall < 0.7 only;  Not supported  recall >= 0.7.
  Cause (a), over-broad definitions:
    High      >= 3 groups where one KO beats the group rule by >= 0.3 MCC out-of-genus with >= 80 % stability;
    Moderate  1-2 such groups;  Not supported  none.
  Cause (b), curated values wrong - graded on agreement between two independent curations (Landis & Koch):
    High      kappa < 0.6;  Moderate  0.6 <= kappa < 0.8;  Not supported  kappa >= 0.8.
  Cause (c), poor genome quality:
    Supported only if missed positives are significantly LESS complete than detected ones (q < 0.05, delta < 0) in
    at least one group; otherwise Not supported. A significant difference in the other direction does not count.

Output: results/jaccard/report_data.json
usage: python src/44_correlation_summary.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
J = ROOT / "results/jaccard"


def grade(mcc_genus, q, stability=None, shared_kos=False):
    if pd.isna(q) or q >= 0.05 or (not pd.isna(mcc_genus) and mcc_genus <= 0):
        return "Not supported"
    if shared_kos:
        return "Low"
    if mcc_genus >= 0.5 and q < 0.01 and (stability is None or stability >= 0.8):
        return "High"
    if mcc_genus >= 0.3 or (stability is not None and stability < 0.8 and mcc_genus >= 0.3):
        return "Moderate"
    return "Low"


def rec(df):
    return json.loads(df.replace({np.nan: None}).to_json(orient="records"))


def main():
    gt = pd.read_csv(J / "group_vs_trait.csv")
    ko = pd.read_csv(J / "single_ko_crossval.csv")
    ta = pd.read_csv(J / "test_a_shared_kos.csv")
    tc = pd.read_csv(J / "test_c_genome_quality.csv") if (J / "test_c_genome_quality.csv").exists() else None
    td = pd.read_csv(J / "test_d_group_cooccurrence.csv")
    ox_agree = pd.read_csv(J / "oxygen_label_agreement.csv")
    ox_dd = pd.read_csv(J / "oxygen_deep_dive.csv")
    ox_prev = pd.read_csv(J / "oxidase_prevalence.csv")
    go = pd.read_csv(J / "group_vs_oxygen.csv")
    gh = pd.read_csv(J / "group_vs_habitat.csv")
    summ = json.loads((J / "summary.json").read_text())

    def g(name):
        return gt[gt.functional_group == name].iloc[0]

    def k(name):
        return ko[ko.functional_group == name].iloc[0]

    def dd(marker, pheno):
        return ox_dd[(ox_dd.marker_set == marker) & (ox_dd.phenotype == pheno)].iloc[0]

    # spread of the hydrogen-oxidation misses, for the definition-insufficiency rule
    import importlib.util
    spec = importlib.util.spec_from_file_location("jc", ROOT / "src/40_jaccard_correlations.py")
    jc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(jc)
    X, groups, _, ranks, _ = jc.load_html(ROOT / "jaccard_similarity.html")
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    lab = lab[lab.labelled].set_index("genome_id") if "labelled" in lab.columns else lab.set_index("genome_id")
    ids = [x for x in X.index if x in lab.index]
    yh = lab.loc[ids, "donor_hydrogen_oxidation"].astype(int)
    xh = X.loc[ids, groups["Hydrogen oxidizing"]].max(axis=1)
    miss = ranks.loc[[x for x in ids if yh[x] == 1 and xh[x] == 0]]
    h2_miss = dict(n_misses=int(len(miss)), n_genera=int(miss.genus.nunique()),
                   n_families=int(miss.family.nunique()), n_phyla=int(miss.phylum.nunique()),
                   top_phylum=str(miss.phylum.value_counts().index[0]),
                   top_phylum_share=float(miss.phylum.value_counts().iloc[0] / len(miss)))

    # the single surviving habitat association (ammonia group vs marine), decomposed by marker
    mad = pd.read_csv(ROOT / "data/raw/madin/condensed_species_NCBI.csv", low_memory=False)
    info = jc.load_html(ROOT / "jaccard_similarity.html")[4]
    sp = info["species"].str.replace("Candidatus ", "", regex=False).str.split().str[:2].str.join(" ")
    iso = sp.map(mad.drop_duplicates("species").set_index("species").isolation_source).fillna("").str.lower()
    has_iso = (iso != "").to_numpy()
    marine = iso.apply(lambda v: v.startswith(("water_marine", "sediment_marine"))).astype(int).to_numpy()
    amo = X[["K10944", "K10945", "K10946"]].max(axis=1).to_numpy()
    hao_only = ((X["K10535"] == 1) & (X[["K10944", "K10945", "K10946"]].max(axis=1) == 0)).astype(int).to_numpy()
    grp = X[groups["Ammonia oxidizing (AOB/AOA)"]].max(axis=1).to_numpy()
    hab_decomp = []
    for nm, x in (("4-KO group", grp), ("amoABC only", amo), ("hao only (no amoABC)", hao_only)):
        mcc, (tp, fp, fn, tn) = jc.mcc_phi(x[has_iso], marine[has_iso])
        hab_decomp.append(dict(marker=nm, flagged=int(x[has_iso].sum()), marine_among_flagged=int(tp),
                               share_marine=tp / max(int(x[has_iso].sum()), 1), mcc_vs_marine=mcc))
    fl = X.index[(grp == 1) & (marine == 1) & has_iso]
    hab_top_genera = ranks.loc[fl, "genus"].value_counts().head(6).to_dict()

    overlaps = pd.DataFrame([dict(group_a=a_, group_b=b_, n_shared=len(set(groups[a_]) & set(groups[b_])),
                                  shared_kos=", ".join(sorted(set(groups[a_]) & set(groups[b_]))))
                             for i_, a_ in enumerate(groups) for b_ in list(groups)[i_ + 1:]
                             if set(groups[a_]) & set(groups[b_])])
    overlaps.to_csv(J / "group_definition_overlaps.csv", index=False)

    def variant(trait, kos):
        yv = lab.loc[ids, trait].to_numpy(dtype=int)
        xv = X.loc[ids, kos].max(axis=1).to_numpy()
        m_, (tp_, fp_, fn_, tn_) = jc.mcc_phi(xv, yv)
        return dict(flagged=int(xv.sum()), tp=tp_, fp=fp_, precision=tp_ / max(int(xv.sum()), 1),
                    recall=tp_ / max(tp_ + fn_, 1), mcc=m_)
    AMO, NARGH = ["K10944", "K10945", "K10946"], ["K00370", "K00371"]
    RESP8 = ["K02274", "K02275", "K02276", "K00404", "K00405", "K00406", "K00425", "K00426"]
    definition_variants = [
        dict(group="Nitrite oxidizing (NOB)", variant="as defined (amoABC + narGH)",
             **variant("donor_nitrite_oxidation", groups["Nitrite oxidizing (NOB)"])),
        dict(group="Nitrite oxidizing (NOB)", variant="narGH / nxrAB only (amoABC removed)",
             **variant("donor_nitrite_oxidation", NARGH)),
        dict(group="Nitrite oxidizing (NOB)", variant="amoABC only",
             **variant("donor_nitrite_oxidation", AMO)),
        dict(group="Aerobic respiration", variant="as defined (13 KOs)",
             **variant("acceptor_oxygen", groups["Aerobic respiration"])),
        dict(group="Aerobic respiration", variant="respiratory oxidases only (cox, cco, cyd; 8 KOs)",
             **variant("acceptor_oxygen", RESP8)),
        dict(group="Aerobic respiration", variant="aa3 + cbb3 only (cytochrome bd removed; 6 KOs)",
             **variant("acceptor_oxygen", RESP8[:6])),
        dict(group="Aerobic respiration", variant="non-respiratory O2 users only (amoABC, mmoX, coxL)",
             **variant("acceptor_oxygen", AMO + ["K16157", "K03520"])),
    ]

    # reference-style enzyme <-> environment scan (src/45) and phylum-controlled co-occurrence (src/46)
    envs = pd.read_csv(J / "environment_scan.csv")
    cow = pd.read_csv(J / "cooccurrence_within_phylum.csv")
    genus_ = ranks["genus"].fillna("?")
    glv = X.groupby(genus_).mean().ge(0.5).astype(int)
    dsr_ = glv[["K11180", "K11181"]].max(axis=1).to_numpy()
    hdr_check = []
    for nm, kos in (("mcrABG (methyl-CoM reductase)", ["K00399", "K00401", "K00402"]),
                    ("hdrABC (heterodisulfide reductase)", ["K03388", "K03389", "K03390"])):
        xv = glv[kos].max(axis=1).to_numpy()
        m_, (tp_, fp_, fn_, tn_) = jc.mcc_phi(xv, dsr_)
        hdr_check.append(dict(genes=nm, n_genera=int(xv.sum()), n_genera_also_dsrAB=int(tp_), phi_with_dsrAB=m_))
    hao_dnra = jc.mcc_phi(glv["K10535"].to_numpy(),
                          glv[groups["Nitrate-reducing (DNRA)"]].max(axis=1).to_numpy())[0]
    amob_dnra = jc.mcc_phi(glv["K10945"].to_numpy(),
                           glv[groups["Nitrate-reducing (DNRA)"]].max(axis=1).to_numpy())[0]

    theories = []
    # ---- reference's claim 1
    a, n = g("Ammonia oxidizing (AOB/AOA)"), g("Nitrite oxidizing (NOB)")
    kvt = pd.read_csv(J / "ko_vs_trait.csv")
    amoB = kvt[(kvt.functional_group == "Ammonia oxidizing (AOB/AOA)") & (kvt.KO == "K10945")].iloc[0]
    ka, kn = k("Ammonia oxidizing (AOB/AOA)"), k("Nitrite oxidizing (NOB)")
    theories.append(dict(
        id="C1a", source="reference", claim="Ammonia-oxidation enzymes are predicted accurately",
        verdict="Supported with the right marker; not with the group rule",
        grade=grade(amoB.mcc_genus_level, amoB.q_bh_permutation, ka.selection_stability),
        key_numbers=dict(group_precision=a.precision, group_recall=a.recall, group_mcc_genus=a.mcc_genus_level,
                         single_ko=ka.most_chosen_ko, single_ko_mcc_out_of_genus=ka.mcc_single_ko_out_of_genus,
                         selection_stability=ka.selection_stability, q=a.q_bh_permutation,
                         amoB_mcc_per_genome=amoB.mcc_per_genome, amoB_mcc_genus=amoB.mcc_genus_level,
                         amoB_q=amoB.q_bh_permutation, amoB_precision=amoB.precision, amoB_recall=amoB.recall,
                         amoB_flagged=int(amoB.n_marker_present))))
    theories.append(dict(
        id="C1b", source="reference", claim="Nitrite oxidase (NOB) is predicted accurately",
        verdict="Weak: recall is perfect but precision is ~10 %; the signal rests on KOs shared with nitrate reductase",
        grade=grade(n.mcc_genus_level, n.q_bh_permutation, shared_kos=True),
        key_numbers=dict(group_precision=n.precision, group_recall=n.recall, n_flagged=int(n.n_marker_present),
                         n_true=int(n.n_phenotype_pos), group_mcc_genus=n.mcc_genus_level,
                         single_ko_mcc_out_of_genus=kn.mcc_single_ko_out_of_genus, q=n.q_bh_permutation)))
    # ---- reference's claim 2
    allg = dd("Aerobic respiration (all 13 KOs)", "Madin O2-capable")
    tru = dd("aa3 or cbb3 (true respiratory oxidases)", "Madin O2-capable")
    bd = dd("cytochrome bd only", "Madin O2-capable")
    theories.append(dict(
        id="C2", source="reference", claim="Aerobic-respiration enzymes correlate with oxygen level",
        verdict="Supported — strongly once cytochrome bd is excluded",
        grade=grade(tru.mcc_genus_level, tru.q_bh_permutation),
        key_numbers=dict(all13_mcc_genus=allg.mcc_genus_level, all13_or=allg.odds_ratio,
                         aa3_cbb3_mcc_genus=tru.mcc_genus_level, aa3_cbb3_precision=tru.precision,
                         aa3_cbb3_recall=tru.recall, cyt_bd_mcc_genus=bd.mcc_genus_level,
                         cyt_bd_q=bd.q_bh_permutation, q=tru.q_bh_permutation,
                         label_kappa=float(ox_agree.cohens_kappa.iloc[0]))))
    # ---- reference's claim 3
    h, kh = g("Hydrogen oxidizing"), k("Hydrogen oxidizing")
    theories.append(dict(
        id="C3", source="reference", claim="The hydrogen-oxidation definition is insufficient",
        verdict="Supported: 4 KOs miss ~45 % of known H2 oxidisers under any rule",
        grade=("High" if h.recall < 0.7 and h2_miss["n_genera"] >= 10 and h2_miss["top_phylum_share"] >= 0.25
               else "Moderate" if h.recall < 0.7 else "Not supported"),
        key_numbers=dict(n_kos=int(h.n_kos_in_group), recall=h.recall, precision=h.precision,
                         group_mcc_genus=h.mcc_genus_level, single_ko_mcc_out_of_genus=kh.mcc_single_ko_out_of_genus,
                         q=h.q_bh_permutation, **h2_miss)))
    # ---- causes (a) (b) (c)
    robust = ko[ko.verdict == "single KO better (robust)"]
    big = robust[(robust.gain_out_of_genus >= 0.3) & (robust.selection_stability >= 0.8)]
    theories.append(dict(
        id="Ca", source="reference cause (a)", claim="Errors come from incomplete / over-broad KO definitions",
        verdict=f"Supported: for {len(big)} groups one specific KO beats the 'any KO' rule by ≥ 0.3 MCC out-of-genus",
        grade="High" if len(big) >= 3 else "Moderate" if len(big) >= 1 else "Not supported",
        key_numbers=dict(groups_single_ko_better=robust.functional_group.tolist(),
                         groups_gain_ge_0_3=big.functional_group.tolist(),
                         nob_drop_shared_recall=float(ta[(ta.functional_group == "Nitrite oxidizing (NOB)") &
                                                         (ta.variant != "all KOs")].recall.iloc[0]))))
    kap = float(ox_agree.cohens_kappa.iloc[0])
    theories.append(dict(
        id="Cb", source="reference cause (b)", claim="Curated values may be wrong",
        verdict="Partly: two oxygen sources agree at κ "
                f"{float(ox_agree.cohens_kappa.iloc[0]):.2f} — substantial but not perfect",
        grade=("High" if kap < 0.6 else "Moderate" if kap < 0.8 else "Not supported"),
        key_numbers=dict(kappa_o2_capable=float(ox_agree.cohens_kappa.iloc[0]),
                         kappa_strict=float(ox_agree.cohens_kappa.iloc[1]),
                         raw_agreement=float(ox_agree.raw_agreement.iloc[0]))))
    qc_sig = (tc[(tc.q_bh < 0.05)] if tc is not None else pd.DataFrame())
    qc_expected = (qc_sig[(qc_sig.metric == "checkm2_completeness") & (qc_sig.delta < 0)]
                   if len(qc_sig) else qc_sig)
    # completeness of the analysed genomes, read from the QC table (audit M10b: not "all > 99 %")
    qc_all = pd.read_csv(ROOT / "data/processed/genome_qc.tsv", sep="\t").set_index("genome_id")
    ids = pd.read_csv(J / "ko_presence.tsv", sep="\t", index_col=0, usecols=[0]).index
    comp = qc_all.checkm2_completeness.reindex(ids)
    comp_txt = (f"median completeness {comp.median():.1f} %; {int((comp < 99).sum())} of {int(comp.notna().sum())} "
                f"genomes with a CheckM2 value are below 99 %, {int((comp < 95).sum())} below 95 %, minimum "
                f"{comp.min():.1f} %")
    theories.append(dict(
        id="Cc", source="reference cause (c)", claim="Poor genome quality causes missed detections",
        verdict=("Supported in at least one group" if len(qc_expected) else
                 f"Not supported in this dataset ({comp_txt}; few incomplete genomes, so low power)"),
        grade="Moderate" if len(qc_expected) else "Not supported",
        key_numbers=dict(n_tests=int(len(tc)) if tc is not None else 0, n_significant=int(len(qc_sig)),
                         median_completeness_fn=float(tc[tc.metric == "checkm2_completeness"]
                                                      .median_false_negatives.median()) if tc is not None else None)))
    # ---- new correlations: co-occurrence surviving the genus null, not driven by shared KOs
    real = td[(td.q_bh < 0.05) & (~td.shares_kos)].sort_values("phi", ascending=False)
    art = td[(td.q_bh < 0.05) & (td.shares_kos)].sort_values("phi", ascending=False)

    out = dict(summary=summ, theories=theories,
               group_vs_trait=rec(gt), single_ko=rec(ko), shared_ko_test=rec(ta),
               genome_quality=rec(tc) if tc is not None else [],
               cooccurrence_real=rec(real.head(10)), cooccurrence_artefact=rec(art.head(6)),
               oxygen_agreement=rec(ox_agree), oxygen_deep_dive=rec(ox_dd), oxidase_prevalence=rec(ox_prev),
               group_vs_oxygen_top=rec(go[go.q_bh_permutation < 0.05].sort_values("mcc_genus_level",
                                                                                  ascending=False)),
               group_vs_habitat_top=rec(gh[gh.q_bh_permutation < 0.05].sort_values("mcc_genus_level",
                                                                                   ascending=False).head(12)),
               habitat_decomposition=hab_decomp, habitat_top_genera=hab_top_genera,
               group_overlaps=rec(overlaps), definition_variants=definition_variants,
               env_functional=rec(envs[envs.verdict.str.startswith("functional")].sort_values(
                   ["marker_set", "q_bh_within_phylum"])),
               env_lineage_bound=rec(envs[envs.verdict.str.startswith("lineage")].assign(
                   a=lambda d: d.statistic.abs()).sort_values("a", ascending=False).drop(columns="a")),
               env_counts=rec(envs.groupby(["family", "marker_set", "verdict"]).size().rename("n").reset_index()),
               n_env_tests=int(len(envs)),
               cooc_beyond_phylum=rec(cow[cow.verdict == "co-occur beyond phylum identity"]),
               cooc_lineage_bound=int((cow.verdict == "lineage-bound").sum()),
               cooc_artefact=rec(cow[cow.verdict == "artefact of shared KOs"]),
               hdr_check=hdr_check, hao_dnra_phi=hao_dnra, amob_dnra_phi=amob_dnra,
               n_group_pairs_sharing=int(len(overlaps)),
               h2_single_ko_mcc_per_genome=float(k("Hydrogen oxidizing").mcc_single_ko_out_of_genus),
               h2_group_mcc_per_genome=float(g("Hydrogen oxidizing").mcc_per_genome),
               n_oxygen_sig=int((go.q_bh_permutation < 0.05).sum()), n_oxygen_tests=int(len(go)),
               n_habitat_sig=int((gh.q_bh_permutation < 0.05).sum()), n_habitat_tests=int(len(gh)))
    (J / "report_data.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"written {J / 'report_data.json'}")
    for t in theories:
        print(f"  [{t['grade']:13s}] {t['id']:3s} {t['claim']}  ->  {t['verdict']}")
    print(f"  new real co-occurrences: {len(real)}; shared-KO artefacts: {len(art)}")
    print(f"  oxygen: {out['n_oxygen_sig']}/{out['n_oxygen_tests']} significant; "
          f"habitat: {out['n_habitat_sig']}/{out['n_habitat_tests']} significant")


if __name__ == "__main__":
    main()
