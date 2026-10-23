"""Test the three explanations offered for the marker-set errors, one experiment each.

The collaborator's note lists three possible causes when a marker group disagrees with the curated phenotype:
  (a) the definition of the detected KO is incomplete,
  (b) the curated values are wrong (too high a share of "experimental" data, wrong strain),
  (c) detection fails when genome quality is poor.
Each is a testable claim, and this script tests it rather than arguing about it.

Test A - shared KOs (cause a). Several KOs sit in two functional groups: narGH (K00370/K00371) is both nitrate
  reductase (denitrifiers) and nitrite oxidoreductase (NOB); amoABC (K10944-6) is both ammonia monooxygenase and
  particulate methane monooxygenase (pmoABC). The group is re-scored without its ambiguous KOs and the change in
  precision / MCC is reported. A large precision gain means the errors are definitional, not biological.

Test B - marker breadth (cause a, second form). For each group, the best achievable MCC using the group's own KOs
  as a set is compared with the best single KO in it. A group whose best single KO beats the "any KO present" rule
  is one whose definition is too broad; a group where no KO is specific is one whose definition is too narrow.

Test C - genome quality (cause c). False negatives (phenotype present, no marker detected) are compared with true
  positives on CheckM2 completeness and contamination, and on the number of proteins. If poor genomes drive the
  misses, completeness should be visibly lower for false negatives.

Test D - co-occurrence beyond taxonomy. Pairs of functional groups are tested for co-occurrence with a
  genus-permutation null, which keeps lineage structure. Pairs that survive are candidate real metabolic couplings
  (e.g. sulfur oxidation with denitrification) rather than shared ancestry.

usage: python src/41_marker_definition_tests.py
Outputs (results/jaccard/): test_a_shared_kos.csv, test_b_marker_breadth.csv, test_c_genome_quality.csv,
                            test_d_group_cooccurrence.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/jaccard"
SEED = 20260918

import importlib.util  # noqa: E402

spec = importlib.util.spec_from_file_location("jc", ROOT / "src/40_jaccard_correlations.py")
jc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(jc)

# KOs whose annotation is shared between two physiologies, and the group each one is *not* diagnostic for.
# narGH/nxrAB: the same enzyme run in both directions; amoABC/pmoABC: homologous monooxygenases.
AMBIGUOUS = {
    "Nitrite oxidizing (NOB)": ["K00370", "K00371"],                # also narGH, in every nitrate reducer
    "Methane oxidizing (methanotroph)": ["K10944", "K10945", "K10946"],  # also amoABC, in ammonia oxidisers
    "Nitrate-reducing (denitrifiers)": ["K00362", "K00363"],        # nirBD, assimilatory nitrite reductase
}


def main():
    rng = np.random.default_rng(SEED)
    X, groups, role, ranks, info = jc.load_html(ROOT / "jaccard_similarity.html")
    G = pd.DataFrame({n: X[k].max(axis=1) for n, k in groups.items()})
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    lab = lab[lab.labelled] if "labelled" in lab.columns else lab
    lab = lab.set_index("genome_id")
    shared = [g for g in X.index if g in lab.index]
    genera = ranks["genus"].reindex(shared).fillna("?").to_numpy()

    # ---------------- Test A: drop the ambiguous KOs
    rows = []
    for name, drop in AMBIGUOUS.items():
        trait = jc.GROUP_TRAIT[name]
        if trait not in lab.columns:
            continue
        y = lab.loc[shared, trait].to_numpy(dtype=int)
        keep = [k for k in groups[name] if k not in drop]
        before = X.loc[shared, groups[name]].max(axis=1).to_numpy()
        after = X.loc[shared, keep].max(axis=1).to_numpy() if keep else np.zeros(len(shared), dtype=int)
        b = jc.associate(before, y, genera, rng, variant="all KOs")
        a = jc.associate(after, y, genera, rng, variant=f"without {', '.join(drop)}")
        for r in (b, a):
            r.update(functional_group=name, phenotype=trait, dropped=", ".join(drop),
                     n_kos_used=len(groups[name]) if r["variant"] == "all KOs" else len(keep))
            rows.append(r)
    ta = pd.DataFrame(rows)
    cols = ["functional_group", "phenotype", "variant", "n_kos_used", "n_marker_present", "tp", "fp",
            "precision", "recall", "mcc_per_genome", "mcc_genus_level", "p_genus_permutation"]
    ta[cols].to_csv(OUT / "test_a_shared_kos.csv", index=False, float_format="%.4f")
    print("=== Test A: removing KOs that are shared between two physiologies ===")
    print(ta[cols].round(3).to_string(index=False))

    # ---------------- Test B: is the group definition too broad or too narrow?
    rows = []
    for name, kos in groups.items():
        trait = jc.GROUP_TRAIT[name]
        if trait not in lab.columns:
            continue
        y = lab.loc[shared, trait].to_numpy(dtype=int)
        anyko = X.loc[shared, kos].max(axis=1).to_numpy()
        m_any, _ = jc.mcc_phi(anyko, y)
        best_ko, best_mcc, best_prec = None, -2, np.nan
        for k in kos:
            x = X.loc[shared, k].to_numpy()
            if x.sum() == 0:
                continue
            m, (tp, fp, fn, tn) = jc.mcc_phi(x, y)
            if not np.isnan(m) and m > best_mcc:
                best_ko, best_mcc, best_prec = k, m, tp / (tp + fp) if tp + fp else np.nan
        rows.append(dict(functional_group=name, phenotype=trait, n_kos=len(kos),
                         mcc_any_ko_present=m_any, best_single_ko=best_ko, mcc_best_single_ko=best_mcc,
                         precision_best_single_ko=best_prec,
                         verdict=("definition too broad: one KO beats the group"
                                  if best_mcc > m_any + 0.02 else
                                  "group rule is the better detector" if m_any > best_mcc + 0.02 else "equivalent")))
    tb = pd.DataFrame(rows).sort_values("mcc_best_single_ko", ascending=False)
    tb.to_csv(OUT / "test_b_marker_breadth.csv", index=False, float_format="%.4f")
    print("\n=== Test B: group rule vs the best single KO in the group ===")
    print(tb.round(3).to_string(index=False))

    # ---------------- Test C: are the misses explained by genome quality?
    qcf = ROOT / "data/processed/genome_qc.tsv"
    if qcf.exists():
        qc = pd.read_csv(qcf, sep="\t").set_index("genome_id")
        qcols = [c for c in ("checkm2_completeness", "checkm2_contamination") if c in qc.columns]
        rows = []
        for name, kos in groups.items():
            trait = jc.GROUP_TRAIT[name]
            if trait not in lab.columns or not qcols:
                continue
            y = lab.loc[shared, trait].to_numpy(dtype=int)
            x = X.loc[shared, kos].max(axis=1).to_numpy()
            fn_idx = [g for g, xi, yi in zip(shared, x, y) if yi == 1 and xi == 0]
            tp_idx = [g for g, xi, yi in zip(shared, x, y) if yi == 1 and xi == 1]
            if len(fn_idx) < 5 or len(tp_idx) < 5:
                continue
            for c in qcols:
                a = qc[c].reindex(fn_idx).dropna().to_numpy()
                b = qc[c].reindex(tp_idx).dropna().to_numpy()
                if len(a) < 5 or len(b) < 5:
                    continue
                u, p = mannwhitneyu(a, b, alternative="two-sided")
                rows.append(dict(functional_group=name, phenotype=trait, metric=c,
                                 n_false_negatives=len(a), n_true_positives=len(b),
                                 median_false_negatives=float(np.median(a)), median_true_positives=float(np.median(b)),
                                 delta=float(np.median(a) - np.median(b)), p_mannwhitney=float(p)))
        tc = pd.DataFrame(rows)
        if len(tc):
            tc["q_bh"] = jc.bh(tc.p_mannwhitney)
            tc.to_csv(OUT / "test_c_genome_quality.csv", index=False, float_format="%.4f")
            print("\n=== Test C: genome quality of missed positives vs detected positives ===")
            print(tc.sort_values("p_mannwhitney").round(3).to_string(index=False))
    else:
        print("\n=== Test C skipped: no data/processed/genome_qc.tsv ===")

    # ---------------- Test D: functional groups co-occurring beyond shared ancestry
    names = list(G.columns)
    gen_all = ranks["genus"].reindex(G.index).fillna("?").to_numpy()
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            xa, xb = G[a].to_numpy(), G[b].to_numpy()
            if xa.sum() < 15 or xb.sum() < 15:
                continue
            p, phi = jc.genus_permutation_p(xa, xb, gen_all, rng, n=1000)
            both = int(((xa == 1) & (xb == 1)).sum())
            exp = xa.sum() * xb.sum() / len(xa)
            rows.append(dict(group_a=a, group_b=b, n_both=both, n_expected_independent=exp,
                             enrichment=both / exp if exp else np.nan, phi=phi, p_genus_permutation=p,
                             shares_kos=bool(set(groups[a]) & set(groups[b]))))
    td = pd.DataFrame(rows)
    td["q_bh"] = jc.bh(td.p_genus_permutation)
    td.sort_values("phi", ascending=False).to_csv(OUT / "test_d_group_cooccurrence.csv", index=False,
                                                  float_format="%.4f")
    print("\n=== Test D: group co-occurrence surviving a genus-permutation null (top 12 by phi) ===")
    print(td.sort_values("phi", ascending=False).head(12)[
        ["group_a", "group_b", "n_both", "enrichment", "phi", "p_genus_permutation", "q_bh", "shares_kos"]
    ].round(3).to_string(index=False))

    json.dump({"tests": ["A shared KOs", "B marker breadth", "C genome quality", "D co-occurrence"],
               "seed": SEED}, open(OUT / "definition_tests_summary.json", "w"), indent=2)


if __name__ == "__main__":
    main()
