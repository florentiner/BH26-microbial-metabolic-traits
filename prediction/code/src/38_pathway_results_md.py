"""Render docs/RESULTS_pathways.md from the result CSVs. Every number is read from a file, never typed by hand.

usage: python src/38_pathway_results_md.py
"""
import json
from pathlib import Path

import pandas as pd

from pathway_common import RES, load_pathways

OUT = Path(__file__).resolve().parents[1] / "docs/RESULTS_pathways.md"
ARM_NAME = {"PA": "KEGG module completeness rule (baseline)", "PB1": "852-KO energy panel",
            "PB2": "full KO profile", "PB3": "full KO profile + module completeness",
            "PB4": "multi-label MLP on the full KO profile", "PE": "KO profile + Bacformer embedding",
            "null_prevalence": "prevalence null", "null_taxNN": "taxonomic group-mean null"}


def fmt(x, n=3):
    return "-" if pd.isna(x) else f"{x:.{n}f}"


def main():
    lab, traits, Y = load_pathways()
    macro = pd.read_csv(RES / "pathways/metrics_macro.csv")
    per_trait = pd.read_csv(RES / "pathways/metrics_per_trait.csv")
    paired = pd.read_csv(RES / "pathways/paired_comparisons.csv")
    L = []
    A = L.append

    A("# Results — energy metabolism pathway prediction\n")
    A(f"Generated from `results/pathways/` by `src/38_pathway_results_md.py`. Pre-registration: "
      f"`docs/preregistration_pathways.md`. Labels frozen at `{json.loads((RES / 'pathways/audit.json').read_text())['labels_sha256'][:16]}…`.\n")
    A(f"**Dataset.** {len(lab)} genomes (one RefSeq genome per species), {len(traits)} evaluable pathway traits, "
      f"{lab.group_genus.nunique()} genus groups, {lab.group_family.nunique()} family groups. Features: KofamScan KO "
      f"profiles from the RDF Portal BH26 store (the resource TogoMCP exposes). Labels: Madin et al. 2020, "
      f"literature-curated, independent of the features.\n")

    A("## Macro AUPRC by arm\n")
    A("| arm | features | genus-blocked | family-blocked | Δ vs taxonomic null (genus) | traits > null (q<0.05) |")
    A("|---|---|---|---|---|---|")
    g = macro[macro.split_level == "genus"].set_index("arm")
    f = macro[macro.split_level == "family"].set_index("arm")
    for arm in [a for a in ("PB4", "PB3", "PB2", "PB1", "PA", "null_taxNN", "null_prevalence") if a in g.index]:
        A(f"| {arm} | {ARM_NAME.get(arm, arm)} | {fmt(g.at[arm, 'macro_auprc'])} "
          f"[{fmt(g.at[arm, 'macro_auprc_ci_low'])}, {fmt(g.at[arm, 'macro_auprc_ci_high'])}] | "
          f"{fmt(f.at[arm, 'macro_auprc']) if arm in f.index else '-'} "
          f"[{fmt(f.at[arm, 'macro_auprc_ci_low']) if arm in f.index else '-'}, "
          f"{fmt(f.at[arm, 'macro_auprc_ci_high']) if arm in f.index else '-'}] | "
          f"{fmt(g.at[arm, 'macro_delta_vs_null'])} | "
          f"{int(g.at[arm, 'n_traits_better_than_null_q_lt_0_05'])} / {len(traits)} |")
    A("")

    A("## Pre-registered contrasts (cluster-bootstrap CIs over split groups)\n")
    A("| contrast | level | Δ macro AUPRC | 95 % CI | tied |")
    A("|---|---|---|---|---|")
    for r in paired[paired.trait == "MACRO"].itertuples():
        A(f"| {r.contrast} | {r.split_level} | {fmt(r.delta_auprc)} | [{fmt(r.ci_low)}, {fmt(r.ci_high)}] | "
          f"{'yes' if r.tied else 'no'} |")
    A("")

    best = g.drop(index=[a for a in ("null_taxNN", "null_prevalence", "PA") if a in g.index]).macro_auprc.idxmax()
    A(f"## Per trait — {best} ({ARM_NAME.get(best, best)})\n")
    A("| trait | n pos | prevalence | AUPRC genus | AUPRC family | lift | MCC (family) | rule baseline PA | "
      "beats null q |")
    A("|---|---|---|---|---|---|---|---|---|")
    pg = per_trait[(per_trait.arm == best) & (per_trait.split_level == "genus")].set_index("trait")
    pf = per_trait[(per_trait.arm == best) & (per_trait.split_level == "family")].set_index("trait")
    pa = per_trait[(per_trait.arm == "PA") & (per_trait.split_level == "genus")].set_index("trait")
    for t in pg.sort_values("auprc", ascending=False).index:
        A(f"| {t} | {int(pg.at[t, 'n_pos'])} | {fmt(pg.at[t, 'prevalence'])} | {fmt(pg.at[t, 'auprc'])} "
          f"[{fmt(pg.at[t, 'auprc_ci_low'])}, {fmt(pg.at[t, 'auprc_ci_high'])}] | "
          f"{fmt(pf.at[t, 'auprc']) if t in pf.index else '-'} | {fmt(pg.at[t, 'auprc_lift'], 1)}x | "
          f"{fmt(pf.at[t, 'mcc']) if t in pf.index else '-'} | {fmt(pa.at[t, 'auprc']) if t in pa.index else '-'} | "
          f"{fmt(pg.at[t, 'delta_vs_null_q_bh'], 4)} |")
    A("")

    mag = RES / "pathways/mag_robustness_macro.csv"
    if mag.exists():
        m = pd.read_csv(mag)
        A("## MAG robustness (incomplete proteomes)\n")
        A("| proteome completeness | macro AUPRC | sd over draws | Δ vs complete |")
        A("|---|---|---|---|")
        for r in m.sort_values("completeness", ascending=False).itertuples():
            A(f"| {int(r.completeness * 100)} % | {fmt(r.macro_auprc)} | {fmt(r.macro_auprc_sd_over_draws)} | "
              f"{fmt(r.delta_vs_complete)} |")
        A("")

    cc = RES / "pathways/crosscheck_donor_acceptor.csv"
    if cc.exists():
        c = pd.read_csv(cc)
        c = c[c.split_level == "genus"]
        A("## Cross-source check against the independently curated donor/acceptor labels\n")
        A("Out-of-fold pathway predictions (genus held out of pathway training) scored against a different "
          "curation of the same physiology. Agreement here cannot come from label leakage.\n")
        A("| pathway trait | donor/acceptor category | n matched | AUPRC vs independent labels | lift | "
          "label agreement |")
        A("|---|---|---|---|---|---|")
        for r in c.itertuples():
            A(f"| {r.pathway_trait} | {r.donor_acceptor_trait} | {int(r.n_matched)} | "
              f"{fmt(r.auprc_vs_external)} | {fmt(r.auprc_lift, 1)}x | {fmt(r.label_agreement)} |")
        A("")

    for sub, title in (("no_faprotax", "SP1 — FAPROTAX-derived labels excluded from the evaluation set"),
                       ("informative_negatives", "SP2 — genomes with ≥ 2 annotated pathway terms only")):
        p = RES / f"pathways/sensitivity_{sub}/metrics_macro.csv"
        if p.exists():
            s = pd.read_csv(p)
            s = s[s.split_level == "genus"].set_index("arm")
            A(f"## Sensitivity {title}\n")
            A("| arm | macro AUPRC | 95 % CI |")
            A("|---|---|---|")
            for arm in [a for a in ("PB4", "PB3", "PB2", "PB1", "PA", "null_taxNN") if a in s.index]:
                A(f"| {arm} | {fmt(s.at[arm, 'macro_auprc'])} | [{fmt(s.at[arm, 'macro_auprc_ci_low'])}, "
                  f"{fmt(s.at[arm, 'macro_auprc_ci_high'])}] |")
            A("")

    rel = RES / "models/pathways/RELEASED.json"
    if rel.exists():
        r = json.loads(rel.read_text())
        A("## Released model\n")
        A(f"`{r['path']}` (arm {r['arm']}), {r['n_reliable_traits']} of {r['n_traits']} traits pass the reliability "
          f"gate (family-level MCC CI lower bound > 0.2 and beats the taxonomic null at q < 0.05). "
          f"`python src/predict.py <accession>` returns donors, acceptors and pathways from one genome input.\n")

    audit = json.loads((RES / "pathways/audit.json").read_text())
    A("## Integrity\n")
    A(f"`src/37_pathway_audit.py`: {'all checks pass' if not audit['failures'] else 'FAILURES: ' + str(audit['failures'])}. "
      f"Labels hashed before training; genus and family groups never span folds; features are KO ids only "
      f"(no taxonomy, no organism name); FAPROTAX share reported per trait in "
      f"`results/pathways/qc_faprotax_share.csv`.\n")

    OUT.write_text("\n".join(L))
    print(f"written {OUT} ({len(L)} lines)")


if __name__ == "__main__":
    main()
