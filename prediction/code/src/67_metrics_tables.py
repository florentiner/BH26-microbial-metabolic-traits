"""One file with every metric table, split only by task (docs/all_metrics_tables.md).

Collects the classical-ML, deep-learning, LLM-baseline, rule and null arms of all three tasks from the result files
they were written to, and renders them as markdown tables. No new computation: every number is read from
results/, results/trees_xai/, results/pathways/ and validation_new_organisms/.

usage: python src/67_metrics_tables.py
"""
import importlib
from pathlib import Path

import pandas as pd

tx = importlib.import_module("47_trees_xai")
RES = tx.RES
ROOT = RES.parent
OUT = ROOT / "docs/all_metrics_tables.md"

# method class and one-line definition; verified against src/06_run_arms.py, 32_pathway_arms.py, 47_trees_xai.py
ARMS = {
    "null_prevalence": ("null", "Predict each category's overall frequency"),
    "null_taxNN": ("null", "Predict the mean label of the nearest taxonomic group in training"),
    "A": ("rule", "KEGG module completeness rule, threshold picked on the training fold"),
    "PA": ("rule", "KEGG module completeness rule (pathway task)"),
    "DT": ("ML", "Decision tree, depth 2-6"),
    "DT3": ("ML", "Decision tree, depth <= 3"),
    "FIGS": ("ML", "Fast interpretable greedy-tree sums"),
    "L1": ("ML", "L1 logistic regression"),
    "RF": ("ML", "Random forest"),
    "ET": ("ML", "Extra-trees"),
    "AM": ("ML", "AutoML (FLAML searches LightGBM/XGBoost/RF/extra-trees/logistic, 30 s per fit)"),
    "B1": ("ML", "Provided marker KOs"),
    "B1r": ("ML", "Marker KOs re-annotated by KofamScan"),
    "B2": ("ML", "852-KO energy panel, logistic or XGBoost (released model)"),
    "B2cnt": ("ML", "Same panel as KO counts instead of presence"),
    "B2tax": ("ML", "B2 + leakage-free taxonomic prior"),
    "B2RP": ("ML", "B2 + probe on a random projection of the KO panel (embedding-shaped control)"),
    "B2B1r": ("ML", "B2 + marker-subset model (no embeddings; ensembling control)"),
    "B2XL": ("ML", "B2 with logistic and XGBoost averaged"),
    "LG": ("ML", "LithoGenie metal/energy-cycling HMM features"),
    "B2LG": ("ML", "B2 + LithoGenie"),
    "SEL": ("ML", "Feature set chosen per category by inner CV from {B2, LG, B2+LG}"),
    "STACK": ("ML", "Cross-category stacking on B2 + LithoGenie"),
    "ENS": ("ML", "Ensemble of the KO arms"),
    "ENSH": ("ML", "Ensemble + multiheme (CXXCH) motif counts"),
    "PB1": ("ML", "852-KO energy panel (pathway task)"),
    "PB2": ("ML", "Full significant-KO profile"),
    "PB3": ("ML", "Full KO profile + KEGG module completeness (released pathway model)"),
    "C": ("DL", "ESM-2 protein-embedding probe"),
    "C_large": ("DL", "ESM-C 300M embedding probe"),
    "D": ("DL", "Attention-pooled ESM-2 (ABMIL)"),
    "E": ("DL", "Bacformer-base genome-embedding probe"),
    "E_MAG": ("DL", "Bacformer-base (MAG weights) probe"),
    "E_large": ("DL", "Bacformer-large genome-embedding probe"),
    "E_large_MAG": ("DL", "Bacformer-large (MAG weights) probe"),
    "NN": ("DL", "Neural network (MLP) on the KO profile"),
    "PB4": ("DL", "Multi-label MLP on the KO profile, all traits jointly"),
    "TPFN": ("DL", "TabPFN v2 tabular foundation model on the KO panel"),
    "MT": ("DL", "Multi-task network: 34 labels, donor/acceptor + pathway genomes together"),
    "ST": ("DL", "Same network, donor/acceptor labels only (transfer control)"),
    "FM_BFL_MLP": ("DL", "Bacformer-large vector -> MLP"),
    "FT_BFL0": ("DL", "Bacformer-large frozen + trained head"),
    "FT_BFL2": ("DL", "Bacformer-large fine-tuned (top 2 layers)"),
    "FTP_BFL0": ("DL", "Bacformer-large frozen + head (pathway task)"),
    "FTP_BFL2": ("DL", "Bacformer-large fine-tuned (pathway task)"),
    "PE_BFL": ("DL", "Bacformer-large probe (pathway task)"),
    "PE_MLP": ("DL", "Multi-label MLP on the Bacformer-large vector (pathway task)"),
    "B2C": ("ML+DL", "B2 + ESM-2 probe"),
    "B2E": ("ML+DL", "B2 + Bacformer-base probe"),
    "B2E_large": ("ML+DL", "B2 + Bacformer-large probe"),
    "B2XLC": ("ML+DL", "B2XL + ESM-2 probe"),
    "B2XLE": ("ML+DL", "B2XL + Bacformer-base probe"),
    "B2XLE_large": ("ML+DL", "B2XL + Bacformer-large probe"),
    "B2LGE": ("ML+DL", "B2 + LithoGenie + Bacformer-base probe"),
    "PE": ("ML+DL", "PB2 + Bacformer-large probe (pathway task)"),
    "FTP_FUSE": ("ML+DL", "PB2 + fine-tuned Bacformer-large (pathway task)"),
}


def collect(files):
    """Merge macro tables, keeping the first file's numbers when an arm appears twice."""
    frames = []
    for f in files:
        if Path(f).exists():
            d = pd.read_csv(f)
            d["_src"] = str(f)
            frames.append(d)
    d = pd.concat(frames, ignore_index=True).drop_duplicates(["arm", "split_level"], keep="first")
    w = d.pivot(index="arm", columns="split_level")
    out = pd.DataFrame(index=w.index)
    for lvl in ("genus", "family"):
        out[f"{lvl}_auprc"] = w[("macro_auprc", lvl)]
        out[f"{lvl}_lo"] = w[("macro_auprc_ci_low", lvl)]
        out[f"{lvl}_hi"] = w[("macro_auprc_ci_high", lvl)]
        out[f"{lvl}_mcc"] = w[("macro_mcc", lvl)]
        if ("macro_auroc", lvl) in w.columns:
            out[f"{lvl}_auroc"] = w[("macro_auroc", lvl)]
    return out.reset_index()


def table(df, note=""):
    rows = ["| Arm | Type | Macro AUPRC, held-out genera [95% CI] | MCC | Macro AUPRC, held-out families [95% CI] | MCC | What it is |",
            "|---|---|---|---|---|---|---|"]
    df = df.sort_values("genus_auprc", ascending=False)
    for _, r in df.iterrows():
        kind, desc = ARMS.get(r["arm"], ("", ""))
        def cell(lvl):
            if pd.isna(r.get(f"{lvl}_auprc")):
                return "—", "—"
            m = r.get(f"{lvl}_mcc")
            return (f"{r[f'{lvl}_auprc']:.3f} [{r[f'{lvl}_lo']:.3f}, {r[f'{lvl}_hi']:.3f}]",
                    "—" if pd.isna(m) else f"{m:.3f}")
        g, gm = cell("genus")
        f, fm = cell("family")
        rows.append(f"| {r['arm']} | {kind} | {g} | {gm} | {f} | {fm} | {desc} |")
    return "\n".join(rows) + ("\n\n" + note if note else "")


def main():
    parts = ["# All metrics, by task",
             "",
             "Every number is read from the result files; nothing is recomputed here. Type: ML = classical machine "
             "learning on gene content, DL = deep learning / foundation models, ML+DL = fusion of both, rule = "
             "deterministic marker rule, null = baseline, LLM = large language model. All cross-validation is "
             "grouped and frozen: no genus (or family) appears in both training and test.",
             "",
             "**Reading the MCC columns.** AUPRC is threshold-free and comparable across every arm. MCC is not: it "
             "depends on how each arm picked its cut-off. Most arms use a 4-fold grouped inner out-of-fold "
             "estimate; TPFN, MT and ST use a single inner split, which costs them roughly 0.10 MCC (at a "
             "well-chosen threshold those three are the highest-MCC arms in the table, not the lowest). Rank arms "
             "by AUPRC, and treat their MCC as a property of the threshold estimator rather than of the model.",
             "",
             "---",
             "",
             "## Task 1 — Electron donors and acceptors",
             "14 categories, 358 labelled genomes. Primary metric: macro AUPRC over the 14 categories.",
             ""]
    da = collect([RES / "metrics_macro.csv", RES / "trees_xai/da/metrics_macro.csv"])
    parts += ["### 1.1 Cross-validated comparison of every arm", "", table(da), ""]

    llm = pd.read_csv(RES / "baseline2_llm_summary.csv")
    NAME = {"H_genome": "Claude, genome accession only", "H_named": "Claude, organism name given",
            "H_gpt_genome": "GPT, genome accession only", "H_gpt_named": "GPT, organism name given"}
    rows = ["| LLM arm | Genomes | LLM macro AUPRC | Compared with | Model macro AUPRC | LLM − model [95% CI] |",
            "|---|---|---|---|---|---|"]
    for a, g in llm.groupby("llm_arm", sort=False):
        self_row = g[g.model == a].iloc[0]
        for _, r in g[g.model != a].iterrows():
            d = "—" if pd.isna(r.llm_minus_model) else f"{r.llm_minus_model:+.3f} [{r.ci_low:+.3f}, {r.ci_high:+.3f}]"
            rows.append(f"| {NAME.get(a, a)} | {int(r.n_genomes)} | {self_row.macro_auprc:.3f} | {r.model} | "
                        f"{r.macro_auprc:.3f} | {d} |")
    parts += ["### 1.2 Large language models vs the released models",
              "Held-out subset of genomes, same categories. A negative difference means the LLM is worse.", "",
              "\n".join(rows), ""]

    ext = pd.read_csv(RES / "external_validation_macro.csv")
    ext = ext[ext.subset == "all"].sort_values("macro_auprc", ascending=False)
    rows = ["| Arm | Type | Macro AUPRC | Macro AUROC | Macro MCC | Lift over prevalence |", "|---|---|---|---|---|---|"]
    for _, r in ext.iterrows():
        kind = ARMS.get(r.arm, ("", ""))[0]
        mcc = "—" if pd.isna(r.macro_mcc) else f"{r.macro_mcc:.3f}"
        rows.append(f"| {r.arm} | {kind} | {r.macro_auprc:.3f} | {r.macro_auroc:.3f} | {mcc} | {r.macro_lift:.1f}× |")
    parts += ["### 1.3 External validation — 562 lineage-novel genomes, 9 categories", "", "\n".join(rows), ""]

    v = ROOT / "validation_new_organisms"
    if (v / "scored_B2.csv").exists():
        rows = ["| Model | Scored calls | Correct | Tested positives found | Tested negatives respected | "
                "The 8 categories with released yes/no calls |", "|---|---|---|---|---|---|"]
        REL = {"acceptor_CO2_methanogen_acetogen", "acceptor_nitrate", "acceptor_oxygen", "acceptor_sulfate",
               "acceptor_sulfur", "donor_ammonia_oxidation", "donor_hydrogen_oxidation", "donor_sulfur_oxidation"}
        for arm in ("B2", "B2XLE_large"):
            d = pd.read_csv(v / f"scored_{arm}.csv")
            s = d[d.expected.notna()]
            pos, neg = s[s.expected == 1], s[s.expected == 0]
            r = s[s.category.isin(REL)]
            rows.append(f"| {arm} | {len(s)} | {int(s.agree.sum())}/{len(s)} ({s.agree.mean():.0%}) | "
                        f"{int(pos.agree.sum())}/{len(pos)} | {int(neg.agree.sum())}/{len(neg)} | "
                        f"{int(r.agree.sum())}/{len(r)} ({r.agree.mean():.0%}) |")
        parts += ["### 1.4 Prospective test — the first 10 newly described organisms", "", "\n".join(rows), ""]
    vs = v / "validation_summary.csv"
    if vs.exists():
        d = pd.read_csv(vs)
        rows = ["| Analysis set | Model | Calls | Accuracy [95% CI] | MCC | Balanced acc. | Sensitivity | "
                "Specificity | Deployable baseline |", "|---|---|---|---|---|---|---|---|---|"]
        for _, r in d.iterrows():
            rows.append(f"| {r.label} | {r.arm} | {int(r.n)} | {r.accuracy:.1%} [{r.ci_low:.1%}, {r.ci_high:.1%}] | "
                        f"{r.mcc:.3f} | {r.balanced_accuracy:.1%} | {r.sensitivity:.1%} | {r.specificity:.1%} | "
                        f"{r.baseline_train_prior:.1%} |")
        parts += ["### 1.5 Prospective test — 162 newly described organisms (scaled, audited)",
                  "Primary analysis = blind evidence only, excluding the two catch-all categories. The deployable "
                  "baseline is a constant rule using training prevalence; confidence intervals are cluster "
                  "bootstraps over description papers. See validation_new_organisms/REPORT_prospective_167.md.",
                  "", "\n".join(rows), ""]

    parts += ["---", "", "## Task 2 — Energy-metabolism pathways",
              "20 traits, 3 378 labelled genomes. Primary metric: macro AUPRC over the 20 traits.", ""]
    pw = collect([RES / "pathways/metrics_macro.csv", RES / "trees_xai/pw/metrics_macro.csv"])
    parts += ["### 2.1 Cross-validated comparison of every arm", "", table(pw), ""]

    pl = RES / "pathway_llm_summary.csv"
    if pl.exists():
        sm = pd.read_csv(pl).set_index("arm").macro_auprc
        ct = pd.read_csv(RES / "pathway_llm_contrasts.csv")
        full = collect([RES / "pathways/metrics_macro.csv", RES / "trees_xai/pw/metrics_macro.csv"]) \
            .set_index("arm").genus_auprc
        NAMEP = {"P_H_genome": "Claude, gene evidence only", "P_H_named": "Claude, organism name given"}
        rows = ["| Arm | Macro AUPRC on the 96-genome subset | Same arm, all 3 378 genomes |", "|---|---|---|"]
        for a, v in sm.sort_values(ascending=False).items():
            rows.append(f"| {NAMEP.get(a, a)} | {v:.3f} | "
                        f"{('%.3f' % full[a]) if a in full.index and pd.notna(full.get(a)) else '—'} |")
        rows += ["", "| LLM arm | Compared with | LLM − model [95% CI] |", "|---|---|---|"]
        for _, r in ct.iterrows():
            rows.append(f"| {NAMEP.get(r.llm_arm, r.llm_arm)} | {r.model} | "
                        f"{r.delta:+.3f} [{r.ci_low:+.3f}, {r.ci_high:+.3f}] |")
        parts += ["### 2.2 Large language models on the pathway task",
                  "The 96-genome subset is enriched for rare traits, so its numbers sit well above the full-set "
                  "numbers; only the subset column is internally comparable. PB1 is the matched comparator (the "
                  "LLM saw exactly PB1's 852-KO evidence); PB3 is the released pathway model.",
                  "", "\n".join(rows), ""]

    parts += ["---", "", "## Task 3 — Oxygen use",
              "1 trait (aerobic vs anaerobic), 2 927 labelled genomes. Included as the calibration point: this trait "
              "is nearly saturated, and every method reaches it.", ""]
    ox = collect([RES / "trees_xai/ox/metrics_macro.csv"])
    parts += ["### 3.1 Cross-validated comparison of every arm", "", table(ox), ""]

    OUT.write_text("\n".join(parts) + "\n")
    print(OUT, f"{len(da)} donor/acceptor arms, {len(pw)} pathway arms, {len(ox)} oxygen arms")


if __name__ == "__main__":
    main()
