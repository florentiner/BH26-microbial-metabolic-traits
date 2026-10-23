"""Baseline 2: LLM over TogoMCP KofamScan evidence (docs/baseline2_llm_togomcp.md).

usage:
  python src/21_baseline2_llm.py prepare     # subset + anonymised evidence batches -> data/baseline2/batches/
  python src/21_baseline2_llm.py score       # collect data/baseline2/out/*.json -> results/baseline2_*.csv
Arms: H_genome / H_named = Claude sub-agents; H_gpt_genome / H_gpt_named = OpenAI GPT-6 Astra via OpenRouter.
"""
import json
import sys
import numpy as np
import pandas as pd

from redox_common import RES, ROOT, load_labelled, oof_matrix, paired_macro_delta, safe_ap

B2DIR = ROOT / "data/baseline2"
SEED = 20260917
N_SUBSET = 96
BATCH = 12
CATEGORIES = {
    "donor_sulfur_oxidation": "Electron donor: reduced sulfur compounds (sulfur-oxidizing bacteria/archaea)",
    "donor_ammonia_oxidation": "Electron donor: ammonia (ammonia-oxidizing bacteria/archaea)",
    "donor_nitrite_oxidation": "Electron donor: nitrite (nitrite-oxidizing bacteria)",
    "donor_hydrogen_oxidation": "Electron donor: H2 (hydrogen-oxidizing)",
    "donor_CO_oxidation": "Electron donor: carbon monoxide (carboxydotrophs)",
    "donor_iron_oxidation": "Electron donor: Fe(II) (iron-oxidizing)",
    "donor_other": "Electron donor: other inorganic/C1 donors not listed above (e.g. formate, methane, methanol, arsenite, Mn(II))",
    "acceptor_oxygen": "Electron acceptor: O2 (aerobic respiration)",
    "acceptor_nitrate": "Electron acceptor: nitrate/nitrite (denitrification or DNRA)",
    "acceptor_sulfate": "Electron acceptor: sulfate (sulfate-reducing)",
    "acceptor_sulfur": "Electron acceptor: elemental sulfur / thiosulfate / polysulfide (sulfur-reducing)",
    "acceptor_CO2_methanogen_acetogen": "Electron acceptor: CO2 (methanogens or acetogens)",
    "acceptor_iron": "Electron acceptor: Fe(III) (iron-reducing)",
    "acceptor_other": "Electron acceptor: other (e.g. arsenate, selenate, Mn(IV), fumarate, DMSO, perchlorate, nitrite by anammox)",
}
PROMPT = """You are an expert microbial physiologist. For each genome below you get ONLY genome annotation evidence:
significant KofamScan hits (KEGG Orthology) in energy metabolism, retrieved from the RDF Portal BH26 Microbes
database (the data behind TogoMCP). Each line: KO | KEGG definition | number of proteins | best score/threshold.

Task: for every genome and every category, give the probability (0.00-1.00) that the organism uses that electron
donor / electron acceptor for energy metabolism, as it would be recorded in a curated electron donor/acceptor
trait database. Reason from the gene content (complete pathways/complexes, marker specificity, ambiguous KOs such as
narGH vs nxrAB or amoABC vs pmoABC). Do NOT use any tool, file, web search or outside lookup: judge from the evidence
given. Genomes are anonymised; do not guess identities.

Categories (key: meaning):
{categories}

Return ONLY a JSON object: {{"<genome_code>": {{"<category_key>": probability, ...}}, ...}} with all 14 keys for every
genome. Write it to the output path given by the orchestrator.
"""


def pick_subset(lab, traits, Y):
    """Only genomes present in the BH26 Microbes endpoint (70/586 genomes have no TogoMCP statements)."""
    rng = np.random.default_rng(SEED)
    log = pd.read_csv(ROOT / "data/processed/togo_fetch_log_internal.tsv", sep="\t")
    available = set(log.loc[log.status == "ok", "genome_id"])
    eligible = set(np.where(lab.genome_id.isin(available))[0])
    chosen = []
    order = np.argsort(Y.sum(0))  # rarest trait first
    for j in order:
        pos = [i for i in np.where(Y[:, j] == 1)[0] if i not in chosen and i in eligible]
        need = max(0, 5 - int(Y[chosen, j].sum())) if chosen else 5
        chosen += list(rng.choice(pos, min(need, len(pos)), replace=False))
    rest = [i for i in sorted(eligible) if i not in chosen]
    chosen += list(rng.choice(rest, N_SUBSET - len(chosen), replace=False))
    return sorted(int(i) for i in chosen)


def prepare():
    lab, traits, Y = load_labelled()
    idx = pick_subset(lab, traits, Y)
    panel = pd.read_csv(ROOT / "data/processed/ko_panel.tsv", sep="\t").set_index("KO").definition
    togo = pd.read_parquet(ROOT / "data/processed/togo_kofam_internal.parquet")
    togo = togo[togo.significant & togo.KO.isin(panel.index)]
    rng = np.random.default_rng(SEED + 1)
    codes = {lab.genome_id[i]: f"G{c:03d}" for i, c in zip(idx, rng.permutation(len(idx)))}
    B2DIR.mkdir(parents=True, exist_ok=True)
    (B2DIR / "batches").mkdir(exist_ok=True)
    (B2DIR / "out").mkdir(exist_ok=True)
    key = lab.loc[idx, ["genome_id", "organism_name"]].assign(code=lambda d: d.genome_id.map(codes))
    key.to_csv(B2DIR / "subset_key.tsv", sep="\t", index=False)  # never given to the LLM in H-genome
    cats = "\n".join(f"- {k}: {v}" for k, v in CATEGORIES.items())
    for variant in ("genome", "named"):
        ordered = key.sort_values("code").reset_index(drop=True)
        for b in range(0, len(ordered), BATCH):
            lines = [PROMPT.format(categories=cats)]
            for _, r in ordered.iloc[b:b + BATCH].iterrows():
                g = togo[togo.genome_id == r.genome_id]
                agg = (g.assign(ratio=g.score / g.threshold).groupby("KO")
                       .agg(n=("protein", "nunique"), best=("ratio", "max")).reset_index().sort_values("KO"))
                head = f"### genome_code: {r.code}"
                if variant == "named":
                    head += f"  (organism: {r.organism_name})"
                lines.append(head + f"\n{len(agg)} significant energy-metabolism KOs:")
                lines += [f"{a.KO} | {panel.get(a.KO, '')} | {a.n} | {a.best:.2f}" for a in agg.itertuples()]
            (B2DIR / "batches" / f"{variant}_batch{b // BATCH:02d}.txt").write_text("\n".join(lines) + "\n")
    print("subset", len(idx), "positives per trait:", dict(zip(traits, Y[idx].sum(0).astype(int))))
    print("batches written:", sorted(p.name for p in (B2DIR / "batches").glob("*.txt")))


def score():
    lab, traits, Y = load_labelled()
    key = pd.read_csv(B2DIR / "subset_key.tsv", sep="\t")
    gi = {g: i for i, g in enumerate(lab.genome_id)}
    idx = np.array([gi[g] for g in key.genome_id])
    rows, summary = [], []
    P = {}
    for variant in ("genome", "named", "gpt_genome", "gpt_named"):
        files = sorted(f for f in (B2DIR / "out").glob(f"{variant}_batch*.json") if not f.name.endswith(".usage.json"))
        if not files:
            continue
        pred = {}
        for f in files:
            pred.update(json.loads(f.read_text()))
        M = np.full((len(lab), len(traits)), np.nan)
        for _, r in key.iterrows():
            d = pred.get(r.code, {})
            for j, t in enumerate(traits):
                if t in d:
                    M[gi[r.genome_id], j] = float(d[t])
        missing = int(np.isnan(M[idx]).sum())
        P[f"H_{variant}"] = M
        for j, t in enumerate(traits):
            for i in idx:
                rows.append(dict(arm=f"H_{variant}", genome_id=lab.genome_id[i], trait=t, y=Y[i, j], prob=M[i, j]))
        print(variant, "files", len(files), "missing predictions", missing)
    for arm in ("B2", "B2XL", "B2XLE_large", "B2E", "ENS", "null_taxNN"):
        P[arm] = oof_matrix(RES / "oof", arm, "genus", lab, traits)
    llm_arms = [a for a in P if a.startswith("H_")]
    for arm in llm_arms:  # each LLM arm is scored on the genomes it returned predictions for (paired refs)
        have = ~np.isnan(P[arm][idx]).any(axis=1)
        ii = idx[have]
        sub = lab.iloc[ii].reset_index(drop=True)
        Ys = Y[ii]
        for name in [arm, "B2", "B2XL", "B2XLE_large", "B2E", "ENS", "null_taxNN"]:
            aps = [safe_ap(Ys[:, j], P[name][ii, j]) for j in range(len(traits))]
            row = dict(llm_arm=arm, model=name, n_genomes=int(have.sum()), macro_auprc=np.nanmean(aps),
                       **{f"auprc_{t}": v for t, v in zip(traits, aps)})
            if name != arm:
                d, lo, hi = paired_macro_delta(sub, traits, Ys, P[arm][ii], P[name][ii], "genus")
                row.update(llm_minus_model=d, ci_low=lo, ci_high=hi)
            summary.append(row)
    s = pd.DataFrame(summary)
    s.to_csv(RES / "baseline2_llm_summary.csv", index=False, float_format="%.4f")
    pd.DataFrame(rows).to_csv(RES / "baseline2_llm_predictions.csv", index=False, float_format="%.3f")
    print(s[["llm_arm", "model", "n_genomes", "macro_auprc", "llm_minus_model", "ci_low", "ci_high"]].round(4).to_string())


if __name__ == "__main__":
    {"prepare": prepare, "score": score}[sys.argv[1]]()
