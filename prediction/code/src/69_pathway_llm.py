"""LLM baseline for the pathway task — the donor/acceptor protocol of src/21_baseline2_llm.py, 20 pathway traits.

Same evidence as the released pathway model sees: significant KofamScan hits from the RDF Portal BH26 Microbes
endpoint (TogoMCP), restricted to the 852-KO energy panel so the prompt stays readable. Two arms:
  P_H_genome : anonymised genome codes, evidence only
  P_H_named  : the same evidence plus the organism name (tests literature recall rather than genome reading)
Claude sub-agents answer the batches; no paid API is used.

usage:
  python src/69_pathway_llm.py prepare   # subset + batches -> data/pathway_llm/batches/
  python src/69_pathway_llm.py score     # data/pathway_llm/out/*.json -> results/pathway_llm_*.csv
"""
import json
import sys

import numpy as np
import pandas as pd

from pathway_common import RES, load_pathways
from redox_common import ROOT, paired_macro_delta, safe_ap

OUTDIR = ROOT / "data/pathway_llm"
SEED = 20260921
N_SUBSET = 96
BATCH = 12
TRAITS = {
    "denitrification": "Denitrification: nitrate/nitrite reduced to N2O or N2",
    "dnra": "Dissimilatory nitrate reduction to ammonium (DNRA)",
    "nitrate_reduction": "Dissimilatory nitrate reduction (any route, includes denitrification and DNRA)",
    "nitrogen_fixation": "N2 fixation to ammonia (nitrogenase)",
    "methanogenesis": "Methanogenesis (any substrate)",
    "acetogenesis": "Acetogenesis: reductive acetyl-CoA (Wood-Ljungdahl) pathway to acetate",
    "sulfate_reduction": "Dissimilatory sulfate reduction",
    "sulfur_reduction": "Dissimilatory reduction of elemental sulfur, thiosulfate or polysulfide",
    "sulfur_oxidation": "Oxidation of reduced sulfur compounds for energy",
    "iron_reduction": "Dissimilatory Fe(III) reduction",
    "iron_oxidation": "Fe(II) oxidation for energy",
    "manganese_reduction": "Dissimilatory Mn(IV) reduction",
    "arsenate_reduction": "Dissimilatory arsenate reduction",
    "fumarate_reduction": "Fumarate respiration (fumarate as terminal electron acceptor)",
    "hydrogen_oxidation": "H2 oxidation as an energy source",
    "fermentation": "Fermentative energy metabolism (substrate-level phosphorylation, no external acceptor)",
    "aerobic_chemoheterotrophy": "Aerobic growth on organic carbon",
    "photoautotrophy": "Photoautotrophy (light energy, CO2 as carbon source)",
    "photoheterotrophy": "Photoheterotrophy (light energy, organic carbon)",
    "methylotrophy": "Growth on C1 compounds such as methanol, methylamines or methane",
}
PROMPT = """You are an expert microbial physiologist. For each genome below you get ONLY genome annotation evidence:
significant KofamScan hits (KEGG Orthology) in energy metabolism, retrieved from the RDF Portal BH26 Microbes
database (the data behind TogoMCP). Each line: KO | KEGG definition | number of proteins | best score/threshold.

Task: for every genome and every trait, give the probability (0.00-1.00) that the organism uses that energy
metabolism pathway, as it would be recorded in a curated microbial trait database. Reason from the gene content
(complete pathways and complexes, marker specificity, ambiguous KOs such as narGH vs nxrAB, sdh vs frd, or
dsrAB used in both oxidative and reductive directions). Do NOT use any tool, file, web search or outside lookup:
judge from the evidence given. {anon}

Traits (key: meaning):
{traits}

Return ONLY a JSON object: {{"<genome_code>": {{"<trait_key>": probability, ...}}, ...}} with all 20 keys for every
genome. Write it to the output path given by the orchestrator.
"""


def pick_subset(lab, traits, Y):
    """96 genomes: at least 5 positives for the rarest traits first, then a random fill."""
    rng = np.random.default_rng(SEED)
    chosen = []
    for j in np.argsort(Y.sum(0)):
        have = int(Y[chosen, j].sum()) if chosen else 0
        pos = [i for i in np.where(Y[:, j] == 1)[0] if i not in chosen]
        need = max(0, 5 - have)
        if need and pos:
            chosen += list(rng.choice(pos, min(need, len(pos)), replace=False))
    rest = [i for i in range(len(lab)) if i not in chosen]
    chosen += list(rng.choice(rest, max(0, N_SUBSET - len(chosen)), replace=False))
    return sorted(int(i) for i in chosen[:N_SUBSET])


def prepare():
    lab, traits, Y = load_pathways()
    idx = pick_subset(lab, traits, Y)
    panel = pd.read_csv(ROOT / "data/processed/ko_panel.tsv", sep="\t").set_index("KO").definition
    togo = pd.read_parquet(ROOT / "data/pathways/togo_kofam.parquet")
    ids = set(lab.genome_id[idx])
    togo = togo[togo.genome_id.isin(ids) & togo.significant & togo.KO.isin(panel.index)]
    rng = np.random.default_rng(SEED + 1)
    codes = {lab.genome_id[i]: f"P{c:03d}" for i, c in zip(idx, rng.permutation(len(idx)))}
    (OUTDIR / "batches").mkdir(parents=True, exist_ok=True)
    (OUTDIR / "out").mkdir(parents=True, exist_ok=True)
    key = lab.loc[idx, ["genome_id", "species"]].assign(code=lambda d: d.genome_id.map(codes))
    key.to_csv(OUTDIR / "subset_key.tsv", sep="\t", index=False)     # never shown to the anonymised arm
    tlist = "\n".join(f"- {k}: {v}" for k, v in TRAITS.items())
    for variant in ("genome", "named"):
        ordered = key.sort_values("code").reset_index(drop=True)
        anon = ("Genomes are anonymised; do not guess identities." if variant == "genome"
                else "The organism name is given with each genome.")
        for b in range(0, len(ordered), BATCH):
            lines = [PROMPT.format(traits=tlist, anon=anon)]
            for _, r in ordered.iloc[b:b + BATCH].iterrows():
                g = togo[togo.genome_id == r.genome_id]
                agg = (g.assign(ratio=g.score / g.threshold).groupby("KO")
                       .agg(n=("protein", "nunique"), best=("ratio", "max")).reset_index().sort_values("KO"))
                head = f"### genome_code: {r.code}"
                if variant == "named":
                    head += f"  (organism: {r.species})"
                lines.append(head + f"\n{len(agg)} significant energy-metabolism KOs:")
                lines += [f"{a.KO} | {panel.get(a.KO, '')} | {a.n} | {a.best:.2f}" for a in agg.itertuples()]
            (OUTDIR / "batches" / f"{variant}_batch{b // BATCH:02d}.txt").write_text("\n".join(lines) + "\n")
    print("subset", len(idx), "positives per trait:", {t: int(n) for t, n in zip(traits, Y[idx].sum(0))})
    print("batches:", sorted(p.name for p in (OUTDIR / "batches").glob("*.txt")))


def score():
    lab, traits, Y = load_pathways()
    key = pd.read_csv(OUTDIR / "subset_key.tsv", sep="\t")
    gi = {g: i for i, g in enumerate(lab.genome_id)}
    idx = np.array([gi[g] for g in key.genome_id])
    P, rows = {}, []
    for variant in ("genome", "named"):
        files = sorted((OUTDIR / "out").glob(f"{variant}_batch*.json"))
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
        P[f"P_H_{variant}"] = M
        print(f"P_H_{variant}: {int(np.isfinite(M[idx]).sum())}/{len(idx) * len(traits)} answers")
    # the model arms, restricted to the same genomes
    for arm in ("PB3", "PB2", "PB1", "PE", "ET", "RF", "PA", "null_taxNN"):
        f = RES / f"pathways/oof/{arm}.parquet"
        if not f.exists():
            f = RES / f"trees_xai/pw/oof/{arm}.parquet"
        if not f.exists() and arm == "null_taxNN":
            f = RES / "pathways/oof/nulls.parquet"
        if not f.exists():
            continue
        d = pd.read_parquet(f)
        d = d[(d.level == "genus") & (d.arm == arm)] if "arm" in d.columns else d[d.level == "genus"]
        if d.empty:
            continue
        M = np.full((len(lab), len(traits)), np.nan)
        tj = {t: j for j, t in enumerate(traits)}
        for t, g in d.groupby("trait"):
            M[g.genome_id.map(gi).to_numpy(), tj[t]] = g.prob.to_numpy()
        P[arm] = M
    for a, M in P.items():
        auprc = [safe_ap(Y[idx, j], M[idx, j]) for j in range(len(traits)) if np.isfinite(M[idx, j]).all()
                 and 0 < Y[idx, j].sum() < len(idx)]
        rows.append(dict(arm=a, n_genomes=len(idx), n_traits=len(auprc), macro_auprc=float(np.mean(auprc))))
    out = pd.DataFrame(rows)
    comp = []
    for llm in [a for a in P if a.startswith("P_H_")]:
        for arm in [a for a in P if not a.startswith("P_H_")]:
            sub = lab.loc[idx].reset_index(drop=True)
            d, lo, hi = paired_macro_delta(sub, traits, Y[idx], P[llm][idx], P[arm][idx], "genus")
            comp.append(dict(llm_arm=llm, model=arm, delta=d, ci_low=lo, ci_high=hi))
    out.to_csv(RES / "pathway_llm_summary.csv", index=False, float_format="%.4f")
    pd.DataFrame(comp).to_csv(RES / "pathway_llm_contrasts.csv", index=False, float_format="%.4f")
    print(out.to_string(index=False))
    print(pd.DataFrame(comp).to_string(index=False))


if __name__ == "__main__":
    {"prepare": prepare, "score": score}[sys.argv[1]]()
