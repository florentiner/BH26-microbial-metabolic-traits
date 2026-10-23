"""X1 (docs/preregistration_additional_experiments.md): our saved models versus tools a reader could use instead,
on the primary 100 newly described species (255 blind growth tests).

Comparators: FAPROTAX 1.2.12 (taxonomy), the taxonomic nearest-group baseline (358 training genomes), the marker
rule (arm A) on the cached 852-KO profiles, and BacDive-AI's aerobic model (acceptor_oxygen only; its scores are
produced by src/87b_bacdive_ai.py in a separate scikit-learn 1.4.0 environment and merged here if present).
No model is refitted. Writes results/additional/x1_*.csv.

usage: python src/87_tool_comparison.py
"""
import fnmatch
import hashlib
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / "validation_new_organisms"
OUT = ROOT / "results/additional"
OUT.mkdir(parents=True, exist_ok=True)
FAPRO = ROOT / "data/additional/faprotax/FAPROTAX_1.2.12/FAPROTAX.txt"
PREREG = ROOT / "docs/preregistration_additional_experiments.md"
SEED = 20261014
TAX = ["domain", "phylum", "class", "order", "family", "genus"]

FAPRO_MAP = {
    "acceptor_oxygen": ["aerobic_chemoheterotrophy", "aerobic_ammonia_oxidation", "aerobic_nitrite_oxidation",
                        "nitrification"],
    "acceptor_nitrate": ["nitrate_respiration", "nitrate_reduction", "nitrate_denitrification", "denitrification",
                         "nitrogen_respiration"],
    "acceptor_sulfate": ["sulfate_respiration"],
    "acceptor_sulfur": ["sulfur_respiration", "thiosulfate_respiration", "sulfite_respiration"],
    "acceptor_iron": ["iron_respiration"],
    "acceptor_CO2_methanogen_acetogen": ["@methanogenesis*", "@hydrogenotrophic_methanogenesis*", "acetogenesis",
                                         "reductive_acetogenesis"],  # last one: deviation D1 (group renamed)
    "donor_sulfur_oxidation": ["dark_oxidation_of_sulfur_compounds", "dark_sulfur_oxidation",
                               "dark_sulfide_oxidation", "dark_sulfite_oxidation", "dark_thiosulfate_oxidation"],
    "donor_hydrogen_oxidation": ["dark_hydrogen_oxidation"],
    "donor_iron_oxidation": ["dark_iron_oxidation"],
    "donor_ammonia_oxidation": ["aerobic_ammonia_oxidation"],
    "donor_nitrite_oxidation": ["aerobic_nitrite_oxidation"],
    "donor_CO_oxidation": [],
}


# ---------------------------------------------------------------- FAPROTAX (data file only; its script is not run)
def parse_faprotax(path):
    groups, cur = {}, None
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        if not raw[0].isspace() and not line.startswith("*") and ":" in line and "\t" in raw \
                and not line.startswith(("add_group", "subtract_group", "intersect_group")):
            cur = line.split()[0]
            groups[cur] = {"patterns": [], "add": [], "sub": []}
            continue
        if cur is None:
            continue
        tok = line.strip().split()[0]
        if tok.startswith("add_group:"):
            groups[cur]["add"].append(tok.split(":", 1)[1])
        elif tok.startswith("subtract_group:"):
            groups[cur]["sub"].append(tok.split(":", 1)[1])
        elif "*" in tok or tok[:1].isalpha():
            groups[cur]["patterns"].append(tok)
    return groups


def faprotax_members(groups, name, lineage, seen=()):
    """True if the lineage string belongs to group `name` (patterns, then add_group, minus subtract_group)."""
    if name in seen or name not in groups:
        return False
    g = groups[name]
    hit = any(fnmatch.fnmatchcase(lineage, p) for p in g["patterns"]) or \
        any(faprotax_members(groups, a, lineage, seen + (name,)) for a in g["add"])
    if hit and any(faprotax_members(groups, s, lineage, seen + (name,)) for s in g["sub"]):
        return False
    return hit


def resolve_map(groups):
    resolved, missing = {}, []
    for cat, names in FAPRO_MAP.items():
        out = []
        for n in names:
            if n.startswith("@"):
                out += [g for g in groups if fnmatch.fnmatchcase(g, n[1:])]
            elif n in groups:
                out.append(n)
            else:
                missing.append(n)
        resolved[cat] = sorted(set(out))
    return resolved, sorted(set(missing))


# ---------------------------------------------------------------- NCBI lineage (cached)
def ncbi_json(url, tries=6):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except Exception:
            time.sleep(2 ** k)
    raise RuntimeError(url)


def lineages(accs, cache):
    cache.parent.mkdir(parents=True, exist_ok=True)
    have = pd.read_csv(cache, sep="\t") if cache.exists() else pd.DataFrame(columns=["accession"])
    todo = [a for a in accs if a not in set(have.accession)]
    rows = []
    for i in range(0, len(todo), 50):
        chunk = todo[i:i + 50]
        rep = ncbi_json("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
                        + "%2C".join(chunk) + "/dataset_report?page_size=100")
        tid = {r["accession"]: r["organism"]["tax_id"] for r in rep.get("reports", [])}
        ids = sorted(set(tid.values()))
        cls = {}
        for j in range(0, len(ids), 50):
            tx = ncbi_json("https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/"
                           + "%2C".join(map(str, ids[j:j + 50])) + "/dataset_report")
            for r in tx.get("reports", []):
                t = r["taxonomy"]
                c = t.get("classification", {})
                cls[t["tax_id"]] = {
                    "domain": c.get("domain", {}).get("name") or c.get("superkingdom", {}).get("name"),
                    "phylum": c.get("phylum", {}).get("name"), "class": c.get("class", {}).get("name"),
                    "order": c.get("order", {}).get("name"), "family": c.get("family", {}).get("name"),
                    "genus": c.get("genus", {}).get("name"), "species": c.get("species", {}).get("name"),
                    "organism_name": t.get("current_scientific_name", {}).get("name")}
        for a in chunk:
            if a in tid:
                rows.append({"accession": a, "tax_id": tid[a], **cls.get(tid[a], {})})
    out = pd.concat([have, pd.DataFrame(rows)], ignore_index=True)
    out.to_csv(cache, sep="\t", index=False)
    return out.set_index("accession").loc[accs]


# ---------------------------------------------------------------- metrics
def point(y, c):
    y, c = np.asarray(y, int), np.asarray(c, int)
    tp, tn = int(((y == 1) & (c == 1)).sum()), int(((y == 0) & (c == 0)).sum())
    fp, fn = int(((y == 0) & (c == 1)).sum()), int(((y == 1) & (c == 0)).sum())
    return dict(n=len(y), accuracy=(tp + tn) / len(y), mcc=matthews_corrcoef(y, c) if len(set(c)) > 1 else 0.0,
                sensitivity=tp / max(tp + fn, 1), specificity=tn / max(tn + fp, 1))


def macro_auroc(df, col):
    v = [roc_auc_score(g.expected, g[col]) for _, g in df.groupby("category")
         if g.expected.nunique() == 2 and g[col].notna().all()]
    return float(np.mean(v)) if v else np.nan


def cluster_boot(df, fn, n=1000, seed=SEED):
    rng = np.random.default_rng(seed)
    papers = df.paper.unique()
    idx = {p: np.where(df.paper.to_numpy() == p)[0] for p in papers}
    out = []
    for _ in range(n):
        take = np.concatenate([idx[p] for p in rng.choice(papers, len(papers))])
        try:
            out.append(fn(df.iloc[take]))
        except ValueError:
            continue
    return np.nanpercentile(out, [2.5, 97.5])


def main():
    prereg_sha = hashlib.sha256(PREREG.read_bytes()).hexdigest()
    b2 = pd.read_csv(VAL / "scored178_B2_balanced100.csv")
    fu = pd.read_csv(VAL / "scored178_B2XLE_large_balanced100.csv")
    calls = b2[b2.blind & ~b2.category.isin(["donor_other", "acceptor_other"])].copy()
    calls = calls.rename(columns={"score": "score_B2", "call": "call_B2"})
    fu = fu.rename(columns={"score": "score_B2XLE_large", "call": "call_B2XLE_large"})
    calls = calls.merge(fu[["accession", "category", "score_B2XLE_large", "call_B2XLE_large"]],
                        on=["accession", "category"], how="left")
    assert len(calls) == 255, len(calls)
    accs = sorted(calls.accession.unique())

    lin = lineages(accs, ROOT / "data/additional/lineage_new_species.tsv")

    # FAPROTAX
    groups = parse_faprotax(FAPRO)
    cmap, missing = resolve_map(groups)
    fap = {}
    for a in accs:
        r = lin.loc[a]
        names = [r.get(k) for k in ["domain", "phylum", "class", "order", "family", "genus", "species"]]
        s = ";".join([str(x) for x in names if isinstance(x, str)])
        assigned = sorted(g for g in groups if faprotax_members(groups, g, s))
        fap[a] = {"lineage": s, "groups": assigned}
    calls["faprotax_covered"] = calls.accession.map(lambda a: len(fap[a]["groups"]) > 0)
    calls["call_FAPROTAX"] = [int(any(g in fap[a]["groups"] for g in cmap.get(c, [])))
                              for a, c in zip(calls.accession, calls.category)]

    # taxonomic nearest-group baseline from the 358 training genomes
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    lab = lab[lab.labelled == True]  # noqa: E712
    def taxnn(a, cat):
        r = lin.loc[a]
        for rank in reversed(TAX):
            v = r.get(rank)
            if isinstance(v, str) and (lab[rank] == v).any():
                return float(lab.loc[lab[rank] == v, cat].mean()), rank
        return float(lab[cat].mean()), "global"
    tn = [taxnn(a, c) for a, c in zip(calls.accession, calls.category)]
    calls["score_taxNN"] = [t[0] for t in tn]
    calls["taxNN_rank"] = [t[1] for t in tn]
    prev = lab[[c for c in calls.category.unique()]].mean()
    calls["call_taxNN"] = (calls.score_taxNN >= calls.category.map(prev)).astype(int)

    # marker rule (arm A) on cached KO profiles
    ko = pd.read_parquet(VAL / "ko_panel_validation.parquet")
    ko = (ko > 0).astype(float)
    A = joblib.load(ROOT / "results/models/A/model.joblib")
    def rule(a, cat):
        b = A[cat]
        if a not in ko.index:
            return np.nan
        if not b["groups"]:  # no marker group for this category: score 0, as src/predict.py does
            return 0.0
        cols = [[k for k in grp if k in ko.columns] for grp in b["groups"].values()]
        return float(max(ko.loc[a, c].mean() if c else 0.0 for c in cols))
    calls["score_A"] = [rule(a, c) for a, c in zip(calls.accession, calls.category)]
    calls["call_A"] = [int(s >= A[c]["threshold"]) if pd.notna(s) else np.nan
                       for s, c in zip(calls.score_A, calls.category)]

    # BacDive-AI (oxygen only), produced by src/87b_bacdive_ai.py
    bd = OUT / "x1_bacdive_ai_aerobic.csv"
    if bd.exists():
        d = pd.read_csv(bd)
        m = d.set_index("accession")
        calls["score_BacDiveAI"] = [m.at[a, "p_aerobic"] if (c == "acceptor_oxygen" and a in m.index) else np.nan
                                    for a, c in zip(calls.accession, calls.category)]
        calls["call_BacDiveAI"] = (calls.score_BacDiveAI >= 0.5).astype(float).where(calls.score_BacDiveAI.notna())

    # METABOLIC (src/95) and DRAM (src/94), added by Addenda A and B, merged if present
    for tool, f in [("METABOLIC", "x1_metabolic_calls.csv"), ("DRAM", "x1_dram_calls.csv"),
                    ("MicrobeDecoder", "x1_microbe_decoder_calls.csv")]:
        if (OUT / f).exists():
            m = pd.read_csv(OUT / f)
            keep = [c for c in m.columns if c in (f"call_{tool}", f"score_{tool}")]
            calls = calls.drop(columns=keep, errors="ignore").merge(
                m[["accession", "category"] + keep], on=["accession", "category"], how="left")

    calls.to_csv(OUT / "x1_calls.csv", index=False)

    # summaries
    rows = []
    tools = ["B2", "B2XLE_large", "A", "taxNN", "FAPROTAX"] + [t for t in ["METABOLIC", "DRAM", "MicrobeDecoder", "BacDiveAI"]
                                                                if f"call_{t}" in calls]
    for subset, df in [("all 255 calls", calls), ("FAPROTAX-covered species", calls[calls.faprotax_covered]),
                       ("acceptor_oxygen", calls[calls.category == "acceptor_oxygen"]),
                       ("Microbe Decoder categories (oxygen, nitrate, methanogen)",
                        calls[calls.category.isin(["acceptor_oxygen", "acceptor_nitrate",
                                                   "acceptor_CO2_methanogen_acetogen"])])]:
        for t in tools:
            d = df[df[f"call_{t}"].notna()]
            if d.empty:
                continue
            r = dict(subset=subset, tool=t, organisms=d.accession.nunique(), **point(d.expected, d[f"call_{t}"]))
            sc = f"score_{t}"
            r["macro_auroc"] = macro_auroc(d, sc) if sc in d else np.nan
            if t != "B2" and len(d) == len(df):
                f = lambda x: (matthews_corrcoef(x.expected, x.call_B2)
                               - matthews_corrcoef(x.expected, x[f"call_{t}"]))
                lo, hi = cluster_boot(d, f)
                r["mcc_B2_minus_tool"] = f(d)
                r["mcc_B2_minus_tool_ci"] = f"[{lo:+.3f}, {hi:+.3f}]"
                if sc in d and d[sc].notna().all():
                    g = lambda x: macro_auroc(x, "score_B2") - macro_auroc(x, sc)
                    lo, hi = cluster_boot(d, g)
                    r["auroc_B2_minus_tool"] = g(d)
                    r["auroc_B2_minus_tool_ci"] = f"[{lo:+.3f}, {hi:+.3f}]"
            rows.append(r)
    summ = pd.DataFrame(rows)
    summ.to_csv(OUT / "x1_summary.csv", index=False)
    cov = pd.DataFrame([{"accession": a, "faprotax_groups": ";".join(v["groups"]), "lineage": v["lineage"]}
                        for a, v in fap.items()])
    cov.to_csv(OUT / "x1_faprotax_assignments.csv", index=False)
    per_cat = []
    for c, d in calls.groupby("category"):
        for t in tools:
            dd = d[d[f"call_{t}"].notna()]
            if len(dd):
                per_cat.append(dict(category=c, tool=t, **point(dd.expected, dd[f"call_{t}"])))
    pd.DataFrame(per_cat).to_csv(OUT / "x1_per_category.csv", index=False)
    meta = dict(prereg_sha256=prereg_sha, faprotax="1.2.12", faprotax_groups_missing=missing,
                faprotax_category_groups=cmap, faprotax_coverage=float(np.mean([len(v["groups"]) > 0
                                                                                for v in fap.values()])),
                taxNN_rank_counts=calls.drop_duplicates("accession").taxNN_rank.value_counts().to_dict())
    (OUT / "x1_meta.json").write_text(json.dumps(meta, indent=2))
    pd.set_option("display.width", 200)
    print(summ.round(3).to_string(index=False))
    print(json.dumps({k: meta[k] for k in ["faprotax_groups_missing", "faprotax_coverage", "taxNN_rank_counts"]}))


if __name__ == "__main__":
    sys.exit(main())
