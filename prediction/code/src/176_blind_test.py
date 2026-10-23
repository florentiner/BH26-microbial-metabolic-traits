"""Pre-registration BL (docs/preregistration_BL_blind_test.md): audit the blind evidence, freeze the test labels, then
score them once with the frozen RedoxKO (and RedoxKO-LX).

audit : every tested_positive / tested_negative status needs a quote found verbatim (after whitespace and quote-mark
        normalization) in the abstract, the Europe PMC full text or the cited source page; otherwise not_tested.
        Statements that are not growth experiments (stimulated growth, medium components, fermentative reduction,
        genome-based statements) are set to not_tested, as in src/72_fix_evidence.py. The two catch-all categories are
        dropped. Writes tests.csv and its SHA-256 (the labels are frozen before any model output exists).
score : checks the frozen model hashes and the tests hash, predicts every genome from its RefSeq accession with
        src/predict.py (arms RF and B2), and computes the registered endpoints.
usage: python src/176_blind_test.py audit|score [donor]   (donor = pre-registration BL2, pool_blind_donor/)
Output: validation_new_organisms/pool_blind/{tests.csv,tests.sha256,audit_log.csv,predictions_RF.csv,
        predictions_B2.csv,results.json}
"""
import hashlib
import html
import json
import re
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
B = ROOT / "validation_new_organisms/pool_blind"
FROZEN = {"results/models/RF/model.joblib": "310d84ef335579bace4d5118af33c8e46a9fe2810327171742e2177227b6362a",
          "results/models/RF/metadata.json": "624b1ee156fae2a07796aa0b6508b6778d16a6b2fa64f371fd96d6f0131cce6e",
          "src/predict.py": "f0bb8b53df71da6a52b2350bab28926a3d7c4e09729355099b9093718559c7af"}
CATS = ["donor_sulfur_oxidation", "donor_ammonia_oxidation", "donor_nitrite_oxidation", "donor_hydrogen_oxidation",
        "donor_CO_oxidation", "donor_iron_oxidation", "acceptor_oxygen", "acceptor_nitrate", "acceptor_sulfate",
        "acceptor_sulfur", "acceptor_CO2_methanogen_acetogen", "acceptor_iron"]
NOT_GROWTH = re.compile(r"stimulat|genom|\bgenes?\b|encod|in the medium|medium contain|fermentative reduc|during fermentation",
                        re.I)
# amendment A2: a nitrate result needs growth coupling (src/72: an API-strip or kit test is a biochemical test)
# BL2 amendment A2: "stimulated/enhanced growth" is not a growth test even though it contains the word growth
STIMULATED = re.compile(r"stimulat|enhanc", re.I)
NITRATE_GROWTH = re.compile(r"\bgrow|growth\b|denitrif|electron acceptor", re.I)
SEED, NBOOT = 20261022, 2000


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def norm(s):
    s = html.unescape(re.sub(r"<[^>]+>", " ", s or ""))
    s = s.replace("’", "'").replace("‘", "'").replace("“", '"').replace("”", '"')
    s = s.replace("−", "-").replace("–", "-").replace(" ", " ").replace("\xa0", " ")
    return re.sub(r"\s+", " ", s).strip().lower()


def flat(s):
    """amendment A2: compare without any whitespace, because tag stripping splits e.g. AT8<sup>T</sup> into 'at8 t'"""
    return re.sub(r"[\s\u2009\u200a\u202f]+", "", s)


def get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (research script)"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "ignore")
        except Exception:
            time.sleep(3 * (i + 1))
    return ""


def audit():
    pool = {r["accession"]: r for r in json.load(open(B / POOL))}
    ev = [r for k in range(1, NEV + 1) for r in json.load(open(B / f"evidence_{k}.json"))]
    cache, rows, log = {}, [], []
    for rec in ev:
        acc = rec["refseq_accession"]
        p = pool[acc]
        texts = [norm(p.get("title", "")), norm(p.get("abstract", ""))]
        if p.get("pmcid"):
            u = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{p['pmcid']}/fullTextXML"
            cache.setdefault(u, norm(get(u)))
            texts.append(cache[u])
        for cat in CATS:
            c = rec["categories"].get(cat, {})
            st = c.get("status", "not_tested")
            if st not in ("tested_positive", "tested_negative"):
                continue
            q = norm(c.get("quote", "")).strip('"').strip("'").strip()
            q = re.sub(r"\.\.\.|…", " ", q)
            parts = [x.strip() for x in re.split(r"\s{2,}| \[.*?\] ", q) if len(x.strip()) > 15] or [q]
            src = c.get("source", "")
            # amendment A2: also read every Europe PMC full text named in the source or other_sources
            for u in sorted(set(re.findall(r"https?://www\.ebi\.ac\.uk/\S+?fullTextXML", " ".join([src, *rec.get("other_sources", [])])))):
                cache.setdefault(u, norm(get(u)))
                texts.append(cache[u])
            ft = [flat(t) for t in texts]
            found = all(any(flat(part) in t for t in ft) for part in parts)
            if not found and src.startswith("http") and "ncbi.nlm.nih.gov/pubmed" not in src:
                cache.setdefault(src, norm(get(src)))
                found = all(flat(part) in flat(cache[src]) for part in parts)
            reason = ""
            if not found:
                reason = "quote not found verbatim in an accessible source"
            elif DONOR and STIMULATED.search(q):
                reason = "stimulated or enhanced growth is not a growth test"
            elif NOT_GROWTH.search(q) and not re.search(r"\bgrow|growth\b", q):
                reason = "statement is not a growth experiment"
            elif cat == "acceptor_nitrate" and not NITRATE_GROWTH.search(q):
                reason = "nitrate result from a biochemical test without growth coupling"
            if reason:
                log.append(dict(accession=acc, category=cat, status=st, reason=reason, quote=c.get("quote", "")[:300]))
                continue
            rows.append(dict(accession=acc, species=rec["species"], paper=p.get("pmid"), tier=p.get("tier", 1), category=cat,
                             expected=int(st == "tested_positive"), quote=c.get("quote", ""), source=src))
    T = pd.DataFrame(rows)
    T.to_csv(B / "tests.csv", index=False)
    pd.DataFrame(log).to_csv(B / "audit_log.csv", index=False)
    (B / "tests.sha256").write_text(f"{sha(B / 'tests.csv')}  tests.csv\n")
    print(f"{len(ev)} species, {len(T)} tests kept ({int(T.expected.sum())} positive), {len(log)} demoted")
    print(T.groupby("category").expected.agg(["count", "sum"]).to_string())


def macro_ci(x):
    from sklearn.metrics import roc_auc_score

    def macro(y):
        v = [roc_auc_score(g.expected, g.score) for _, g in y.groupby("category") if g.expected.nunique() == 2]
        return float(np.mean(v)) if len(v) >= 2 else np.nan
    cats = [c for c, g in x.groupby("category") if g.expected.nunique() == 2]
    rng = np.random.default_rng(SEED)
    papers = x.paper.unique()
    idx = {q: np.where(x.paper.to_numpy() == q)[0] for q in papers}
    boots = [macro(x.iloc[np.concatenate([idx[q] for q in rng.choice(papers, len(papers))])]) for _ in range(NBOOT)] if len(papers) else []
    boots = np.array(boots, dtype=float)
    return {"tests": int(len(x)), "positives": int(x.expected.sum()), "categories": cats,
            "per_category_auroc": {c: float(roc_auc_score(g.expected, g.score)) for c, g in x.groupby("category")
                                   if g.expected.nunique() == 2},
            "macro_auroc": macro(x),
            "macro_auroc_ci": np.nanpercentile(boots, [2.5, 97.5]).tolist() if np.isfinite(boots).any() else None}


def predict(arm):
    out = B / f"predictions_{arm}.csv"
    if not out.exists():
        accs = sorted(pd.read_csv(B / "tests.csv").accession.unique())
        subprocess.run([sys.executable, str(ROOT / "src/predict.py"), *accs, "--arm", arm, "--pathways", "none",
                        "--out", str(out)], check=True)
    return pd.read_csv(out)


def score():
    for f, h in FROZEN.items():
        assert sha(ROOT / f) == h, f"{f} changed after registration"
    assert (B / "tests.sha256").read_text().split()[0] == sha(B / "tests.csv"), "tests changed after freezing"
    from sklearn.metrics import matthews_corrcoef, roc_auc_score
    T = pd.read_csv(B / "tests.csv")
    res = {"tests": int(len(T)), "positives": int(T.expected.sum()), "species": int(T.accession.nunique()),
           "papers": int(T.paper.nunique()), "per_category": T.groupby("category").expected.agg(["count", "sum"]).to_dict("index")}
    for arm in ("RF", "B2"):
        P = predict(arm)
        gcol = "genome" if "genome" in P else P.columns[0]
        s = P.set_index([gcol, "trait"])
        d = T.copy()
        d["score"] = [s.loc[(a, c), "score"] for a, c in zip(d.accession, d.category)]
        d["call"] = [s.loc[(a, c), "predicted"] for a, c in zip(d.accession, d.category)]
        d.to_csv(B / f"scored_{arm}.csv", index=False)

        def macro(x):
            v = [roc_auc_score(g.expected, g.score) for _, g in x.groupby("category") if g.expected.nunique() == 2]
            return float(np.mean(v)) if len(v) >= 2 else np.nan
        cats = [c for c, g in d.groupby("category") if g.expected.nunique() == 2]
        rng = np.random.default_rng(SEED)
        papers = d.paper.unique()
        idx = {p: np.where(d.paper.to_numpy() == p)[0] for p in papers}
        boots = []
        for _ in range(NBOOT):
            x = d.iloc[np.concatenate([idx[p] for p in rng.choice(papers, len(papers))])]
            try:
                boots.append(macro(x))
            except ValueError:
                pass
        # released yes/no calls only in the categories that RedoxKO's release rule passes (RedoxKO-LX has no rule of
        # its own in its metadata, so the same categories are used for it); predict.py outputs 0/1 for every category
        rel = d[d.category.isin(json.load(open(ROOT / "results/models/RF/metadata.json"))["reliable_traits"])]
        r = {"per_category_auroc": {c: float(roc_auc_score(g.expected, g.score)) for c, g in d.groupby("category")
                                    if g.expected.nunique() == 2},
             "macro_auroc": macro(d), "macro_auroc_ci": (np.nanpercentile(boots, [2.5, 97.5]).tolist()
                                if np.isfinite(boots).any() else None),
             "categories_in_auroc": cats,
             "released_calls": int(len(rel)),
             "released_accuracy": float((rel.expected == rel.call).mean()) if len(rel) else None,
             "released_mcc": float(matthews_corrcoef(rel.expected, rel.call)) if rel.expected.nunique() == 2 else None}
        res["RedoxKO" if arm == "RF" else "RedoxKO-LX"] = r
    if DONOR:   # BL2 primary endpoint: donor categories only; also each tier
        for arm, name in (("RF", "RedoxKO"), ("B2", "RedoxKO-LX")):
            d = pd.read_csv(B / f"scored_{arm}.csv")
            out = {}
            for label, x in (("donor_all_tiers", d[d.category.str.startswith("donor_")]),
                             ("donor_tier1", d[d.category.str.startswith("donor_") & (d.tier == 1)]),
                             ("donor_tier2", d[d.category.str.startswith("donor_") & (d.tier == 2)]),
                             ("all_tier1", d[d.tier == 1]), ("all_tier2", d[d.tier == 2])):
                out[label] = macro_ci(x)
            res[name]["bl2"] = out
        dn = T[T.category.str.startswith("donor_")]
        both = [c for c, g in dn.groupby("category") if g.expected.nunique() == 2]
        res["descriptive_only"] = bool(len(dn) < 20 or len(both) < 2)
    else:
        res["descriptive_only"] = bool(res["tests"] < 20 or len(res["RedoxKO"]["categories_in_auroc"]) < 2)
    (B / "results.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=lambda v: round(float(v), 3)))


DONOR = len(sys.argv) > 2 and sys.argv[2] == "donor"
if DONOR:
    B, POOL, NEV = ROOT / "validation_new_organisms/pool_blind_donor", "pool_all.json", 6
else:
    POOL, NEV = "pool.json", 4

if __name__ == "__main__":
    {"audit": audit, "score": score}[sys.argv[1]]()
