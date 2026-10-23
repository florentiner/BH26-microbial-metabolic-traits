"""D (docs/preregistration_label_recuration_D.md): re-curate the training labels of the 358 genomes (355 species)
with Claude Haiku 5.5 extraction + verbatim-quote check + Qwen 3.8 Flash confirmation. Hard stop at $6.00.

usage: python src/101_label_recuration.py fetch | extract | labels
"""
import json
import random
import re
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
pilot = __import__("100_llm_pilot")
W = ROOT / "data/additional/label_recuration"
AUDIT = Path("~/work/hack_res_fedotmas/final/phenotype_data_GRADED_OPEN_ACCESS.csv")
SEED, MAX_CHARS, BUDGET = 20261015, 40000, 6.00
EXTRACT, CHECK = "anthropic/claude-haiku-5.5", "qwen/qwen3.8-flash"
CATS = list(pilot.DESC)
LOCK = threading.Lock()
USE_CHECK = "--with-check" in sys.argv  # deviation D12: off by default
SPENT = [0.0]


def epmc(query):
    js = pilot.get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=" +
                   urllib_quote(query) + "&resultType=core&format=json&pageSize=5")
    return json.loads(js).get("resultList", {}).get("result", []) if js else []


def urllib_quote(s):
    import urllib.parse
    return urllib.parse.quote(s)


def fulltext(pmcid):
    xml = pilot.get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML")
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", xml)) if xml else ""


def fetch():
    W.mkdir(parents=True, exist_ok=True)
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    species = sorted(lab[lab.labelled == True].species.dropna().unique())  # noqa: E712
    p = pd.read_csv(AUDIT, low_memory=False)
    p["binomial"] = p.organism.astype(str).str.replace(r"\(.*?\)", "", regex=True).str.split().str[:2].str.join(" ")
    cited = {}
    for _, r in p.iterrows():
        ids = re.split(r";\s*", f"{r.get('electron_donors_proof_dois')};{r.get('electron_acceptors_proof_dois')}")
        cited.setdefault(r.binomial, []).extend(i for i in ids if i.isdigit())
    out_f = W / "evidence.jsonl"
    done = set()
    if out_f.exists():
        done = {json.loads(l)["species"] for l in out_f.read_text().splitlines()}

    def one(sp):
        parts, used, abstracts = [], [], []
        for pm in list(dict.fromkeys(cited.get(sp, [])))[:6]:
            res = epmc(f"EXT_ID:{pm} AND SRC:MED")
            if not res:
                continue
            r = res[0]
            if r.get("isOpenAccess") == "Y" and r.get("pmcid"):
                t = fulltext(r["pmcid"])
                if t:
                    parts.append(t); used.append(f"PMID:{pm} full text")
                    continue
            if r.get("abstractText"):
                abstracts.append(re.sub(r"<[^>]+>", " ", r["abstractText"])); used.append(f"PMID:{pm} abstract")
        if sum(map(len, parts)) < MAX_CHARS:
            for r in epmc(f'"{sp}" AND ("sp. nov." OR "gen. nov.") AND OPEN_ACCESS:y'):
                if r.get("pmcid") and f"PMID:{r.get('pmid')} full text" not in used:
                    t = fulltext(r["pmcid"])
                    if t:
                        parts.append(t); used.append(f"{r['pmcid']} description full text")
                        break
        text = (" ".join(parts + abstracts))[:MAX_CHARS]
        return {"species": sp, "sources": used, "text": text}

    todo = [s for s in species if s not in done]
    with ThreadPoolExecutor(6) as ex:
        for i, rec in enumerate(ex.map(one, todo), 1):
            with open(out_f, "a") as f:
                f.write(json.dumps(rec) + "\n")
            if i % 25 == 0:
                print(i, "/", len(todo), flush=True)
    recs = [json.loads(l) for l in out_f.read_text().splitlines()]
    print("species", len(recs), "with text", sum(bool(r["text"]) for r in recs))


def norm(s):
    return re.sub(r"\s+", " ", s or "").strip().lower()


def check_prompt(sp, text, claims):
    lines = "\n".join(f'- "{c}" ({pilot.DESC[c]}): {v["status"]}; quote: "{v["quote"]}"' for c, v in claims.items())
    return (f"Organism: {sp}\nA first reader extracted the claims below from the paper text. For each, answer "
            f"\"confirm\" only if the text supports it for this organism as an experimental result (growth test, "
            f"activity assay); otherwise \"reject\".\n{lines}\n\nAnswer with JSON only: "
            f"{{\"<category>\": \"confirm\" | \"reject\", ...}}\n\nPAPER TEXT:\n{text}")


def extract():
    key = Path.home().joinpath(".openrouter_key").read_text().strip()
    recs = [json.loads(l) for l in (W / "evidence.jsonl").read_text().splitlines()]
    recs = [r for r in recs if r["text"]]
    random.Random(SEED).shuffle(recs)
    out_f = W / "curation.jsonl"
    done = set()
    if out_f.exists():
        for l in out_f.read_text().splitlines():
            d = json.loads(l)
            done.add(d["species"]); SPENT[0] += d["cost"]

    def one(r):
        with LOCK:
            if SPENT[0] >= BUDGET:
                return None
        it = {"organism": r["species"], "text": r["text"]}
        try:
            a = pilot.call(EXTRACT, pilot.prompt(it), key)
        except Exception as e:  # noqa: BLE001
            print("extract error", r["species"], e, flush=True)
            return None
        cost = (a.get("usage") or {}).get("cost", 0) or 0
        try:
            ans = json.loads(re.search(r"\{.*\}", a["choices"][0]["message"]["content"], re.S).group(0))
        except Exception:  # noqa: BLE001
            ans = {}
        claims, quote_fail = {}, 0
        for c in CATS:
            v = ans.get(c) if isinstance(ans.get(c), dict) else {}
            st, q = v.get("status", "not_tested"), v.get("quote", "") or ""
            if st in ("tested_positive", "tested_negative"):
                if q and norm(q) in norm(r["text"]):
                    claims[c] = {"status": st, "quote": q}
                else:
                    quote_fail += 1
        verdict = {}
        if claims and USE_CHECK:
            try:
                b = pilot.call(CHECK, check_prompt(r["species"], r["text"], claims), key)
                cost += (b.get("usage") or {}).get("cost", 0) or 0
                verdict = json.loads(re.search(r"\{.*\}", b["choices"][0]["message"]["content"], re.S).group(0))
            except Exception as e:  # noqa: BLE001
                print("check error", r["species"], e, flush=True)
        if USE_CHECK:
            final = {c: (claims[c]["status"] if c in claims and str(verdict.get(c, "")).lower() == "confirm"
                         else "unknown") for c in CATS}
        else:  # deviation D12: Haiku claims with a verbatim quote are final
            final = {c: (claims[c]["status"] if c in claims else "unknown") for c in CATS}
        rec = {"species": r["species"], "sources": r["sources"], "claims": claims, "verdict": verdict,
               "final": final, "quote_failures": quote_fail, "cost": cost}
        with LOCK:
            SPENT[0] += cost
            with open(out_f, "a") as f:
                f.write(json.dumps(rec) + "\n")
            print(f"{r['species']}: {sum(v != 'unknown' for v in final.values())} cells kept, "
                  f"${cost:.4f}, total ${SPENT[0]:.3f}", flush=True)
        return rec

    todo = [r for r in recs if r["species"] not in done]
    with ThreadPoolExecutor(8) as ex:
        list(ex.map(one, todo))
    print(f"done, spent ${SPENT[0]:.3f}")


def labels():
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    cur = {json.loads(l)["species"]: json.loads(l)["final"] for l in (W / "curation.jsonl").read_text().splitlines()}
    new = lab.copy()
    stats = {"new_positive": 0, "tested_negative": 0, "unknown": 0, "kept_positive": 0}
    for i, r in lab[lab.labelled == True].iterrows():  # noqa: E712
        f = cur.get(r.species, {})
        for c in CATS:
            if r[c] == 1:
                stats["kept_positive"] += 1
                continue
            s = f.get(c, "unknown")
            if s == "tested_positive":
                new.at[i, c] = 1.0; stats["new_positive"] += 1
            elif s == "tested_negative":
                new.at[i, c] = 0.0; stats["tested_negative"] += 1
            else:
                new.at[i, c] = np.nan; stats["unknown"] += 1
    new.to_csv(ROOT / "data/processed/labels_curated_D.tsv", sep="\t", index=False)
    print(stats, "species curated", len(cur))


if __name__ == "__main__":
    {"fetch": fetch, "extract": extract, "labels": labels}[sys.argv[1]]()
