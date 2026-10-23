"""Pilot for option D: which LLM can turn papers into donor/acceptor labels (tested positive / tested negative /
not tested) as well as the reference? 30 training organisms (seed 20261015) whose audited donor/acceptor values
cite open-access papers; texts from Europe PMC. Gold = the audited table (hack_res_fedotmas/final, open-access
version) mapped to the 14 categories with src/92_label_audit_impact.py's compound map; reported-not-used -> negative.
Hard spending stop: $1 total, using the cost OpenRouter reports for every call. Key read from ~/.openrouter_key.

usage: python src/100_llm_pilot.py fetch | run | score
"""
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
lab92 = __import__("92_label_audit_impact")
W = ROOT / "data/additional/llm_pilot"
OUT = ROOT / "results/additional/llm_pilot"
AUDIT = Path("~/work/hack_res_fedotmas/final/phenotype_data_GRADED_OPEN_ACCESS.csv")
SEED, N_ORG, MAX_CHARS, BUDGET = 20261015, 30, 40000, 1.00
CATS = lab92.CATS
MODELS = ["google/gemini-3.8-flash:batch", "anthropic/claude-haiku-5.5", "qwen/qwen3.8-flash", "qwen/qwen3.7-flash"]
FALLBACK = {"google/gemini-3.8-flash:batch": "google/gemini-3.8-flash"}
DESC = {
    "donor_sulfur_oxidation": "uses reduced sulfur compounds (sulfide, S0, thiosulfate, tetrathionate, sulfite) as electron donor",
    "donor_ammonia_oxidation": "oxidizes ammonia/ammonium for energy", "donor_nitrite_oxidation": "oxidizes nitrite for energy",
    "donor_hydrogen_oxidation": "uses H2 as electron donor", "donor_CO_oxidation": "uses carbon monoxide as electron donor",
    "donor_iron_oxidation": "oxidizes Fe(II) for energy", "acceptor_oxygen": "respires O2 (aerobic growth with O2 as acceptor)",
    "acceptor_nitrate": "reduces nitrate as electron acceptor", "acceptor_sulfate": "reduces sulfate as electron acceptor",
    "acceptor_sulfur": "reduces S0, thiosulfate or sulfite as electron acceptor",
    "acceptor_CO2_methanogen_acetogen": "reduces CO2 to methane or acetate (methanogen or acetogen)",
    "acceptor_iron": "reduces Fe(III) as electron acceptor"}


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:  # noqa: BLE001
            time.sleep(3 * (i + 1))
    return None


def gold_row(r):
    pos = lab92.categorise(r.get("electron_donors"), lab92.DONOR, "donor_other") | \
        lab92.categorise(r.get("electron_acceptors"), lab92.ACCEPTOR, "acceptor_other")
    neg_d = re.sub(r"\[[^\]]*\]", "", str(r.get("reported_not_used_as_donor") or ""))
    neg_a = re.sub(r"\[[^\]]*\]", "", str(r.get("reported_not_used_as_acceptor") or ""))
    neg = lab92.categorise(neg_d if neg_d.strip() not in ("", "nan") else None, lab92.DONOR, "donor_other") | \
        lab92.categorise(neg_a if neg_a.strip() not in ("", "nan") else None, lab92.ACCEPTOR, "acceptor_other")
    return {c: ("tested_positive" if c in pos else "tested_negative" if c in neg else "not_tested") for c in CATS}


def fetch():
    W.mkdir(parents=True, exist_ok=True)
    p = pd.read_csv(AUDIT, low_memory=False)
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    train = set(lab[lab.labelled == True].species.dropna())  # noqa: E712
    p["binomial"] = p.organism.astype(str).str.replace(r"\(.*?\)", "", regex=True).str.split().str[:2].str.join(" ")
    p = p[p.binomial.isin(train)]
    p = p[p.electron_donors_proof_dois.notna() | p.electron_acceptors_proof_dois.notna()]
    p = p.sample(frac=1, random_state=SEED)
    items = []
    for _, r in p.iterrows():
        ids = re.split(r";\s*", f"{r.electron_donors_proof_dois};{r.electron_acceptors_proof_dois}")
        pmids = [i for i in ids if i.isdigit()][:4]
        texts = []
        for pm in pmids:
            js = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pm}%20AND%20SRC:MED"
                     f"&resultType=core&format=json")
            if not js:
                continue
            res = json.loads(js).get("resultList", {}).get("result", [])
            if not res or res[0].get("isOpenAccess") != "Y" or not res[0].get("pmcid"):
                continue
            xml = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{res[0]['pmcid']}/fullTextXML")
            if xml:
                txt = re.sub(r"<[^>]+>", " ", xml)
                texts.append(re.sub(r"\s+", " ", txt))
            if sum(map(len, texts)) > MAX_CHARS:
                break
        if texts:
            items.append({"organism": r.binomial, "text": " ".join(texts)[:MAX_CHARS], "gold": gold_row(r)})
            print(len(items), r.binomial, len(items[-1]["text"]), flush=True)
        if len(items) >= N_ORG:
            break
    (W / "items.json").write_text(json.dumps(items))


def prompt(it):
    cats = "\n".join(f'- "{c}": {DESC[c]}' for c in DESC)
    return (f"Organism: {it['organism']}\n\nFrom the paper text below only, decide for each category whether the "
            f"organism was experimentally shown to have the capability (tested_positive), shown NOT to have it "
            f"(tested_negative, e.g. 'did not grow with', 'not used'), or the text does not report a test "
            f"(not_tested). Genome-based statements and statements about other organisms are not_tested.\n"
            f"Categories:\n{cats}\n\nAnswer with JSON only: {{\"<category>\": {{\"status\": \"...\", "
            f"\"quote\": \"verbatim sentence or empty\"}}, ...}} covering all {len(DESC)} categories.\n\n"
            f"PAPER TEXT:\n{it['text']}")


def call(model, text, key):
    body = json.dumps({"model": model, "messages": [{"role": "user", "content": text}], "temperature": 0,
                       "response_format": {"type": "json_object"}, "usage": {"include": True}}).encode()
    req = urllib.request.Request("https://openrouter.ai/api/v1/chat/completions", data=body,
                                 headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.load(r)


def run():
    key = Path.home().joinpath(".openrouter_key").read_text().strip()
    items = json.loads((W / "items.json").read_text())
    out_f = W / "responses.jsonl"
    done = set()
    spent = 0.0
    if out_f.exists():
        for l in out_f.read_text().splitlines():
            d = json.loads(l)
            done.add((d["model"], d["organism"]))
            spent += d.get("cost", 0) or 0
    for model in MODELS:
        use = model
        for it in items:
            if (model, it["organism"]) in done:
                continue
            if spent >= BUDGET:
                print(f"budget stop at ${spent:.3f}")
                return
            try:
                r = call(use, prompt(it), key)
            except Exception as e:  # noqa: BLE001
                if use in FALLBACK and use == model:
                    print(model, "unavailable synchronously -> using", FALLBACK[model], e)
                    use = FALLBACK[model]
                    r = call(use, prompt(it), key)
                else:
                    print("error", model, it["organism"], e)
                    continue
            cost = (r.get("usage") or {}).get("cost", 0) or 0
            spent += cost
            content = r["choices"][0]["message"]["content"]
            with open(out_f, "a") as f:
                f.write(json.dumps({"model": model, "served_by": use, "organism": it["organism"], "content": content,
                                    "cost": cost, "usage": r.get("usage")}) + "\n")
            print(f"{model} {it['organism']} ${cost:.4f} total ${spent:.3f}", flush=True)
    print(f"done, spent ${spent:.3f}")


def score():
    from sklearn.metrics import cohen_kappa_score
    items = {it["organism"]: it["gold"] for it in json.loads((W / "items.json").read_text())}
    R = {}
    costs = {}
    for l in (W / "responses.jsonl").read_text().splitlines():
        d = json.loads(l)
        try:
            ans = json.loads(re.search(r"\{.*\}", d["content"], re.S).group(0))
        except Exception:  # noqa: BLE001
            ans = {}
        R.setdefault(d["model"], {})[d["organism"]] = {c: (ans.get(c) or {}).get("status", "not_tested")
                                                       if isinstance(ans.get(c), dict) else "not_tested" for c in CATS}
        costs[d["model"]] = costs.get(d["model"], 0) + (d.get("cost") or 0)
    ref = "google/gemini-3.8-flash:batch"
    rows = []
    for m, preds in R.items():
        orgs = [o for o in preds if o in items]
        g = [items[o][c] for o in orgs for c in CATS]
        p = [preds[o][c] for o in orgs for c in CATS]
        gp = np.array([x == "tested_positive" for x in g]); pp = np.array([x == "tested_positive" for x in p])
        row = dict(model=m, organisms=len(orgs), kappa_vs_audit=cohen_kappa_score(g, p),
                   positive_recall=(gp & pp).sum() / max(gp.sum(), 1), positive_precision=(gp & pp).sum() / max(pp.sum(), 1),
                   negatives_found=sum(x == "tested_negative" for x in p), cost_usd=costs.get(m, 0))
        if m != ref and ref in R:
            common = [o for o in orgs if o in R[ref]]
            row["kappa_vs_gemini"] = cohen_kappa_score([R[ref][o][c] for o in common for c in CATS],
                                                       [preds[o][c] for o in common for c in CATS])
        rows.append(row)
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "pilot_scores.csv", index=False)
    print(df.round(3).to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "run":  # optional: restrict to listed models
        MODELS = sys.argv[2].split(",")
    {"fetch": fetch, "run": run, "score": score}[sys.argv[1]]()
