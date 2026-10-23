"""Score model calls for the 10 new organisms against the published tests (selection.md, scoring rule).

usage: python validation_new_organisms/score.py predictions_B2.csv [predictions_B2XLE_large.csv ...]
Writes scored_<model>.csv (one row per organism x category with published status, call, score) and prints a summary.
"""
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
cands = {c["refseq_accession"]: c for c in json.load(open(HERE / "candidates.json"))}
sel = [a.strip() for a in open(HERE / "selected_accessions.txt") if a.strip()]
EXPECT = {"tested_positive": 1, "tested_negative": 0}
# traits for which the released predictor returns yes/no calls (out-of-clade skill; results/predictions_unlabelled_B2.csv)
RELIABLE = {"acceptor_CO2_methanogen_acetogen", "acceptor_nitrate", "acceptor_oxygen", "acceptor_sulfate", "acceptor_sulfur",
            "donor_ammonia_oxidation", "donor_hydrogen_oxidation", "donor_sulfur_oxidation"}

for f in sys.argv[1:]:
    p = pd.read_csv(HERE / f)
    arm = p.arm.iloc[0]
    rows = []
    for acc in sel:
        c = cands[acc]
        for t, ev in c["categories"].items():
            q = p[(p.genome_id == acc) & (p.trait == t)]
            if q.empty:
                continue
            q = q.iloc[0]
            status = ev.get("status", "not_tested")
            exp = EXPECT.get(status)
            rows.append(dict(accession=acc, organism=c["species"], category=t, published=status,
                             substrates=ev.get("substrates", ""), expected=exp, call=int(q.predicted),
                             score=round(float(q.score), 3), calibrated=round(float(q.calibrated_score), 3),
                             threshold=round(float(q.threshold), 3),
                             agree=(None if exp is None else int(exp == int(q.predicted))),
                             source=ev.get("source", "")))
    d = pd.DataFrame(rows)
    d.to_csv(HERE / f"scored_{arm}.csv", index=False)
    s = d[d.expected.notna()]
    pos, neg = s[s.expected == 1], s[s.expected == 0]
    print(f"\n== {arm}: {len(s)} scored calls on {s.accession.nunique()} organisms")
    print(f"   correct {int(s.agree.sum())} / {len(s)} ({s.agree.mean():.0%});  tested positives found "
          f"{int(pos.agree.sum())} / {len(pos)} ({pos.agree.mean():.0%});  tested negatives respected "
          f"{int(neg.agree.sum())} / {len(neg)} ({neg.agree.mean():.0%})")
    r = s[s.category.isin(RELIABLE)]
    rp, rn = r[r.expected == 1], r[r.expected == 0]
    print(f"   deployed mode (only the 8 categories with released yes/no calls): correct {int(r.agree.sum())} / {len(r)} "
          f"({r.agree.mean():.0%}); positives {int(rp.agree.sum())}/{len(rp)}, negatives {int(rn.agree.sum())}/{len(rn)}")
    by = s.groupby("category").agree.agg(["sum", "count"])
    print("   per category (correct / scored):", {k: f"{int(v['sum'])}/{int(v['count'])}" for k, v in by.iterrows()})
    mism = s[s.agree == 0][["organism", "category", "published", "call", "score", "substrates"]]
    print(f"   mismatches ({len(mism)}):")
    print(mism.to_string(index=False))
    nt = d[d.expected.isna() & (d.call == 1)]
    print(f"   predicted positive where the literature did not test ({len(nt)}):",
          sorted(set(zip(nt.organism.str.split().str[0], nt.category))))
