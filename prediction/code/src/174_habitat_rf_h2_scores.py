"""X4 habitat tests for RedoxKO after the release-rule amendment of 22 October 2026 (H2 oxidation returns a score only).

Same 300 MAGs, predictions and expectations as src/172_habitat_rf.py. Categories with released yes/no calls are tested
with one-sided Fisher's exact tests on the calls; score-only categories (Fe(II) oxidation, and now H2 oxidation) with
one-sided Mann-Whitney tests on the raw scores, as the registration prescribes for categories without released calls.
All p values are corrected together (BH).
Output: results/additional/rf_habitat/x4_tests.csv (overwritten; the earlier file is kept as x4_tests_calls_h2.csv)
"""
import json
import shutil
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "data/additional/habitat"
OUT = ROOT / "results/additional/rf_habitat"
rel = set(json.loads((ROOT / "results/models/RF/metadata.json").read_text())["reliable_traits"])

sel = pd.read_csv(WORK / "selection.tsv", sep="\t").set_index("mag")
P = pd.concat([pd.read_csv(p) for p in sorted((WORK / "pred_rf").glob("batch_*.csv"))], ignore_index=True)
g = "genome" if "genome" in P else P.columns[0]
sc = P.pivot_table(index=g, columns="trait", values="score")
cl = P.pivot_table(index=g, columns="trait", values="predicted")
sc["habitat"] = cl["habitat"] = sel.habitat.reindex(sc.index)

EXPECT = [("X4-E1", "donor_hydrogen_oxidation", "hydrothermal", ["marine_surface"]),
          ("X4-E2", "donor_sulfur_oxidation", "hydrothermal", ["human_gut"]),
          ("X4-E3", "donor_iron_oxidation", "acid_mine", ["hydrothermal", "human_gut", "digester", "marine_surface"]),
          ("X4-E4", "acceptor_CO2_methanogen_acetogen", "digester", ["marine_surface"]),
          ("X4-E5", "acceptor_oxygen", "marine_surface", ["human_gut", "digester"])]
rows = []
for name, trait, a, others in EXPECT:
    for b in others:
        if trait in rel:
            x, y = cl[cl.habitat == a][trait], cl[cl.habitat == b][trait]
            t = [[int(x.sum()), int(len(x) - x.sum())], [int(y.sum()), int(len(y) - y.sum())]]
            p, fa, fb, kind = fisher_exact(t, alternative="greater")[1], x.mean(), y.mean(), "Fisher one-sided on released calls"
        else:
            x, y = sc[sc.habitat == a][trait], sc[sc.habitat == b][trait]
            p, fa, fb = mannwhitneyu(x, y, alternative="greater").pvalue, x.median(), y.median()
            kind = "Mann-Whitney one-sided on raw scores (medians shown)"
        rows.append(dict(test=f"{name} vs {b}", trait=trait, higher_in=a, lower_in=b, frac_a=fa, frac_b=fb, p=p, kind=kind))
T = pd.DataFrame(rows)
T["q_bh"] = multipletests(T.p, method="fdr_bh")[1]
T["supported"] = T.q_bh < 0.05
f = OUT / "x4_tests.csv"
if f.exists() and not (OUT / "x4_tests_calls_h2.csv").exists():
    shutil.copy(f, OUT / "x4_tests_calls_h2.csv")
T.to_csv(f, index=False)
print(T.round(4).to_string(index=False))
