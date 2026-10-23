# Addendum C — Microbe Decoder (machine-learning tool) as an X1 comparator

Written 2026-10-14 ~23:30 MSK, before any genome of the evaluation set was submitted. Microbe Decoder
(Hackmann et al., NAR 2026; public API https://api.microbe-decoder.org, no login) predicts traits from KO-annotated
genes (`POST /compute/ml`, `gene_functions`). Input per genome: every KO present in its full KOfam profile
(10 517 KOs; release 2026-08-25, our thresholds), one row per KO. Organism names are **not** sent, so the tool cannot
use its organism database. A trait is predicted positive at probability ≥ 0.5.

| Our category | Microbe Decoder ML model | Rule |
|---|---|---|
| acceptor_oxygen | Anaerobe (oxygen tolerance) | positive when the anaerobe model is < 0.5 |
| acceptor_nitrate | Nitrate (metabolites utilized) | positive at ≥ 0.5 |
| acceptor_CO2_methanogen_acetogen | Methanogenesis (type of metabolism) | positive at ≥ 0.5 (acetogens not covered) |
| all other categories | not covered | scored negative in the all-call table; primary comparison on the three covered categories |

Primary metric: accuracy and MCC on the calls of the three covered categories, B2 on the same calls, paired cluster
bootstrap over papers; AUROC from the probabilities. The taxonomy tool (FAPROTAX-based) is not run separately
(FAPROTAX is already a comparator); the metabolic-network tool is not run (it predicts reaction pathways, not these
categories).
