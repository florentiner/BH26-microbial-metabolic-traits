# Pre-registration — D: re-curating the training labels with LLM extraction (Haiku + Qwen 3.8 check)

Written 2026-10-15, before any training organism was submitted. Hash in `docs/preregistration_additional_experiments.sha256`.
Model choice from the pilot (`docs/RESULTS_llm_pilot.md`).

## Scope and blinding
* **Organisms:** the 355 species of the 358 labelled training genomes only (labels.tsv). Verified before writing:
  0 overlap with the newly described test species (validation_set.json, scored sets) and 0 with the external 562.
  No test organism is submitted to any model.
* **Categories:** the 12 donor/acceptor categories with a definition (the two catch-all bins are left unchanged).

## Evidence per organism
Up to 40 000 characters of literature text, in this order until the limit: (1) open-access full texts (Europe PMC)
of the PMIDs cited for that species in the audited tables (`hack_res_fedotmas/final/*_OPEN_ACCESS.csv`), (2) the
open-access full text of the species' description paper found by a Europe PMC query
`"<binomial>" AND ("sp. nov." OR "gen. nov.")`, (3) abstracts of the remaining cited PMIDs. Fixed now; no other
sources.

## Extraction and check
1. **Claude Haiku 5.5** (temperature 0) with the pilot prompt: status per category (tested_positive /
   tested_negative / not_tested) and a verbatim quote.
2. **Mechanical check:** a quote that does not occur verbatim (after whitespace normalisation) in the supplied text
   makes that cell not_tested.
3. **Qwen 3.8 Flash** (temperature 0) receives the same text and Haiku's non-"not_tested" cells with quotes and
   answers confirm / reject per cell. Only confirmed cells are kept.

## New label matrix ("curated_D")
For each training genome and category, starting from the current closed-world labels:
* current 1 → stays 1 (no positive is removed; the audit already corrected those);
* current 0 and confirmed tested_positive → 1 (new positive);
* current 0 and confirmed tested_negative → 0 (now a tested negative);
* current 0 otherwise → **unknown (NaN)**: excluded from training and from cross-validated scoring for that cell.
Catch-all categories unchanged.

## Models and evaluation
Retrain B2 (same learners, inner CV, frozen genus/family folds, 358 genomes) and B2XLE_large on `curated_D`.
Evaluate on data whose labels are unchanged and were not used for anything here:
* **primary:** the 43 held-out new species (75 blind growth tests);
* **secondary:** the external 562 (Madin), and the primary 100 (reported only).
Metrics: accuracy, MCC, sensitivity, specificity at each model's own released thresholds; AUROC; paired cluster
bootstrap (papers / genera) for curated − original.

## Predictions
* D-P1: on the held-out 43, B2-curated sensitivity exceeds B2's by ≥ 10 points.
* D-P2: on the held-out 43, B2-curated MCC ≥ B2 MCC (direction; CI reported).
* D-P3: on the external 562, B2-curated macro AUROC is not lower than B2's by more than 0.02.

## Budget
Hard stop at $6.00 total API cost (OpenRouter-reported); organisms are processed in random order (seed 20261015)
with both steps completed per organism, so a stop leaves a usable random subset.
