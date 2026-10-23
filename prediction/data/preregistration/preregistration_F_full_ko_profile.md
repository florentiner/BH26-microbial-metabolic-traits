# Pre-registration: F, full KO profile for the donor and acceptor models

Registered 17 October 2026, before any full-profile feature was computed for the donor and acceptor genomes or the
newly described species, and before any model was trained. SHA-256 in `docs/preregistration_F_full_ko_profile.sha256`.
Released models are not changed.

## Question
RedoxKO uses an 852-KO energy-metabolism panel. The pathway models gained from the full KO profile. Does the full
profile improve RedoxKO's ranking of growth tests on newly described species?

## What was looked at before registration
* RedoxKO results on all sets (as published in the text), experiments A, B-clean, D and H.
* The TogoMCP KO matrix of the 586 donor/acceptor genomes, only to compare it with the local annotation of the
  852-KO panel: 97.1% agreement, but 6,086 local-only hits against 3 TogoMCP-only hits. Because of this asymmetry the
  TogoMCP matrix is **not** used; all features below come from the same local annotation.
* No full-profile feature, model or score for any donor/acceptor genome or new species has been computed.

## Features
Presence/absence of the 10,517 KOs of the released pathway model EnergyKO-M (a list fixed before this experiment),
computed with `src/predict.py` `ko_features_for` (pyhmmer, KOfam profiles and thresholds of the 25 August 2026
release, KofamScan scoring rules) for: the 358 training genomes (proteins from `data/processed/proteins`), and the
143 newly described species of the primary and withheld sets (proteins from `data/predict_cache/proteins`). The 99
species already in `validation_new_organisms/pathway100_ko_matrix.parquet`, computed with the same function and
release, are reused. Columns constant across the 358 training genomes are dropped.

## Model (RedoxKO-full)
The released RedoxKO procedure unchanged (`src/09_train_final.py select_and_fit`, the six `ko_configs()`
configurations, inner grouped cross-validation by genus, inverse-frequency class weights, MCC-optimal threshold on
inner out-of-fold predictions, Platt calibration), for all 14 categories, on the full-profile features instead of the
panel. No other change.

## Endpoints
Paired bootstrap over description papers, 2,000 resamples, seed 20261017; macro AUROC = mean over categories with
both classes, as in `src/87_tool_comparison.py`.
* **HF1 (primary):** macro AUROC on the 255 growth tests of the 100 primary species, RedoxKO-full minus RedoxKO.
  Success: lower 95% bound above 0.
* **HF2:** MCC at each model's own thresholds, same comparison (sensitivity and specificity reported).
* **HF3 (confirmation):** HF1 and HF2 on the 75 tests of the 43 withheld species; same sign counts as confirmation.
* **HF4 (internal):** genus-blocked cross-validation on the frozen folds, macro AUPRC, RedoxKO-full minus the
  released RedoxKO out-of-fold predictions, genus-cluster bootstrap (2,000 resamples). Reported with its interval;
  non-inferiority margin −0.02.
* **Expected result, stated in advance:** about 15–20% chance that HF1 succeeds. More features with 358 genomes may
  add lineage signal rather than gene evidence, so a null or negative result is the most likely outcome.

## Reporting
All endpoints are reported whatever the outcome, in `docs/RESULTS_Bclean_A_D.md` (section F) and the supplement. Any
deviation is logged in `docs/preregistration_additional_experiments_deviations.md`.
