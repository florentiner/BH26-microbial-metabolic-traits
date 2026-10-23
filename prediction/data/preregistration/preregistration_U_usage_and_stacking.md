# Pre-registration: U, usage analyses (top-k, decision curve, selective prediction) and S, learned stacking

Registered 18 October 2026, before any of the statistics below was computed. SHA-256 in
`docs/preregistration_U_usage_and_stacking.sha256`. Released models are not changed.

## What was looked at before registration
Everything reported in `docs/RESULTS_Bclean_A_D.md` (A, B-clean, D, H, F, M1-M4), including RedoxKO's and METABOLIC's
sensitivity, specificity, precision, AUROC and MCC on the external set and new species, and the conformal summary
(singleton accuracy 79% on new species). Not computed so far: precision at fixed k, net benefit, METABOLIC's accuracy
on the subset of tests with a confident RedoxKO label, and any stacked model.

## Data (as in experiment M)
* External: 562 genomes, 9 categories, raw RedoxKO scores (`results/external_validation_predictions.csv`, arm B2),
  RedoxKO thresholds, METABOLIC calls from `src/metabolic_calls.py`; bootstrap over genera.
* New species: 255 tests of the 100 primary species (`x1_calls.csv`, real METABOLIC calls) and 75 tests of the 43
  withheld species (`hybrid/holdout_calls.csv`); bootstrap over description papers.
* Bootstrap: 2,000 resamples, seed 20261017. Macro = mean over categories (with both classes where required).

## U1: top-k screening precision
For each category, precision among the k genomes ranked highest by RedoxKO's score, against METABOLIC's expected
precision for k genomes drawn at random from its positive calls (if k exceeds its positive calls, the remainder is
drawn at random from its negative calls). **Primary (external):** k = 20, macro over the 9 categories, RedoxKO minus
METABOLIC; success: lower 95% bound above 0. Secondary (new species, primary set): k = the top quarter of tests in
each category with at least 8 tests.

## U2: decision-curve analysis
Net benefit NB(p_t) = TP/N − FP/N · p_t/(1 − p_t), per category, then macro. RedoxKO calls positive when its
**calibrated** probability (the released Platt calibrator) is at least p_t; METABOLIC's calls do not depend on p_t.
Thresholds p_t = 0.05, 0.10, …, 0.90. **Reported (external, descriptive with intervals):** the p_t values at which the
95% interval of NB(RedoxKO) − NB(METABOLIC) lies above 0, and those at which it lies below 0. Same on new species
(secondary).

## U3: selective prediction
RedoxKO's Mondrian conformal sets as in `src/89_conformal.py` (quantiles `results/additional/x3_quantiles.json`, 90%
target). On the tests where RedoxKO gives a single confident label: accuracy of that label against the accuracy of
METABOLIC's calls on the same tests. **Primary (new species, primary set):** difference RedoxKO minus METABOLIC;
success: lower 95% bound above 0. Secondary: external set; fraction of tests covered.

## S: learned stacking of RedoxKO and METABOLIC
For each of the 12 categories that METABOLIC maps, a logistic regression (L2, C = 1, class weights balanced) on two
features, logit of RedoxKO's genus-held-out out-of-fold score and METABOLIC's call (0/1), fitted on the 358 training
genomes (METABOLIC calls from `data/processed/genome_features_metabolic_hmm.parquet`). Decision threshold: maximizing
MCC on out-of-fold predictions of the stacker from 5-fold cross-validation grouped by genus on the training genomes.
Applied unchanged to new species and external genomes (feature: logit of the released RedoxKO score). Categories
without a METABOLIC mapping keep RedoxKO.
* **HS1 (primary):** pooled MCC on the 255 primary tests (as in the tool comparison), stack minus METABOLIC.
  Success: lower 95% bound above 0.
* **HS2 (confirmation):** the same on the 75 withheld tests; same sign counts as confirmation.
* **HS3:** external macro MCC, stack minus METABOLIC, and stack minus RedoxKO.
* **HS4:** macro AUROC, stack minus RedoxKO, on new species and external (does stacking also help ranking?).

## Expected results, stated in advance
U1: ~80% that RedoxKO wins on external; U2: RedoxKO better at higher p_t, METABOLIC at low p_t; U3: ~50%;
S: HS1 25-35%, HS3 60-70%.

## Reporting
All endpoints are reported in `docs/RESULTS_Bclean_A_D.md` (sections U and S) and the supplement.
