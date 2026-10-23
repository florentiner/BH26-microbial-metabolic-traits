# Pre-registration: D, positive-unlabeled bagging for the KO model

Registered 17 October 2026, before any PU model was trained or scored. SHA-256 in
`docs/preregistration_D_pu_bagging.sha256`. Released models are not changed by this experiment.

## Question
Training labels are closed-world: a trait not reported for a genome counts as negative. An earlier correction
(experiment C) shifted decision thresholds by an estimated reporting rate. Shifting thresholds cannot change the ranking
of genomes. Does positive-unlabeled (PU) training, which changes how the model is fitted, improve the ranking of growth
tests on newly described species?

## Method
Bagging PU learning (Mordelet and Vert): for each of the 14 categories, positives P are the genomes labeled positive and
U the genomes labeled negative (treated as unlabeled). For each of K = 20 bags, a random subset of U of size
min(|P|, |U|), drawn without replacement (seed 20261017 + bag), serves as negatives, and the base classifier is fitted on
P plus that subset. Scores are averaged over bags.
* **Base classifier:** for each category, the configuration selected for the released B2 (logistic regression with its
  C, or XGBoost with its depth), from `results/models/B2/model.joblib`; the same 852-KO panel features (columns varying
  in training), the same class weighting.
* **Threshold:** maximizes MCC on inner out-of-fold predictions of the bagged model, using the B2 inner grouped folds
  (`inner_folds`, grouped by genus), never on test data. Platt calibration as for B2.
* Training data: the 358 labeled genomes only (no B-clean genomes), so the comparison isolates the training method.

## Endpoints
* **Primary (HD1):** macro AUROC (mean over categories with both classes, as in `src/87_tool_comparison.py`) on the 255
  growth tests of the 100 primary new species, PU-bagged minus released B2. Paired bootstrap over description papers,
  2,000 resamples, seed 20261017. Success: lower 95% bound above 0.
* **HD2:** MCC, sensitivity and specificity at each model's own thresholds, same comparison.
* **HD3 (confirmation):** HD1 and HD2 on the 75 tests of the 43 withheld species; same sign counts as confirmation.
* **HD4 (non-inferiority, internal):** genus-blocked cross-validation on the frozen folds; macro AUPRC of PU-bagged minus
  the released B2 out-of-fold predictions, genus-cluster bootstrap; margin −0.02.
* **Expected result, stated in advance:** sensitivity up, macro AUROC within ±0.02 of B2, MCC within noise. A
  significant gain is considered unlikely (about 10–15%).

## Reporting
All endpoints are reported whatever the outcome, in `docs/RESULTS_Bclean_A_D.md` and the supplement.
