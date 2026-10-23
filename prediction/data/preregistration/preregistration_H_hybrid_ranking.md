# Pre-registration: H, combined ranking of METABOLIC calls and RedoxKO scores

Registered 17 October 2026, before any combined score was computed. SHA-256 in
`docs/preregistration_H_hybrid_ranking.sha256`. Released models are not changed. Runs after experiment D in the queue.

## Origin of the idea (stated so that readers can judge it)
The idea follows from results already seen: (i) the exploratory analysis on the 100 primary species, in which RedoxKO
scores ranked true above false positives among METABOLIC's positive calls (AUROC 0.66 [0.57, 0.75],
`src/108_ranking_vs_metabolic.py`), and (ii) experiment A, which confirmed this on the 43 withheld species (AUROC 0.74
[0.58, 0.91]). Call-level hybrids were tried earlier and did not improve on METABOLIC: METABOLIC HMMs as RedoxKO
features (B2M) and METABOLIC calls for categories RedoxKO does not release (`results/additional/hybrid/hybrid_results.csv`).

## What was looked at before registration
* Call-level metrics (accuracy, MCC, sensitivity, specificity) and macro AUROC of RedoxKO, and call-level metrics of
  METABOLIC, on both sets (`x1_summary.csv`, `hybrid_results.csv`), plus experiment A.
* Missing values only for the inputs of this test (none in `call_METABOLIC` or `score_B2`; scores lie strictly
  between 0 and 1).
* No combined score has been computed. Because of experiment A, the comparison with METABOLIC (HH2) is partly
  predictable, so the comparison with RedoxKO (HH1) is primary.

## The rule (one variant only, no tuning)
For each growth test, combined score = `call_METABOLIC` + `score_B2`, where `call_METABOLIC` is 0 or 1 (categories
METABOLIC does not cover count as 0, as in the tool comparison) and `score_B2` is the released, calibrated RedoxKO score
in (0, 1). Tests are thus ordered by METABOLIC's call first and by RedoxKO's score within each call. No weights, no
other tools, no other RedoxKO variant.

## Data
* Primary: the 255 growth tests of the 100 primary new species (`results/additional/x1_calls.csv`).
* Confirmation: the 75 tests of the 43 withheld species (`results/additional/hybrid/holdout_calls.csv`).

## Endpoints
Macro AUROC = mean over categories with both classes in the set, as in `src/87_tool_comparison.py`. METABOLIC's AUROC
uses its 0/1 call as the score. Paired bootstrap over description papers, 2,000 resamples, seed 20261017; categories
lacking both classes in a resample are skipped in that resample.
* **HH1 (primary):** macro AUROC, combined minus RedoxKO, primary set. Success: lower 95% bound above 0.
* **HH2:** macro AUROC, combined minus METABOLIC, primary set. Success: lower 95% bound above 0.
* **HH3 (confirmation):** HH1 and HH2 on the withheld set. Same sign as on the primary set counts as confirmation;
  significance is reported but not required (75 tests).
* **Secondary, descriptive:** pooled AUROC over all tests; per-category AUROC.
* **Expected result, stated in advance:** HH2 likely positive (consequence of experiment A); HH1 about 55-65% likely
  to succeed on the primary set and 30-40% likely to be significant on the withheld set.

## Reporting
All endpoints are reported whatever the outcome, in `docs/RESULTS_Bclean_A_D.md`, the supplement and, if HH1 succeeds,
the Results, described as a combination rule motivated by earlier results, not as a new model.
