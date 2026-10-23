# Pre-registration — C: positive-unlabeled threshold correction for the released models

Written 2026-10-15, before any corrected threshold was computed or applied. Hash in
`docs/preregistration_additional_experiments.sha256`.

## Motivation and fairness
On the 100 primary new species B2's sensitivity was 52 %: the training labels are closed-world (an unreported
trait counts as negative), so thresholds learned on them are conservative. Because this observation came from the
primary 100, **claims about the correction rest only on data not used to motivate it**:
* **Primary test:** the 43 blind-evidence new species outside the primary 100 (75 growth tests, 41 papers;
  catch-all categories excluded as in the main analysis), never used for any comparison or decision.
* **Secondary test:** the external 562 lineage-novel genomes (Madin labels; also closed-world, so a fall in measured
  precision there is expected and is not by itself evidence against the correction).
* The primary 100 is reported for completeness only.
No model is refitted; the training data are the same 358 genomes; nothing from the test sets is used to set
anything below.

## Method (Elkan & Noto 2008, "e1" estimator; SCAR assumption)
For each model (B2 primary; B2XLE_large secondary) and each category k, using the **family-held-out out-of-fold**
scores of the 358 training genomes:
1. Fit a Platt calibrator g_k (logistic regression of the observed label on the raw OOF score).
2. Label frequency c_k = mean of g_k(score) over training genomes labelled positive, clipped to [0.3, 1].
3. Original decision: raw score ≥ t_k (the released MCC-optimal threshold). Its calibrated value is τ_k = g_k(t_k).
4. Corrected decision: g_k(score) / c_k ≥ τ_k, i.e. positive when the calibrated probability exceeds c_k · τ_k.
Nothing else changes (released / withheld categories, scores, ranking — AUROC is unaffected by construction).

## Metrics
Accuracy, MCC, sensitivity, specificity for original and corrected calls on each test set; paired cluster
bootstrap over papers (new species) or genera (external) for corrected − original; c_k reported per category.

## Predictions
* C-P1: on the 43 held-out species, B2's sensitivity rises by at least 10 points.
* C-P2: on the 43 held-out species, B2's MCC does not decrease (corrected − original ≥ 0; CI reported).
* C-P3: on the external set, specificity falls (direction only; expected under closed-world Madin labels).

## Deviations
(none yet)
