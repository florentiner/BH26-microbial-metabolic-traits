# Pre-registration CB: class balancing for the donor and acceptor model

Registered 20 October 2026 before any of these models was fitted. SHA-256 in
`docs/preregistration_CB_class_balance.sha256`. Descriptive sensitivity analysis requested in review preparation; the
final models are not changed by the outcome.

**Question.** The 14 donor and acceptor categories are imbalanced (3-58% positive genomes among the 358 labeled
genomes). RedoxKO handles this by weighting classes by inverse frequency (logistic regression `class_weight=
"balanced"`, XGBoost `scale_pos_weight = negatives / positives`). Would a different treatment of the imbalance give
better predictions?

**Known before registration.** Only the existing RedoxKO results (arm `B2`): macro AUPRC 0.576 genus level, 0.560
family level. None of the variants below has been run.

**Arms** (everything else identical to `B2`: 852-KO panel, frozen genus- and family-level outer folds, inner grouped
CV over the same six configurations, MCC threshold from inner out-of-fold predictions, `redox_common.nested_cv`):
* `CB_none`: no class weighting (`class_weight=None`, `scale_pos_weight=1`).
* `CB_under`: random undersampling of negatives to the number of positives in every training set, no weighting;
  probabilities averaged over 5 undersampled draws.
* `CB_over`: random oversampling of positives (with replacement) to the number of negatives in every training set,
  no weighting.
Resampling is applied only to training data (inner and outer), never to validation or test genomes.

**Endpoints.** Primary: macro AUPRC over 14 categories in genus-level CV, paired difference arm − `B2` with a
genus-cluster bootstrap (`redox_common.paired_macro_delta`, 1,000 resamples, seed 20260917). Secondary: the same at
family level; macro MCC of the yes/no calls at the inner-CV thresholds; macro AUPRC over the categories with fewer
than 10% positives.

**Expected.** No arm clearly better than `B2` (all paired genus-level intervals include zero or lie below it),
probability about 75%. `CB_none` lower on the rare categories (about 60%), because without weighting the rare
positives contribute little to the loss. `CB_under` lower overall (about 65%), because it discards most negatives
of a 358-genome set. A variant "clearly better" means a paired genus-level interval entirely above zero.

**Reporting.** All three arms go into Supplementary Table S3 regardless of their outcome, and the Methods sentence on
class weighting states the result.
