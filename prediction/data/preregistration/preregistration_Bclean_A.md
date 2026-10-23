# Pre-registration: B-clean (more non-taxonomic labels) and A (ranking within METABOLIC calls)

Registered 17 October 2026, before any genome was annotated, any model was trained or any new score was computed.
The SHA-256 of this file is recorded in `docs/preregistration_Bclean_A.sha256`. The released models
(`results/models/B2`, `B2XLE_large`, pathway models) are not changed by this experiment, whatever the outcome.

## What was looked at before registration
* Counts only: genomes available, clean positives per category, number of METABOLIC-positive calls (44) and
  their labels (38 true, 6 false) on the 43 withheld species, used for a power simulation (26% power for A).
* No KO-model score of the 43 withheld species was inspected for this purpose. No B-clean model exists yet.

## B-clean: question and design
**Question.** Does adding genomes whose trait records are not derived from taxonomy (FAPROTAX) improve the KO model
on newly described species?

**Added genomes.** Genomes of the pathway dataset (`data/pathways/labels.tsv`, Madin et al. 2020 compilation) with
`faprotax_source == False`, after removing (i) every genus in the external evaluation set, (ii) every genus of any
newly described candidate species (`validation_new_organisms/candidates.json`, scored and pathway sets; 173 genera),
and (iii) every genome already in the training set. This leaves 546 genomes (list:
`results/additional/bclean/added_genomes.tsv`).

**Labels for added genomes.** Only categories with at least 20 clean positives among the added genomes receive
labels: `acceptor_nitrate` (Madin nitrate_reduction or denitrification or dnra; 457 positives) and
`acceptor_CO2_methanogen_acetogen` (methanogenesis or acetogenesis; 36 positives). A trait not listed in the record
is negative (closed-world, as in the training labels). All other categories are missing for added genomes, so their
models are not retrained and equal the released B2.

**Features.** The 852-KO panel computed for the added genomes with the training pipeline (NCBI proteins, pyhmmer,
KOfam profiles and thresholds of the 25 August 2026 release), through `src/predict.py` functions
`prepare_proteins` and `ko_features`. Genomes that fail to download or have no NCBI annotation are dropped and
listed.

**Model.** For the two categories, the released B2 procedure unchanged (`src/09_train_final.py select_and_fit`):
the same six configurations, inner grouped cross-validation by genus, inverse-frequency class weights, MCC-optimal
threshold on inner out-of-fold predictions, Platt calibration, fitted on the 358 training genomes plus the added
genomes. Added genomes form their own genus groups (training genera keep their groups).

**Endpoints.**
* **Primary (H1).** AUROC of `acceptor_nitrate` scores on the 57 nitrate growth tests of the 100 primary new species,
  B-clean minus released B2, paired bootstrap over description papers (2,000 resamples, seed 20261017).
  Success: lower 95% bound above 0.
* **H2.** Same comparison for MCC at the released thresholds (B-clean uses its own threshold).
* **H3 (confirmation).** H1 and H2 on the nitrate tests of the 43 withheld species; success if the difference has
  the same sign (no significance required, given the size).
* **H4 (non-inferiority, internal).** Genus-blocked cross-validation on the frozen folds of the 358 genomes, with
  added genomes used only in training folds: nitrate AUPRC of B-clean minus B2, margin −0.02 (lower bound above
  −0.02 counts as non-inferior).
* **Descriptive only.** Methanogen/acetogen results (5 new-species tests; too few for inference); macro AUROC over
  all 255 tests (other categories unchanged by design). The external 562-genome set is not used as evidence,
  because it comes from the same compilation as the added labels.

## A: ranking within METABOLIC's positive calls (confirmation of an exploratory result)
* **Exploratory result being confirmed.** On the 100 primary species, among tests METABOLIC called positive, the
  released B2 scores ranked true above false positives with AUROC 0.66 [0.57, 0.75] (`src/108_ranking_vs_metabolic.py`).
* **Test (HA).** The same statistic for the released B2 on the 43 withheld species (44 METABOLIC-positive tests,
  6 false positives), paper-cluster bootstrap (3,000 resamples, seed 20261017). Success: lower 95% bound above 0.5.
  Expected power about 26%, so an inconclusive result is the most likely outcome and will be reported as such.
* **Descriptive.** The same statistic for B-clean.

## Reporting
All endpoints are reported whatever their outcome, in `docs/RESULTS_Bclean_A.md` and the paper's supplement.
Any deviation from this plan is logged with its reason in `docs/preregistration_additional_experiments_deviations.md`.
