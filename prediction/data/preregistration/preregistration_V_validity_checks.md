# Pre-registration: V, two validity checks requested by an internal review

Registered 18 October 2026 before either analysis was run. SHA-256 in `docs/preregistration_V_validity_checks.sha256`.
Released models are not changed.

## What was looked at before registration
All results in `docs/RESULTS_Bclean_A_D.md`. For V1, only label counts of the subset (170 genomes; positives per
category: methanogen/acetogen 35, nitrate 100, oxygen 98 of 155 labelled, sulfate 5, sulfur reduction 26, sulfur
oxidation 1, iron and hydrogen categories 0). For V2, only that all 100 primary species belong to NCBI families absent
from training (so family membership cannot stratify). No score, call or distance has been computed for these checks.

## V1: RedoxKO against METABOLIC on external labels without FAPROTAX
**Question.** Does RedoxKO's advantage over METABOLIC on the external genomes hold when labels derived from FAPROTAX
(taxonomy-based) are removed?
* **Subset.** External genomes whose Madin `data_source` does not include FAPROTAX (the definition used in
  `src/22_external_eval.py`, subset `no_faprotax_source`), restricted to genomes with METABOLIC calls.
* **Categories.** Those with at least 5 positives and 5 negatives in the subset (expected: methanogen/acetogen,
  nitrate, oxygen, sulfate, sulfur reduction).
* **Endpoints.** Macro AUROC (primary) and macro MCC, RedoxKO minus METABOLIC; sensitivity and specificity; genus
  bootstrap, 2,000 resamples, seed 20261018. Success for "the external advantage holds": lower 95% bound of the macro
  AUROC difference above 0.
* **Expected in advance:** uncertain (about 50%). If the advantage disappears, the text must attribute the external
  result partly to label source.

## V2: accuracy against distance to the nearest training genome
**Question.** Does accuracy on newly described species fall with distance from the training genomes, and does it fall
more for RedoxKO (trained) than for METABOLIC (rule-based)?
* **Distance.** For each of the 100 primary species, 1 minus the Jaccard similarity of its full KO profile (the
  10,517 KOs of experiment F, same annotation) to the most similar of the 358 training genomes. A functional, not a
  phylogenetic, distance; stated as such.
* **Strata.** Tertiles of species by distance (about 33 species each).
* **Endpoints (descriptive, with paper-bootstrap 95% intervals, 2,000 resamples, seed 20261018).**
  Balanced accuracy of calls pooled over tests in each tertile, for RedoxKO and METABOLIC (primary); RedoxKO macro
  AUROC per tertile (categories with both classes); the difference far minus near tertile for each tool and the
  difference of these differences (RedoxKO minus METABOLIC). Spearman correlation between distance and per-test
  correctness for each tool.
* **Expected in advance:** RedoxKO declines with distance (about 65%); METABOLIC declines less (about 55%).

## Reporting
Both reported whatever the outcome, in `docs/RESULTS_Bclean_A_D.md` (section V), the appendix and, where they change a
claim, the main text.
