# Pre-registration — additional compute-only experiments for the BH26 paper

Written 2026-10-14, **before** any comparator tool was run, any MAG was downloaded, or any conformal or habitat
analysis was computed. The SHA-256 of this file is stored in `docs/preregistration_additional_experiments.sha256`
and every script below records it in its output. Deviations will be logged at the end of this file with a date.

Models evaluated are the saved, already-released models only (`results/models/B2`, `B2XL`, `B2XLE_large`,
`pathways/PB3`); **no model is refitted or re-tuned** in any experiment below. Unless stated otherwise, intervals are
95 % cluster bootstraps (1 000 resamples) over the unit named, seed 20261014.

---

## X1. Comparison with existing tools on the newly described species

**Question.** On the primary set of 100 newly described species (255 blind growth tests, 9 categories;
`validation_new_organisms/balanced_prospective_primary_100.txt`), how do our models compare with tools that a
reader could use instead?

**Comparators.**
1. **FAPROTAX** (latest release, Louca et al. 2016): taxonomy-based. Lineage = current NCBI taxonomy of each
   species (domain … genus, species). Matching uses FAPROTAX's own `collapse_table.py` if distributed with the
   release, otherwise wildcard matching of each FAPROTAX taxon pattern against the lineage string. Group → category
   mapping (fixed now):
   * acceptor_oxygen: aerobic_chemoheterotrophy, aerobic_ammonia_oxidation, aerobic_nitrite_oxidation, nitrification
   * acceptor_nitrate: nitrate_respiration, nitrate_reduction, nitrate_denitrification, denitrification, nitrogen_respiration
   * acceptor_sulfate: sulfate_respiration
   * acceptor_sulfur: sulfur_respiration, thiosulfate_respiration, sulfite_respiration
   * acceptor_iron: iron_respiration
   * acceptor_CO2_methanogen_acetogen: every group whose name starts with `methanogenesis` or `hydrogenotrophic_methanogenesis`, plus `acetogenesis` if present
   * donor_sulfur_oxidation: dark_oxidation_of_sulfur_compounds, dark_sulfur_oxidation, dark_sulfide_oxidation, dark_sulfite_oxidation, dark_thiosulfate_oxidation
   * donor_hydrogen_oxidation: dark_hydrogen_oxidation
   * donor_iron_oxidation: dark_iron_oxidation
   * donor_ammonia_oxidation: aerobic_ammonia_oxidation; donor_nitrite_oxidation: aerobic_nitrite_oxidation
   * donor_CO_oxidation: none (FAPROTAX has no such group; scored as "not covered")
   A category is predicted positive if any mapped group is assigned, otherwise negative. Group names absent from
   the downloaded release are dropped and listed.
2. **Taxonomic nearest-group baseline** on the new species: mean training label (358 genomes) of the lowest shared
   rank (genus → family → order → class → phylum → domain → prevalence), same rule as `null_taxNN`. Lineage as above.
3. **Marker rule (arm A)**: the saved rule model applied to the cached 852-KO profiles.
4. **BacDive-AI** (Koblitz et al. 2025; github.com/JKoblitz/bacdive-AI), `aerobic` model, for **acceptor_oxygen
   only**. Pfam annotation with pyhmmer `hmmsearch` against Pfam-A using gathering thresholds (InterProScan Pfam is
   not installable locally; this is logged as a deviation from the tool's documented input). Run in a separate
   environment with scikit-learn 1.4.0 as the tool requires.

**Not run:** Microbe Decoder (web server only), DRAM / METABOLIC (database size; arm A is the rule-based stand-in).

**Metrics.** Per tool: coverage (fraction of organisms with ≥ 1 assignment; FAPROTAX only), accuracy, MCC,
sensitivity, specificity on all 255 calls, and on the calls of covered organisms only; AUROC where the tool gives
a score (taxonomic baseline, B2, B2XLE_large, BacDive-AI). Paired differences B2 − tool on identical calls,
cluster bootstrap over description papers.

**Predictions.**
* X1-P1: FAPROTAX covers fewer than 50 % of the 100 species (they come from genera absent from our data, and
  FAPROTAX is largely genus-indexed).
* X1-P2: On all 255 calls, B2's MCC exceeds FAPROTAX's MCC (CI excluding 0).
* X1-P3: B2's AUROC exceeds the taxonomic baseline's AUROC (CI excluding 0).
* X1-P4: For acceptor_oxygen, BacDive-AI's accuracy is within 10 points of B2 (no directional prediction; it is a
  dedicated oxygen model trained on far more strains).

---

## X2. Real metagenome-assembled genomes (MAGs)

**Question.** Do predictions made from real MAGs agree with predictions from the isolate genome of the same
species, and how does agreement degrade with completeness?

**Pairs.** Species clusters of the latest GTDB release that contain (a) one of our evaluation isolate genomes
(external 562 set or the 143 blind new species) and (b) at least one genome flagged by GTDB/NCBI as derived from a
metagenome. Up to 3 MAGs per species, chosen at random (seed 20261014) after stratifying by GTDB CheckM2
completeness bins [50, 70), [70, 90), [90, 100]. MAG nucleotide FASTA from NCBI; genes called with Pyrodigal as
in `src/predict.py --fasta`. Model: **B2** (CPU path; the fusion needs GPU embeddings and is excluded).

**Metrics.** (1) Call agreement MAG vs isolate per category and macro, by completeness bin; (2) mean absolute
score difference; (3) where the isolate has labels (external Madin labels or new-species growth tests): AUROC and
MCC from MAG predictions vs from the paired isolate predictions on the same label set, paired bootstrap over
species.

**Predictions.**
* X2-P1: macro call agreement ≥ 0.90 for MAGs ≥ 90 % complete.
* X2-P2: agreement decreases monotonically across the three completeness bins.
* X2-P3: errors are asymmetric: MAG-only negatives (isolate positive, MAG negative) outnumber MAG-only positives.

---

## X3. Per-prediction uncertainty (conformal prediction)

**Method.** Mondrian (class-conditional, per category) split-conformal prediction for B2 with non-conformity
1 − p(true class). Calibration = family-held-out out-of-fold scores of the 358 training genomes (`results/oof`).
Target coverage 90 %. Each prediction becomes {1}, {0}, {0,1} (uncertain) or {} .

**Evaluation sets.** External 562 (Madin labels) and the 100 new species (growth tests).

**Metrics.** Empirical coverage per class and macro; fraction of singleton sets; accuracy of singleton
predictions; same split by category.

**Predictions.**
* X3-P1: empirical coverage on the new species ≥ 0.85 (some loss expected from distribution shift).
* X3-P2: accuracy of singleton predictions on the new species exceeds B2's overall accuracy (64.7 %) by ≥ 5 points.

---

## X4. Habitat case study on uncultured lineages

**Question.** Are predicted donor/acceptor profiles of MAGs from different habitats consistent with the
chemistry of those habitats?

**Genomes.** GTDB species representatives that are MAGs (derived from metagenome), CheckM2 completeness ≥ 90 %
and contamination ≤ 5 %, whose NCBI isolation source / biosample text matches exactly one habitat class by the
keyword lists below (case-insensitive), at most 60 per class chosen at random (seed 20261014), at most one per
GTDB genus. Classes and keywords (fixed now):
* hydrothermal: "hydrothermal", "vent", "black smoker", "chimney"
* acid_mine: "acid mine", "mine drainage", "AMD"
* human_gut: "human gut", "human feces", "human faeces", "human stool"
* digester: "anaerobic digester", "biogas", "digestor"
* marine_surface: "seawater", "surface water", "ocean water" (excluding any hydrothermal keyword)
Model: B2 (raw scores and released calls).

**Pre-registered expectations** (Fisher exact test on predicted positives, one-sided, BH across the list):
* X4-E1: H2 oxidation more frequent in hydrothermal than in marine_surface.
* X4-E2: sulfur oxidation more frequent in hydrothermal than in human_gut.
* X4-E3: Fe(II) oxidation more frequent in acid_mine than in every other class.
* X4-E4: CO2 reduction (methanogen/acetogen) more frequent in digester than in marine_surface.
* X4-E5: O2 respiration less frequent in human_gut and digester than in marine_surface.
Iron-oxidation calls are not released by B2; for X4-E3 the raw score is compared (Mann–Whitney, one-sided).
A result contrary to an expectation is reported as such; habitats are coarse proxies, so this is a
plausibility check, not a validation.

---

## X5. Human check of the LLM-extracted growth tests (sheet only)

A stratified random sample of 50 calls from the 255 primary calls (seed 20261014; at least 3 per category where
available) is exported with paper, quote and extracted label, with blank columns for two independent human
reviewers. Analysis (agreement, Cohen's κ, error rate of the extraction) runs only after both columns are filled
by people; nothing is imputed.

---

## Outputs
Scripts `src/87_tool_comparison.py`, `src/88_real_mags.py`, `src/89_conformal.py`, `src/90_habitat_case.py`,
`src/91_label_check_sheet.py`; results in `results/additional/`; report `docs/RESULTS_additional.md`.

## Deviations
(none yet)
