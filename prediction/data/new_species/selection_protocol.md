# Scaled prospective validation — protocol fixed before any prediction (written 2026-09-21, ~02:30 JST)

No model has been run on any organism in this set at the time of writing. The first round's protocol
(selection.md, 10 organisms) is kept; only the sampling is widened and the hand-picking is removed.

## 1. Candidate harvest (automated, fixed query)
PubMed, species-description papers 2015-2026: `nov[ti]` AND a fixed list of donor/acceptor chemistry terms
("electron acceptor", "electron donor", thiosulfate, "elemental sulfur", "sulfate-reducing", "sulfur-oxidizing",
"iron-reducing", "nitrate reduction", methanogen*, chemolithoautotroph*, chemolithotroph*, "hydrogen oxidation",
acetogen*, "Fe(III)") -> 553 papers -> 194 proposed binomials after removing every genus, species and accession
that appears anywhere in the project's data (exclusions.json: 1 746 genera, 3 776 species, 3 695 accessions)
-> 178 with an annotated RefSeq genome (src/68_harvest_candidates.py).

## 2. Evidence extraction (blind)
Six independent agents read the papers (abstract, PMC full text where open access, LPSN, BacDive) with no access to
genome content, model output or predictions, and recorded for each of the 14 categories:
tested_positive / tested_negative / not_tested, with a verbatim quote and its source (pool/SCHEMA.md).
Strict rules: assimilatory use does not count; "strictly anaerobic" without an oxic growth test is not_tested;
a capability stated only from genes is not_tested; ambiguous wording ("potential electron acceptors",
"stimulated growth", a nitrate-reduction biochemical test without growth coupling) is not_tested.

## 3. Inclusion rule (no curation)
**Every organism with at least one tested_positive or tested_negative is included.** Nothing is dropped for being
inconvenient, unbalanced or hard; the first round's 10 hand-curated organisms are merged in and de-duplicated by
accession, keeping the richer record. Balance is reported (validation_balance.md), not engineered.

## 4. Models (both trained on all 358 labelled genomes, unchanged)
* **B2** — released model, 852-KO panel.
* **B2XLE_large** — B2XL + Bacformer-large probe; genome vectors computed on Kaggle exactly as in training.
Binary calls at the stored MCC cut-offs; scores and calibrated scores reported alongside.

## 5. Scoring rule (fixed now)
* Only categories with a published test are scored: tested_positive -> expected 1, tested_negative -> expected 0.
* Reported per model: scored calls correct, sensitivity on tested positives, specificity on tested negatives,
  the same restricted to the 8 categories that get released yes/no calls ("deployed mode"), and per-category
  counts with Wilson confidence intervals.
* Every mismatch is classified against the organism's own marker genes and further sources as model error,
  feature gap (no KEGG orthology for the chemistry), genotype-phenotype gap (gene present, growth not observed),
  unreliable category, or literature ambiguity — the classification used in the first round.
* Both models are scored on exactly the same calls; the paired difference is reported with a bootstrap CI.
