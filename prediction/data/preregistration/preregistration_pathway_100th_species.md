# Pre-registration: a 100th newly described species for the pathway test set

Registered 18 October 2026, before the PubMed search below was run and before any new candidate was seen. SHA-256 in
`docs/preregistration_pathway_100th_species.sha256`. Released models are not changed.

## Why
The pathway test set holds 99 species: every species from the original search (run 19-20 September 2026) that had an
annotated RefSeq genome, a blind growth test mappable exactly onto a pathway trait, and a genus of its own. A 100th
species can only come from new data. This document fixes how it is found before anything is looked at.

## What was known before registration
All results on the 99 species (AUROC 0.67-0.70 for the pathway models). Nothing about any new candidate.

## Search (identical to the original)
`src/68_harvest_candidates.py` functions, unchanged: the same PubMed query (`nov[ti]` + the same chemistry terms,
publication years 2015:2026), binomials parsed from titles, an annotated RefSeq genome required (same preference
order). Output goes to `validation_new_organisms/pool100/`, so the original pool is untouched.

## Exclusions
Every genus, species and accession in `exclusions.json` (training, external and pathway data), plus every genus and
species in the original candidate pool (`pool/pool.json`, 178 species) and in `validation_set.json` (167 species).
This keeps one species per genus across the whole new-species collection.

## Selection rule (one species, fixed now)
1. Order the remaining candidates at random with seed 20261018.
2. In that order, extract the growth tests of each candidate from its description paper, following
   `pool/SCHEMA.md` exactly (verbatim quote and source for each test; full text where open access, tier A, otherwise
   the abstract, tier B; never infer from genome content), without running any model on its genome.
3. The first candidate with at least one test that maps exactly onto a pathway trait (`EXACT_PW` in
   `src/78_balanced_subsets.py`: sulfate, sulfur and nitrate reduction, sulfur, iron and hydrogen oxidation, iron
   reduction) and status `tested_positive` or `tested_negative` becomes species 100. Extraction stops there.
4. If no candidate qualifies, the set stays at 99 and the text says "all 99 eligible species".

## Scoring
The 100th genome is annotated with the same pipeline (`src/79_pathway_validation_predict.py`), its Bacformer-large
vector is computed with `src/74_local_embed.py`, and all saved pathway models are scored on the 100 species with
`src/84_novel_pathway_comparison.py`, unchanged. The paper reports the 100-species results; the 99-species values stay
in the repository.

## Reporting
The paper states that species 100 was added after the original selection, by this rule, from a search run on
18 October 2026.
