# Pre-registration BL: blind test of the frozen final model on recently described species

Registered on 22 October 2026, 12:51 MSK, before any candidate species was searched, extracted or scored. The SHA-256
of this file is stored in `docs/preregistration_BL_blind_test.sha256`.

## Why

The final donor and acceptor model RedoxKO (random forest) was chosen on 21 October 2026, after the tests on the
external genomes and on the 143 newly described species had been run, and partly using them. Its results on those
sets are therefore not blind. This test scores the frozen model on species that played no part in any choice: not in
training, not in the external set, not in any earlier new-species set, and not seen by the authors before this
registration.

## Frozen model (no change after this registration)

| File | SHA-256 |
|---|---|
| results/models/RF/model.joblib | 310d84ef335579bace4d5118af33c8e46a9fe2810327171742e2177227b6362a |
| results/models/RF/metadata.json (thresholds, release rule with four conditions) | 624b1ee156fae2a07796aa0b6508b6778d16a6b2fa64f371fd96d6f0131cce6e |
| src/predict.py | f0bb8b53df71da6a52b2350bab28926a3d7c4e09729355099b9093718559c7af |

RedoxKO-LX (results/models/B2), unchanged since 17 September 2026, is scored as a secondary model.

## Candidates

- PubMed species-description papers with "nov" in the title and the same donor/acceptor terms as
  `src/68_harvest_candidates.py`, with PubMed entry date from 15 September 2026 to 22 October 2026. The earlier
  new-species search was run on about 19 September 2026, so this window covers papers that could not have been in it
  or appeared only after it; genus exclusion below removes any overlap.
- Excluded: every genus, species and genome accession in `validation_new_organisms/exclusions.json`,
  `validation_set.json`, `pool/`, `pool_auto/` and `pool100/`, i.e. everything in the training, external, pathway and
  earlier new-species data. One species per genus. An annotated RefSeq genome is required.

## Extraction

Growth tests are extracted from the description papers by independent LLM agents with the rules of
`validation_new_organisms/pool/SCHEMA.md` (verbatim quote and source for every test; never infer from genes), without
access to gene content or any model output. Quotes are checked against the source text by script. Statements that
are not growth experiments (stimulated growth, medium components, fermentative reduction, genome-based statements)
are set to not_tested with the rules of `src/72_fix_evidence.py`. The two catch-all categories are excluded, as in the
primary analysis.

## Scoring and endpoints

Each genome is scored once with the frozen model from its RefSeq accession (`src/predict.py`).

- Primary: macro AUROC of RedoxKO over the categories that have both positive and negative tests, with a 95% cluster
  bootstrap interval over description papers (2,000 resamples, seed 20261022).
- Secondary: accuracy and MCC of RedoxKO's released yes/no calls; the same values for RedoxKO-LX; the number of tests
  per category.
- Expected: macro AUROC between 0.60 and 0.75 (about 65%), i.e. similar to the 143 species (0.73) and the 43 withheld
  species (0.62). If fewer than two categories have both positive and negative tests, or fewer than 20 tests in total,
  the result is reported descriptively only.

Results are reported regardless of their outcome, in `docs/RESULTS_BL_blind_test.md`.

## Repeat

The same procedure, with the same frozen files, will be repeated on species whose description papers enter PubMed
after 22 October 2026, once at least six months have passed (from 22 April 2027).

## Amendment A1 (22 October 2026, 12:52 MSK, before any extraction or scoring; time corrected after the harvest from a typo of 13:10)

The registered window (PubMed entry 15 September to 22 October 2026) returned 3 description papers and no species
that passed the exclusions; nothing was extracted or scored. Because blindness requires only that a species played
no part in any choice, not a late publication date, the candidate search is widened:

- Query: "nov" in the title and any of nitrate, anaerob*, aerobic, electron acceptor, electron donor, thiosulfate,
  sulfate, Fe(III), methanogen*, hydrogen in title or abstract; PubMed entry date 1 January to 22 October 2026.
- The same exclusions as above (every genus, species and accession in any earlier set or pool), one species per
  genus, annotated RefSeq genome required.
- The 60 newest candidates (by PubMed ID) are extracted. Candidates from the autotroph pools are not used, because
  a gene-based filter selected them.
- All other rules (extraction, frozen model, endpoints, expectations) are unchanged.

## Amendment A2 (22 October 2026, after extraction and the first audit run, before any scoring)

No model output had been produced for any candidate when this amendment was made. The first run of the audit in
`src/176_blind_test.py` showed two faults in the script, not in the registered rules:

- Quote matching. Quotes taken from Europe PMC full texts named in an evidence record were not found, because the
  script read only the full text linked in the harvest record and because tag stripping split strain names
  (e.g. "AT8<sup>T</sup>" became "at8 t"). The audit now also reads every Europe PMC full text named in the record and
  compares text with all whitespace removed.
- Demotion rule. The word "ferment" also matched "isolated from fermented fish" and API-strip result lists. It now
  matches only fermentative reduction, as in `src/72_fix_evidence.py`. As in that file, a nitrate-reduction result
  from an API strip or kit test without growth coupling is a biochemical test and is set to not_tested.

The scoring step also reports the AUROC of each category that has both positive and negative tests, so that the result
can be described when the primary endpoint is not defined. All other rules are unchanged.
