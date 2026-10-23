# Pre-registration BL2: blind test of the frozen final model on lithotrophic species (donor categories)

Registered on 22 October 2026, before any candidate species was searched, extracted or scored. The SHA-256 of this
file is stored in `docs/preregistration_BL2_donor_blind_test.sha256`.

## Why

**Gap left by BL.** Blind test BL (`docs/RESULTS_BL_blind_test.md`) gave tests almost only for O2 respiration, because recent species descriptions rarely test lithotrophic electron donors. Reviewers asked for a blind test of the donor categories that motivate the paper.

**How this test closes it.** It searches description papers of lithotrophs from all years. It keeps only species from genera that occur in none of the project's data and that played no part in any choice.

## Frozen model (unchanged since BL)

| File | SHA-256 |
|---|---|
| results/models/RF/model.joblib | 310d84ef335579bace4d5118af33c8e46a9fe2810327171742e2177227b6362a |
| results/models/RF/metadata.json | 624b1ee156fae2a07796aa0b6508b6778d16a6b2fa64f371fd96d6f0131cce6e |
| src/predict.py | f0bb8b53df71da6a52b2350bab28926a3d7c4e09729355099b9093718559c7af |

RedoxKO-LX (results/models/B2) is scored as a secondary model.

## Candidates

**PubMed query.** "nov" in the title, and in the title or abstract any of:
- sulfur-oxidizing, sulphur-oxidizing, thiosulfate, sulfide-oxidizing;
- hydrogen-oxidizing, "hydrogen oxidation", "Fe(II)", iron-oxidizing;
- ammonia-oxidizing, nitrite-oxidizing, "carbon monoxide";
- chemolithoautotroph*, chemolithotroph*, lithoautotroph*, "electron donor".

There is no date limit; the latest entry date is 22 October 2026.

**Exclusions.**
- Every genus, species and accession in the earlier sets: `exclusions.json`, `validation_set.json`, `pool/`, `pool_auto/`, `pool100/`, `candidates.json` and `pool_blind/`.
- Candidates from the autotroph pools are not used, because a gene-based filter selected them.

**Selection.**
- One species per genus.
- An annotated RefSeq genome is required.
- The 80 newest candidates (by PubMed ID) that pass these rules are extracted.

## Extraction and audit

As in BL, after amendment A2:
- Independent LLM agents extract tests with the rules of `validation_new_organisms/pool/SCHEMA.md`. They have no access to gene content or model output.
- Every quote is checked verbatim against an open source.
- Statements that are not growth experiments are set to not_tested.
- Nitrate results from API-strip or kit tests without growth coupling are set to not_tested.
- The catch-all categories are excluded.
- The labels are frozen with a SHA-256 hash before scoring.

## Endpoints

Each genome is scored once with the frozen model from its RefSeq accession.

**Primary.** Macro AUROC of RedoxKO over the six donor categories, counting those that have both positive and negative tests. The 95% interval comes from a cluster bootstrap over description papers (2,000 resamples, seed 20261022).

**Secondary.**
- Macro AUROC over all 12 categories with both classes.
- Per-category AUROC and test counts.
- Accuracy and MCC of released yes/no calls.
- The same values for RedoxKO-LX.

**Expected.** Donor macro AUROC between 0.60 and 0.85 (about 60%). On the 143 earlier new species, the donor categories ranked about as well as the acceptor categories, except Fe(II) oxidation.

**Rule for too few tests.** If fewer than two donor categories have both classes, or there are fewer than 20 donor tests, the primary result is reported descriptively only.

Results are reported whatever the outcome, in `docs/RESULTS_BL2_donor_blind_test.md`.

## Amendment A1 (22 October 2026, after the harvest, before any extraction or scoring)

The registered search returned 852 description papers. Only 8 species passed the exclusions and had an annotated
RefSeq genome, because the earlier sets already cover nearly every lithotroph genus. That is too few for the primary
endpoint. Nothing has been extracted or scored. A second tier is added and reported separately:

- **Tier 1 (as registered).** Genus absent from all project data. Scope: the 8 species already harvested.
- **Tier 2.** Genus absent from the 358 training genomes of RedoxKO, but present in other project data (external, pathway or earlier new-species sets). This is the same distance to training as the external genomes.
  - The species and its genome accession must be absent from every earlier set.
  - One species per genus.
  - An annotated RefSeq genome is required.
  - Same query; the 80 newest candidates are extracted.

**How tier 2 is reported.**
- It is not blind to genus-level knowledge in the label sources.
- It is blind to model choice, because no tier-2 species was seen before this registration.
- The primary endpoint is computed on tiers 1 and 2 combined, and is also reported for each tier.
- All other rules are unchanged.

## Amendment A2 (22 October 2026, during extraction, before the audit and before any scoring)

The audit script exempted any statement containing the word "growth" from the not-a-growth-test rule. As a result,
statements such as "elemental sulfur stimulated growth" or "thiosulfate enhanced cell growth" would have been kept,
although `src/72_fix_evidence.py` and the schema treat stimulated growth as not a growth test. For BL2 the audit now
sets any tested status whose quote reports stimulated or enhanced growth to not_tested. BL was already scored and is
not changed.
