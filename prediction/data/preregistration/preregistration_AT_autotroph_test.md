# Pre-registration AT: test of the donor and acceptor models on 100 autotroph-like new species

Registered 21 October 2026 before any new species was harvested or scored. SHA-256 in
`docs/preregistration_AT_autotroph_test.sha256`.

**Question.** The training collection consists of autotrophs (80% autotroph-like by the genome proxy below), but only
46 of the 100 primary new species are autotroph-like. How do RedoxKO and RedoxKO-ens-BF perform on newly described
autotroph-like species, with enough species for a usable estimate?

**Autotroph-like** (unchanged from `src/133_autotroph_subset.py`): at least one of the KEGG modules M00165, M00173,
M00377, M00376, M00375, M00374 or M00567 is at least 80% complete in the 852-KO panel.

**Test set.** (1) The autotroph-like species among the 167 already scored new species (81 known so far). (2) New
species harvested with the same procedure as `src/68_harvest_candidates.py` (PubMed `nov[ti]`, description papers,
same exclusion of every genus, species and accession in the project data, annotated RefSeq genome, one species per
genus, genus not among the existing new species), with a broader query that adds autotrophy and carbon-fixation terms
and the years 2010-2014. Growth evidence is extracted blind to any model output by agents following
`validation_new_organisms/pool/SCHEMA.md`, with a quoted string for every call. New species are added, in harvest
order, until 100 autotroph-like species with at least one growth experiment are reached or the candidates run out.

**Endpoints.** Primary: macro AUROC of RedoxKO over the categories with both outcomes, with a 95% bootstrap interval
over papers. Secondary: pooled MCC, sensitivity and specificity; the same for RedoxKO-ens-BF; the difference between
autotroph-like and other species among the 167 already scored species.

**Expected.** Macro AUROC on autotroph-like species at least as high as on all primary species (0.73) (about 70%),
because the training set consists of autotrophs; MCC lower than on all primary species (about 60%), because
autotroph genomes carry more energy-metabolism genes and receive more positive calls.

**Reporting.** Reported whatever the outcome in a supplementary table and in one sentence of the Results.
