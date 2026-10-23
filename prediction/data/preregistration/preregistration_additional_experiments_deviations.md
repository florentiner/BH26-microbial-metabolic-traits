# Deviations from `preregistration_additional_experiments.md` (SHA-256 d66d2c1d…36d3ae)

Kept in a separate file so that the registered document and its hash stay unchanged.

| # | Date | Experiment | Deviation | Reason / effect |
|---|---|---|---|---|
| D1 | 2026-10-14 | X1 | FAPROTAX 1.2.12 has no group named `acetogenesis`; its group `reductive_acetogenesis` was mapped to acceptor_CO2_methanogen_acetogen instead. | Same biology under a different name; the literal name in the protocol does not exist. |
| D2 | 2026-10-14 | X1 | BacDive-AI's documented input is InterProScan Pfam output; InterProScan is not installable locally. Pfam matches were produced with pyhmmer `hmmsearch` on Pfam-A 38.2 with gathering cut-offs, keeping families with a domain independent E-value ≤ 1e-20 (the tool's own threshold), restricted to the 1 586 families the aerobic model uses. | Features may differ slightly from InterProScan's (Pfam release and E-value definition). Reported as an approximation of the tool. |
| D3 | 2026-10-14 | X2 | Paired isolates are predicted from nucleotide FASTA through the same Pyrodigal path as the MAGs, in addition to the existing annotated-protein predictions. | Keeps gene calling identical on both sides, so differences measure assembly incompleteness, not annotation pipelines. |
| D4 | 2026-10-14 | X1 | FAPROTAX's own `collapse_table.py` was not executed (third-party code); matching uses the pre-registered fallback (wildcard matching of each taxon pattern against the lineage string), including `add_group` / `subtract_group`. | Pre-registered fallback; no effect on the mapping. |
| D5 | 2026-10-14 | X4 | The `marine_surface` keyword "ocean water" also matched deep-water samples (e.g. abyssopelagic) and "seawater" matched aquaculture seawater. The class is reported as "seawater" rather than "surface". | Keyword lists were fixed in advance and are kept; the label of the class is corrected in reporting. |
| D6 | 2026-10-14 | — | Added analysis outside the protocol: impact of the 2026-10-05 audit of the phenotype tables on the training labels (`src/92_label_audit_impact.py`). | Requested after the protocol was written; reported separately as exploratory. |
| D7 | 2026-10-14 | Add. A | METABOLIC-G run in its `full` KOfam mode with a KO list (`prokaryote.hal`) restricted to the union of the KOs in its function table and its module list (2 723 KOs; 2 722 exist in KOfam 2026-08-25), instead of all prokaryotic KOs. | KOs outside these lists never enter METABOLIC's function or module calls, so function presence is unchanged; it only saves ~6 GB of pressed HMMs. METABOLIC's own "small" mode would have skipped 79 function-table KOs and was not used. |
| D8 | 2026-10-14 | Add. A | METABOLIC's function summary (worksheet 2) omits its iron-cycling functions; these are read from its per-HMM sheet (worksheet 1), present when any HMM of the function is present. Functions are matched by name because the two sheets label the metal-reduction category differently. | Needed to score the iron categories at all. |
| D9 | 2026-10-14 | Add. A | KO K18513 is referenced by METABOLIC but absent from KOfam 2026-08-25; its search returns nothing. | One KO; treated as absent. |
| D10 | 2026-10-15 | X1 / D2 | InterProScan 6.0.2.2 (Docker, Pfam only, InterPro data fetched by the pipeline) was run to remove the D2 approximation. Its Pfam search container (`interpro/hmmer:3.3`, x86-64) was killed under emulation on Apple silicon; the same image was rebuilt natively for arm64 from EBI's Dockerfile with HMMER 3.4 instead of 3.3 (3.3 has no ARM vector code). Most proteins received precomputed matches from EBI's matches service, which InterProScan uses by default. | Same profiles, cut-offs and post-processing; HMMER 3.4 uses the same scoring model as 3.3. |
| D11 | 2026-10-15 | X2 | The released pathway model (PB3) on the MAG pairs is run in addition to B2 (not in the protocol), locally with the training KOfam release; Kaggle could not be used because the current KOfam release changes thresholds of 6 064 of its 10 517 features. | Extension; reported as such. |
| D12 | 2026-10-15 | D | The Qwen 3.8 Flash confirmation step was dropped before any training organism was submitted; final labels = Claude Haiku 5.5 claims whose quote occurs verbatim in the supplied text. | On the 19 pilot organisms read by both models, Qwen agreed with every Haiku positive/negative call (identical κ, precision, recall), so the check had no measurable effect on the pilot; Haiku's claims and quotes are stored so the check can be added later. Decided by the authors before the run. |
| D13 | 2026-10-15 | D | Retraining on `curated_D` was stopped before it ran. At 54 of 355 species, 98 % of former zeros stayed "unknown" (1 new positive, 7 tested negatives); masking unknowns as pre-registered would remove nearly all negatives and leave a degenerate training set. No curated model was fitted and no test set was scored. | The literature rarely reports tested negatives for donor/acceptor use, so the premise of D (real negatives are recoverable) does not hold. Extraction completes for descriptive statistics only (< $1). |

## 17 October 2026: A, B-clean, D, H
* B-clean: no deviation. All 546 added genomes were annotated; none failed.
* B-clean, reporting file: results are in `docs/RESULTS_Bclean_A_D.md` rather than the `docs/RESULTS_Bclean_A.md`
  named in the pre-registration, so that A, B-clean, D and H are reported together.
* D, H: no deviation.
* Figures revised after these results (no analysis changed): Fig. 2 now shows the released RedoxKO's attributions
  (`src/113_xai_redoxko.py`) beside the random forest; pathway enrichment and the label-noise simulation moved to the
  supplement; Fig. 3 compares RedoxKO with METABOLIC instead of RedoxKO-ens-BF.

## 18 October 2026: F, M
* F: no deviation; 402 genomes annotated, none failed.
* M: no deviation. M4 timing for METABOLIC is read from `full_run.log` (100 genomes) and `holdout_run.log` (43
  genomes), as registered; the RedoxKO times are those already reported in the paper. M3 used all 379 pairs of
  `x2_pair_calls.csv` (594 genomes, none failed). MAG FASTA files were deleted after gene calling; proteins are kept in
  `data/additional/mags/faa/`.

## 18 October 2026: search for a 100th pathway species
Registered in `docs/preregistration_pathway_100th_species.md` (SHA-256 37d4fd58...) before the search. The repeated
PubMed search returned 556 description papers (553 in September); after the exclusions 15 candidate species remained,
of which one had an annotated RefSeq genome (*Hydrogenimesophilus fermentans*, GCF_038792415.1, IJSEM 2026). Its full
text is not open access, and the abstract reports no test of a pathway trait (`validation_new_organisms/pool100/
evidence_100th.json`), so by rule 4 the pathway set stays at all 99 eligible species. No model was run on the genome.

## 18 October 2026: correction to the F pre-registration
`docs/preregistration_F_full_ko_profile.md` states that the TogoMCP KO matrix showed "6,086 local-only hits against
3 TogoMCP-only hits", read as a systematically sparser annotation. That reading was wrong. All 6,085 local-only hits
come from 52 training genomes absent from the shared resource (their rows were empty). On the 306 genomes present in
both, the local search and the resource agree in 99.95% of 260,712 genome-KO cells (132 local-only, 3 resource-only),
consistent with `results/qc/togo_concordance.json`. Experiment F is unaffected: it used the local annotation for every
genome, which was necessary anyway because 52 genomes are missing from the resource.

## 19-21 October 2026: BT, CB, AT, analyses outside the protocols, and the change of final model
* BT (registered 19 October), CB (registered 20 October) and AT (registered 21 October, before harvesting): run as
  registered. AT: the candidates ran out at 99 autotroph-like species, which the protocol allows.
* Not registered (exploratory, reported as such): pathway models trained on autotroph-like genomes only
  (`src/135_pathway_autotroph.py`); external values for the Table 1 rows that had none (`src/138_external_missing_rows.py`);
  comparison of the random forest with the final model on every test (`src/139`-`141`).
* **Change of final model (21 October, after all tests).** The random forest on the 852-KO panel replaced B2 as the
  final donor and acceptor model RedoxKO; B2 is renamed RedoxKO-LX. The forest was better on the external genomes
  (+0.043 macro AUPRC [0.018, 0.069], also without FAPROTAX-derived labels), in genus-level CV (+0.037 [-0.001, +0.059])
  and in yes/no calls on new species (MCC +0.056 [0.006, 0.115]). Post hoc: decided after the results were known.
  All analyses of the final model were repeated with the forest (`src/143`-`155`); these repeats were not registered.
  The registered experiments keep their RedoxKO-LX results.
