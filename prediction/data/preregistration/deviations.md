# Deviation and amendment log

Every change relative to `docs/preregistration.md` and to `redox-model-search-protocol.html`.
Times are local (JST) on 2026-09-17. The project was put under git at commit `1b678a8`
("Snapshot before audit round 1 fixes"). Before that the only evidence of order is file modification
times, listed here so readers can judge them.

## Timeline before git
| time | event |
|---|---|
| 18:06 | labels + frozen splits written (`splits.sha256` = ae0cdef1…) |
| 18:07:52 | final pre-registration text saved (includes arm B1r) |
| 18:09 | nulls and a first B1 run scored (B1 re-run later with identical logic) |
| 18:15–18:17 | refactors, no change to any model, grid, rule or metric: tabular fitting moved to `redox_common.fit_models`; KOfam/ESM/Bacformer steps split into importable functions; `RULE_GROUPS` unchanged since it was created (~18:08) |
| 18:42–18:44 | arms A, B1, B1r, B2 scored (first scoring of A, B1r, B2) |

`04_kofam_scan.py` and `03_embed_esm2.py` were refactored while their jobs were running. So:
* KOfam: re-run with the frozen code; the result must be identical to the first run (check in `results/qc/reproducibility_checks.json`).
* ESM-2: a random sample of stored embeddings is recomputed with the current code and compared (same file).

## Amendments to the pre-registration (after audit round 1, before any embedding arm was scored)
1. **Target definition.** Labels are *database annotations*, not verified phenotypes. The audit found likely
   label errors among the closed-world negatives, e.g. *Methanocaldococcus* genomes with mcrABG labelled
   non-methanogen, and ammonia oxidisers with `acceptor_oxygen = 0`. Added:
   (a) sensitivity analysis S1: when a genome's whole donor (or acceptor) block is all-zero, that block is
   treated as unobserved (NaN) instead of negative;
   (b) a candidate-label-error table (`results/candidate_label_errors.csv`): high-scoring out-of-fold
   "false positives" with their marker-gene evidence.
2. **Taxonomy.** The frozen splits use NCBI taxonomy from the provided table, not GTDB as protocol gate G3 asks.
   The NCBI splits stay primary (they were frozen before scoring). Added sensitivity S2: genus and family folds
   rebuilt from GTDB R232 taxonomy.
3. **Partition variance.** Only one partition per level is frozen, so genus vs family differences are
   confounded with fold-assignment noise. Added sensitivity S3: 5 further fold allocations (seeds 1–5) at both levels.
4. **Genome quality.** Added sensitivity S4: evaluation restricted to genomes with CheckM2 completeness ≥ 90 %
   and contamination ≤ 5 % (GTDB metadata). Arm E is additionally run with the `bacformer-masked-MAG`
   checkpoint (arm E_MAG), because 174/358 labelled assemblies are not complete genomes.
5. **Multiplicity.** Per-trait bootstrap intervals are descriptive. Per-trait tests against the null and the
   paired arm contrasts get two-sided bootstrap p-values with Benjamini–Hochberg q-values within each
   (contrast, split level) family. Primary inferential claims use macro-averaged contrasts only.
6. **Release rule clarified.** The pre-registered rule "arm with highest genus-level macro AUPRC" is limited to
   arms computable from a genome sequence. B1 consumes the provided marker table and cannot be applied to a
   new genome. Arms whose macro-AUPRC paired-contrast CI includes 0 are reported as **tied**, and the release
   choice among tied arms is presented as a tie, not a win.
7. **Arm A threshold.** The registration table says θ ∈ {0.5, 1.0} "on inner CV"; the text and code use
   θ ∈ {>0, 0.5, 1.0}, selected by MCC on the outer training fold. Arm A fits nothing, so inner CV and the full
   training fold would be equivalent up to pooling. The code is authoritative; the table was a typo.
8. **Arm A rule groups.** Non-specific markers are left out of the rules: cynS, rbcL/prkB, nirBD
   (K00362/K00363), hdrABC, Rnf, pta; cooS is kept only under CO oxidation. This was fixed in code before arm A
   was scored but was not written into the pre-registration text.
9. **Arm D.** Hyper-parameter grid reduced to prefilter {off, on} x trunk 128 x instance dropout 0.2 x entropy 0
   (protocol asks trunk {64,128}, dropout {0,0.2,0.3}, entropy {0,0.01,0.05}) for local compute. The selected config
   was prefilter=on in 7/10 fold-levels, which injects curated KO knowledge; arm **D_nopf** (prefilter forced off) is
   therefore run separately as the clean attention-vs-mean-pooling contrast for the protocol's Q3.
   Thresholds are selected on inner out-of-fold predictions of the same 5-seed ensemble that
   is applied to the test fold (the registered text used 1 seed for selection). The prefilter uses the
   energy-metabolism KO panel, so "prefilter on" also injects annotation knowledge (D vs C is then not a pure
   aggregation contrast; both D configurations are reported).

10. **Secondary macro over 12 categories** (added after first results, 2026-09-17): `macro12_*` columns exclude
    `donor_other` and `acceptor_other`, heterogeneous catch-all categories (every genus with a positive has mixed
    labels; results/label_ceiling_summary.csv). The pre-registered 14-trait macro stays the primary metric.
11. **Research loop and external validation** (exploratory, after the primary analysis): docs/research_loop.md,
    docs/external_validation.md; all variants in results/research_ledger.csv.

12. ~~**Release rule after audit round 2.**~~ **SUPERSEDED by #16** (kept for the record). Released model = the simplest model in the set statistically tied with the
    best variant (paired genus-cluster bootstrap CI includes 0; simultaneous max-T CIs in
    results/simultaneous_inference.csv). This is **B2XL**. Supersedes #6 and the research-loop promotion of B2E
    (docs/research_loop.md). Pre-registered-rule choice (B2) and all tied models are reported.
13. **External validation corrections (audit round 2).** The `nitrate_strict` subset (acceptor_nitrate positives
    restricted to denitrification / nitrate_reduction_to_ammonia, i.e. excluding possibly assimilatory
    `nitrate_reduction`) was added at the same time and is declared here. Novelty is lineage-based: 20 of 582 genomes removed (17 were
    among the internal 586, owing to genus renames since Madin et al. 2020). The taxonomic null now uses current NCBI
    lineages. Paired contrasts were added. FAPROTAX-sourced positives are reported per trait: for acceptor_iron,
    donor_hydrogen_oxidation and donor_iron_oxidation, all external positives come from FAPROTAX-sourced species.
14. **LLM baseline provenance limits (audit round 2, M8).** The Claude arms ran as isolated sub-agents of this analysis
    session (model: Claude Opus 5; tools limited by instruction to reading one batch file and writing one output
    file). Prompts are the batch files; per-agent transcripts, token usage and a hard sandbox are not archived.
    These arms are therefore reported as *non-reproducible* comparisons. The GPT-6 Astra arms (OpenRouter, model and
    usage logged per batch) are reproducible up to provider-side model versioning.
15. **predict.py on macOS.** Native runtime order is enforced (XGBoost unpickling, then PyTorch, then pyhmmer); other
    orders segfault. Tested end-to-end on a training genome (scores equal those from precomputed features) and a new
    genome.

16. **Confirmatory external test for release (added 2026-09-18, audit round 2 recommendation).** A research-loop
    variant replaces the pre-registered arm only if it (a) improves in internal nested CV *and* (b) does not lose on
    the external novel-genus set in any reported subset. Internal winners B2XL (+0.026 genus) and ENS (+0.013) both
    failed (b): external B2XL − B2 = −0.007 [−0.022, +0.009] overall and −0.050 [−0.080, −0.014] without
    FAPROTAX-sourced labels; ENS − B2 = −0.008 [−0.022, +0.007]. **Released model: the pre-registered B2.** The
    internal gains of learner ensembling are reported as not transferring to independently curated labels.

17. **Second task added: energy-metabolism pathway prediction (2026-09-18).** The project brief asks for electron
    donors, electron acceptors *and* energy-metabolism pathways. The donor/acceptor labels are organism-level
    categories and cannot supply pathway labels, so a second, separately pre-registered task was added
    (`docs/preregistration_pathways.md`, written and hashed before any pathway model was fitted or scored:
    `93c7a99f2c46820f58dbd8a30e8c4ab633c3992afc40de5be751e34e0f53ea32`). It has its own labels (Madin et al. 2020
    `pathways` field), its own genomes (one RefSeq genome per species with a BH26 KO profile), its own frozen
    genus/family-blocked splits, and its own arms (PA rule / PB1 / PB2 / PB3 / PE). Nothing in the donor/acceptor
    analysis was re-run, re-selected or re-reported because of it. Scripts: `src/30_pathway_dataset.py`,
    `src/32_pathway_arms.py`, `src/33_pathway_metrics.py`, `src/34_pathway_mag_robustness.py`,
    `src/35_pathway_train_final.py`. The KEGG module-completeness readout (`src/29_pathway_modules.py`), previously
    the only pathway output, becomes the rule baseline (arm PA) for this task.
18. **`predict.py` now emits all three outputs from one genome input (2026-09-18).** Donor/acceptor categories and
    pathway probabilities come from trained models; module completeness remains available as deterministic evidence
    (`--pathways model|modules|both|none`). The pathway model scans the KO profiles its own feature columns require,
    which is slower than the 852-KO donor/acceptor panel; this cost is reported rather than avoided by silently
    substituting a smaller feature set.

19. **Pathway label sensitivities, declared before any pathway arm was scored (2026-09-18).** The pathway labels are
    closed-world: a species whose Madin record does not list a pathway is treated as a negative for it. Two
    sensitivity analyses bound what that assumption can buy, both specified before scoring:
    * **SP1 `no_faprotax`** — evaluation restricted to genomes whose Madin record is not FAPROTAX-derived. FAPROTAX
      propagates traits by taxonomy, so FAPROTAX-sourced labels are partly circular with any taxonomically
      structured feature; this subset removes that circularity from the evaluation set (models are unchanged).
    * **SP2 `informative_negatives`** — evaluation restricted to genomes with at least two annotated pathway terms,
      i.e. organisms characterised well enough that the absence of a third term is weak evidence of absence.
    Primary results use the full closed-world set; both sensitivities are reported alongside, and a claim that
    reverses under either is not made.

    **Measured before any pathway arm was scored (label side only).** Of 3 378 genomes, 1 176 have a Madin record
    that does not cite FAPROTAX (SP1) and 998 carry >= 2 annotated pathway terms (SP2). SP1 therefore bounds the
    circularity for only **5 of 20 traits** with >= 10 positives left (nitrate reduction 957, methanogenesis 100,
    sulfur reduction 40, nitrogen fixation 34, sulfate reduction 21); for the other 15 traits essentially every
    positive comes from a record that cites FAPROTAX, and no non-taxonomic bound on those labels is available from
    this source. SP2 keeps all 20 traits. This is a limitation of the only large literature-curated pathway label
    set available, and it is reported as such rather than hidden: the honest statement a reviewer needs is that
    *for 15 of 20 pathway traits the label provenance cannot be separated from FAPROTAX's taxonomy-based
    propagation*, and that the genus- and family-blocked design plus the taxonomic group-mean null are what
    quantify how much of the measured performance is lineage identity rather than gene content.
    `faprotax_source` is a conservative flag: it marks a record citing FAPROTAX among its sources for *any* trait,
    not specifically for the `pathways` field.

20. **Pathway splits re-frozen once, before any model was scored (2026-09-18 13:55).** The first freeze left 1 genome
    with no `group_genus` and 10 with no `group_family` (species absent from GTDB R232 and lacking an NCBI genus or
    family: *Leptothrix ochracea*, *Sphaerotilus natans*, *[Desulfobacterium] catecholicum* and 7 others). A missing
    group is not a cosmetic defect — the grouped splitter would have had to order mixed types, and a NaN group would
    silently merge unrelated genomes into one split unit. Every learned arm failed loudly on it, so **no pathway
    model was ever fitted or scored under the first freeze**; the prevalence/taxonomic nulls and the rule arm PA had
    produced output and were deleted unread. Fix: a genome with no resolvable rank becomes its own singleton group,
    which can never share a fold with another genome. Labels re-frozen:
    `c866e2fe236ed408e0a3a68517e0c4b7e4aff48c240c66e341e10b678ed14492` (was `760398fe58ae…`). Trait counts,
    genome set and all pre-registered decisions are unchanged.
21. **External check made symmetric across pre-registered arms (2026-09-18 14:40).** When the large-encoder arms were
    scored, `E_large` (Bacformer-large) came top of the pre-registered ranking on internal genus-level macro AUPRC by
    **0.0013** (0.5769 vs 0.5756 for B2) — a difference far inside the bootstrap CIs ([0.534, 0.623] vs
    [0.545, 0.630]) — and the release code, which took the internal point-estimate leader, wrote `E_large` into
    `RELEASED.json`. On the external novel-genus set with independently curated labels, `E_large` is *significantly
    worse* than B2: −0.071 [−0.109, −0.029] overall and −0.265 [−0.416, −0.106] on the nitrate-strict subset (tied,
    +0.019 [−0.038, +0.079], without FAPROTAX-sourced labels).
    Deviation #16 already required a research-loop variant to confirm externally before displacing the pre-registered
    arm; there is no principled reason to exempt a pre-registered arm from the same test. The rule is therefore
    stated symmetrically: **the released arm is the highest-ranked pre-registered arm on internal genus-level macro
    AUPRC that is not significantly worse than another pre-registered arm on the external set**, and a research-loop
    variant displaces it only if it also confirms externally. Under this rule `E_large` is excluded and **B2 is
    released**, unchanged from every previously reported result.
    *Disclosure of timing:* this rule extension was written after seeing the 0.0013 internal ordering. The evidence
    that decides it (the external validation) was pre-planned in #16 and already computed; the rule is
    direction-neutral and would have excluded B2 on the same evidence. `RELEASED.json` records the excluded arms in
    `excluded_by_external_check`. No metric, model or reported number changed as a result — only which arm the file
    names.
22. **Arm PB4 (multi-label MLP) runs on the Apple-silicon GPU (MPS) (2026-09-18, before PB4 was run).** The CPU is the
    shared bottleneck, so the pre-registered network is trained on MPS; model, loss, optimiser, hyper-parameter grid and
    mini-batch order are unchanged. Dropout masks come from the device's own random generator, so an MPS run is
    equivalent to a different random seed rather than a bitwise copy of a CPU run (test problem: mean |Δp| 0.016,
    AUPRC Δ ≈ 0.005). The device is recorded in every OOF row's `config`; `PB4_DEVICE=cpu` reproduces on CPU. The
    scikit-learn and XGBoost models cannot use MPS (their GPU back-ends are CUDA-only) and stay on the CPU.

## Deviations from the protocol (declared)
| protocol item | here | reason |
|---|---|---|
| ~~Bacformer-large / ESM-C 300M~~ | ~~Bacformer-base / ESM-2 35M~~ — **superseded**: both encoder families were run (ESM-C 300M + Bacformer-large computed on Kaggle GPUs, arms C_large/E_large/E_large_MAG/B2E_large/B2XLE_large) | initially 34 h locally; later run on 2x T4 |
| G3 GTDB taxonomy | NCBI primary, GTDB sensitivity (S2) | splits frozen before GTDB was obtained |
| full KofamScan; KO+Pfam variant | 852-KO energy panel; no Pfam | local compute |
| B depth {4,6,8} | logistic C grid (4) + XGB depth {2,4} | N = 358; shallower trees |
| C top-k pooling | mean and max only | budget |
| E checkpoint sweep | complete-genomes (primary) + MAG (E_MAG, reported separately) | added after audit |
| F fine-tuning | not run | stopping rule: C did not beat B |
| G retrieval → rules | partly: curated HMM retrieval (MagicLamp LithoGenie/FeGenie) used as features (LG), no rule layer | compute/curation |
| H LLM over KOs (MCP) | run as baseline 2 on 96 genomes: Claude and GPT-6 Astra, evidence from the TogoMCP BH26 RDF endpoint, anonymised vs named | docs/baseline2_llm_togomcp.md |
| phylo_D, constraint_violation_rate | not reported | no tree; no hierarchical label constraints defined |
| decay curve vs phylogenetic distance | AUPRC by lowest taxonomic rank shared with any training genome | no tree |
| pooled OOF AUPRC | pooled over outer folds (also mean per-fold AUPRC in `metrics_per_fold.csv`) | pooled is the primary metric |

## AutoML on the pathway task (2026-09-21, decision of the authors)
Addendum O pre-registered AutoML on the donor/acceptor task and said the pathway task would follow "only if da
shows a gain". It did not (AM 0.5688 vs RF 0.6124, CI excluding zero). We ran the pathway task anyway,
so it is reported as an additional confirmatory run, with the same fixed settings (FLAML 2.7.0, estimator space lgbm/xgboost/rf/extra_tree/lrl2, metric average
precision, 30 s per fit, frozen pathway folds, and the prevalence-filtered 7 594-column matrix that RF/ET use --
not PB2's 10 517-column matrix). Prediction recorded before the run: AM does not beat PB2 (0.5000 genus) or ET (0.5057) by >= 0.02 macro
AUPRC.

## scikit-learn downgrade during the TabPFN experiment (2026-09-21, caught and fixed)
Installing tabpfn 2.2.1 pulled scikit-learn down from 1.9.1 to 1.6.1. The released models are pickled under 1.9.1,
so `predict.py` emitted InconsistentVersionWarning while predicting for the 178 validation genomes. That run was
killed, scikit-learn 1.9.1 was reinstalled, and the predictions were recomputed from scratch; no result in the
repository comes from a model unpickled under the wrong version. Two runs that happened while 1.6.1 was installed
are unaffected because they fit fresh models rather than loading pickles: the TabPFN arm (TPFN, which uses its own
transformer) and, for part of its run, nothing else — the pathway AutoML process had started before the downgrade
and therefore held 1.9.1.

**Outcome (2026-09-21).** AM on the pathway task: macro AUPRC **0.4472 genus / 0.3930 family**, against ET 0.5057,
PB2 0.5000, RF 0.4978. Paired bootstrap: AM − ET −0.060 [−0.078, −0.042] and AM − RF −0.052 [−0.070, −0.033] at
genus, both CIs excluding zero; the same at family. The pre-recorded prediction (no gain of >= 0.02 over PB2 or ET)
is confirmed, and more strongly than for donors/acceptors: automated search is **significantly worse** than the
hand-tuned forests on both tasks. With 3 378 genomes and 7 594 features, 30 s of search per fit is spent
overfitting the searcher's internal splits.


## Deployment model for TogoMCP set to B2XLE_large (2026-09-21, decision of the authors)
The pre-registered release rule selects B2 (results/models/RELEASED.json, `preregistered_rule_choice`). We
chose to deploy the deep-learning fusion B2XLE_large instead, after the prospective head-to-head on 100 novel
organisms placed it first among the five saved models by AUROC (0.746) and AUPRC (0.858), with B2 second (0.732,
0.851) and the difference not statistically resolved (−0.014 [−0.116, +0.084]). This is a post-hoc choice on
prospective evidence and is recorded as such; the paper reports both models. B2 stays the fallback for genomes that
have no Bacformer-large vector. Both final models are trained on all 358 labelled genomes with per-trait class
balancing.


## Deployment model for the pathway task set to PE (2026-09-22, decision of the authors)
The pre-registered rule releases PB3. We chose the deep-learning fusion PE for TogoMCP after the prospective
head-to-head on 99 novel organisms placed it first by AUROC (0.701 vs 0.673) with a confidence interval spanning
zero, while PB3 stays ahead on AUPRC and on the released yes/no calls, and on internal cross-validation
(0.501 vs 0.486). This is a preference for the sequence-informed model under statistical equivalence, not evidence
of superiority; the paper reports both. PB3 remains the fallback for genomes without a Bacformer-large vector.
