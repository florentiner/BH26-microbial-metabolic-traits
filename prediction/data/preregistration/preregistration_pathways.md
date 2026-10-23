# Pre-registration — energy metabolism pathway prediction from proteomes

Written 2026-09-18, **before any pathway model was fitted or scored**. Companion to `docs/preregistration.md`
(electron donors/acceptors). The SHA-256 of this file is recorded in the run manifest of the first scoring run.

## Motivation (project brief)
Metagenome-assembled genomes are accumulating faster than their physiology can be characterised. The aim is to
predict, from a predicted proteome alone, which energy-metabolism pathways an organism uses — to help choose
cultivation conditions for uncultured lineages. Two consequences for the design:
* evaluation must be **out-of-clade** (novel lineages are the use case), and
* evaluation must be **robust to incomplete genomes** (MAGs are partial), which is tested explicitly here.

## Task
Multi-label classification. Input: a genome's predicted proteome → KofamScan KO annotation (from the RDF Portal
BH26 store, the resource exposed through TogoMCP; verified equivalent to a local KofamScan run to 1 cell in 43 344).
Output: one probability per pathway trait.

## Labels
Source: Madin et al. 2020 trait synthesis (`pathways` field), literature-curated, independent of the KO features.
Pathway terms are merged into traits (fixed now):

| trait | Madin terms merged |
|---|---|
| denitrification | denitrification, nitrous_oxide_reduction, nitric_oxide_reduction |
| dnra | nitrate_reduction_to_ammonia, nitrite_reduction_to_ammonia |
| nitrate_reduction | nitrate_reduction, nitrite_reduction (+ the two above) |
| nitrogen_fixation | nitrogen_fixation |
| nitrification | nitrification, ammonia_oxidation |
| methanogenesis | methanogenesis and all its subtypes (H2CO2, acetoclastic, methanol, methylamine, …) |
| acetogenesis | acetogenesis |
| sulfate_reduction | sulfate_reduction |
| sulfur_reduction | sulfur_reduction, thiosulfate_reduction, sulfite_reduction |
| sulfur_oxidation | sulfur_oxidation, sulfur_oxidation_dark, thiosulfate_oxidation*, sulfide_oxidation*, sulfite_oxidation*, anoxygenic_photoautotrophy_sulfur_oxidation |
| iron_reduction | iron_reduction |
| iron_oxidation | iron_oxidation, iron_oxidation_dark |
| manganese_reduction | manganese_reduction |
| arsenate_reduction | arsenate_reduction |
| fumarate_reduction | fumarate_reduction |
| hydrogen_oxidation | hydrogen_oxidation_dark, anoxygenic_photoautotrophy_hydrogen_oxidation |
| fermentation | fermentation |
| aerobic_chemoheterotrophy | aerobic_chemo_heterotrophy, aerobic_heterotrophy |
| photoautotrophy | photoautotrophy, anoxygenic_photoautotrophy, oxygenic_photoautotrophy |
| photoheterotrophy | photoheterotrophy, aerobic_anoxygenic_phototrophy |
| methylotrophy | methylotrophy, methanol_oxidation, methylamine_oxidation |
| methanotrophy | methane_oxidation |

A trait is evaluable with ≥ 25 positives and ≥ 25 negatives. Labels are **closed-world** (absence of a term is
treated as a negative), as in the donor/acceptor task; the FAPROTAX-sourced share of positives is reported per trait,
and a no-FAPROTAX sensitivity analysis is run.

## Genomes and features
One RefSeq genome per species (reference genome if available, else the highest assembly level, then most recent),
restricted to species whose genome has a KO profile in the BH26 store. Features, all genome-derived:
* **P1 (primary)**: significant KofamScan KO presence/absence, full profile.
* **P2**: the 852-KO energy-metabolism panel only (comparable to the donor/acceptor arm B2).
* **P3**: KEGG module completeness vector (71 energy modules) — the rule-based representation.
* **P4**: P1 + P3.
No taxonomy, organism name or literature is given to any model.

## Splits and evaluation
Frozen before training and hashed: 5 outer folds grouped by **genus** and, separately, by **family** (GTDB where
available, NCBI otherwise), assigned by the same balanced greedy allocation as the donor/acceptor task. Inner grouped
CV selects hyper-parameters and thresholds; test folds are never used for selection.
Metrics per trait × arm × level: AUPRC (primary), AUPRC lift, AUROC, MCC and F1 at the inner-CV threshold, Brier;
95 % cluster-bootstrap CIs over split groups; BH q-values per trait and max-T simultaneous intervals across arms.

## Baselines (must be beaten to claim anything)
1. Prevalence null.
2. Taxonomic group-mean null (lowest shared rank with the training fold).
3. **KEGG module completeness rule**: the score is the completeness of the module(s) mapped to that trait — the
   current state of the art in practice (DRAM/METABOLIC-style annotation), and the baseline this project's own
   pathway readout provides.
4. LLM over the same KO evidence (arm H), on a subset, anonymised.

## Arms
| arm | features | model |
|---|---|---|
| PA | module completeness rule | none (threshold on completeness) |
| PB1 | P2 energy panel | logistic regression (C grid 4) / XGBoost (depth 2, 4), 3 seeds |
| PB2 | P1 full KO profile | same |
| PB3 | P4 = KO profile + module completeness | same |
| PE | P1 + Bacformer-large genome embedding (if GPU budget allows) | KO model averaged with a logistic probe |

## MAG-robustness experiment (core, not optional)
Because the application is metagenome-assembled genomes, every released model is additionally evaluated on
**simulated incomplete genomes**: proteins are randomly removed to leave 90 %, 70 % and 50 % of the proteome
(3 draws each, fixed seeds), features recomputed, and the same out-of-fold models applied. Reported as AUPRC versus
completeness per trait. A model that collapses below 90 % completeness is not useful for MAGs regardless of its
score on complete genomes.

## Pre-registered predictions
* PP1. The full KO profile (PB2) beats the 852-KO panel (PB1) — more pathways are in scope than energy metabolism alone.
* PP2. Learned models beat the KEGG module completeness rule by ≥ 0.15 macro AUPRC at genus level.
* PP3. Traits with a single canonical module (methanogenesis, sulfate reduction, nitrogen fixation) are predicted best
  (AUPRC > 0.8); broad physiological categories (fermentation, aerobic chemoheterotrophy) worst, because they are
  defined by absence of alternatives rather than by marker genes.
* PP4. Performance degrades gracefully to 70 % completeness (< 0.05 macro AUPRC loss) and sharply at 50 %.
* PP5. Pathway prediction is *easier* than the donor/acceptor categories (macro AUPRC higher), because the labels
  name a mechanism rather than an organism-level category.

## Addendum — foundation-model arms for pathways (written at 22:33 system time on 2026-09-18, before any pathway embedding was computed)
**Vectors:** ESM-C 300M per protein → Bacformer-large (masked-complete-genomes @ 84f85c37), mean over proteins →
960 numbers per genome, computed on the user's Kaggle GPUs exactly as for the donor/acceptor genomes
(kaggle/esmc_bacformer_large); proteins re-downloaded from NCBI in genome order; ESM-C vectors of the 851 genomes
already embedded are reused. No local data is uploaded for this step.
**Arms:** **PE** as pre-registered (PB2's KO model averaged ½ + ½ with an L2 logistic probe on the Bacformer-large
vector, configuration from the PB2 grid in inner CV); **PE_BFL**: the probe alone (the foundation model without
gene content); **PE_MLP**: one multi-label MLP (20 outputs) on the Bacformer-large vector, the frozen-vector neural
counterpart. Same frozen folds, nested CV, nulls and metrics as PB2. Fine-tuning Bacformer-large is decided only
after these results (it gave no gain for donors/acceptors, docs/preregistration_trees_xai.md Addendum N).
**Predictions.** PE-P1: PE beats PB2 on held-out genera (macro AUPRC, paired CI excluding zero), as B2XL +
Bacformer-large did for donors/acceptors. PE-P2: PE_BFL alone does not beat PB2. PE-P3: the gain of PE over PB2 is
smaller at family than at genus level (the vector partly encodes lineage).

## Addendum — fine-tuning Bacformer-large for pathways (written at 23:46 system time on 2026-09-18, before the kernel was built or run)
**Choice of model (from the variants already evaluated):** Bacformer-large, the best genome foundation model on the
donor/acceptor task (0.577 alone; 0.625 fused with gene content). Fine-tuning there gave no gain over the frozen
model at 358 genomes; the pathway task has 3 378 labelled genomes, so it is retried here.
**Arms:** **FTP_BFL2** — Bacformer-large layers 1–28 frozen and cached; layers 29–30, the final LayerNorm and a head
(mean over proteins → 960 → 256 → 20 outputs, all pathways jointly) trained on the labels; **FTP_BFL0** — the same
head on the frozen 30-layer mean vector (ablation). Training fixed in advance: AdamW, lr 2e-5 (transformer) and 1e-3
(head), weight decay 0.01, BCE with per-pathway pos_weight = n_neg/n_pos capped at 50, one genome per step,
gradient accumulation 16, fp16 autocast, 8 epochs, epoch chosen by macro AUPRC on the next frozen fold of the same
split level (inner validation), per-pathway MCC cut-offs from that fold. Frozen pathway genus and family folds.
**Fusion (computed locally from the out-of-fold predictions, no extra training):** **FTP_FUSE** = ½ PB2 + ½ FTP_BFL2.
**Data on Kaggle:** at the user's request (fine-tuning needs the labels next to the GPU), the 3 378 pathway labels
(derived from the public Madin et al. 2020 trait synthesis) and fold ids are embedded in a private kernel under the
user's account; ESM-C vectors are read from the user's earlier embedding kernels.
**Predictions.** FT-P1: FTP_BFL2 beats FTP_BFL0 by ≥ 0.02 macro AUPRC on held-out genera. FT-P2: FTP_BFL2 alone
does not beat PB2. FT-P3: FTP_FUSE beats PB2 on held-out genera (paired CI excluding zero).

## Outcomes of both foundation-model addenda (written at 13:20 system time on 2026-09-20, after the results were computed)
Vectors: all 3 378 pathway genomes embedded (Kaggle kernels redox-pw-embed-1..4, 0 proteins missing).
Fine-tuning: yahoodan/redox-pw-finetune, all 20 cross-validation jobs completed (10.5 h on 2 T4 GPUs).
Macro AUPRC, held-out genera / families: PB2 0.500 / 0.459, ET 0.506 / 0.433, PE 0.486 / 0.426,
FTP_FUSE 0.483 / 0.425, PE_MLP 0.417 / 0.325, PE_BFL 0.401 / 0.322, FTP_BFL2 0.356 / 0.289,
FTP_BFL0 0.345 / 0.288, taxonomic null 0.340 / 0.281.
Paired bootstrap differences (results/trees_xai/pw/paired.csv), genus and family:
* PE − PB2 −0.013 [−0.023, −0.002] and −0.031 [−0.047, −0.016] → **PE-P1 refuted**: the fusion is worse than the
  KO model, not better, at both levels (the opposite of the donor/acceptor result, where B2XLE_large gained).
* PE_BFL − PB2 −0.098 [−0.120, −0.077] and −0.134 [−0.163, −0.108] → **PE-P2 confirmed** (the probe alone is worse).
* PE-P3 held in direction but not as a gain: the deficit of PE is larger at family level (−0.031 vs −0.013), which
  is consistent with the vector carrying lineage rather than mechanism.
* FTP_BFL2 − FTP_BFL0 +0.012 [+0.002, +0.022] genus → **FT-P1 refuted** (the pre-registered bar was ≥ 0.02;
  training two transformer layers buys almost nothing over the frozen vector). Family +0.003 [−0.005, +0.011], tied.
* FTP_BFL2 − PB2 −0.143 [−0.168, −0.118] → **FT-P2 confirmed**.
* FTP_FUSE − PB2 −0.017 [−0.028, −0.006] genus, −0.032 [−0.048, −0.013] family → **FT-P3 refuted**.
Conclusion: with 3 378 labelled genomes — nine times the donor/acceptor set — a genome foundation model, frozen or
fine-tuned, does not reach gene content for pathway traits, and fusing it with the KO model costs accuracy. The
released pathway model therefore stays the KO-based one. Code: src/62_pathway_fm.py, src/61_finetune_collect.py pw.

## Addendum — LLM baseline for the pathway task (written at 01:40 system time on 2026-09-21, before any answer was collected)
**Why:** the donor/acceptor task has an LLM baseline (docs/baseline2_llm_togomcp.md); the user asked for the same
comparison on pathways. Claude sub-agents only — no paid API, so no GPT arm.
**Protocol, identical to src/21_baseline2_llm.py:** 96 genomes drawn with seed 20260921, rarest trait first so every
trait has at least 5 positives in the subset (achieved: minimum 5, maximum 48); evidence = significant KofamScan
hits from the RDF Portal BH26 Microbes endpoint (TogoMCP), restricted to the 852-KO energy panel, given as
KO | definition | number of proteins | best score/threshold; batches of 12 genomes; probabilities for all 20 traits.
**Arms:** **P_H_genome** (anonymised codes, evidence only) and **P_H_named** (same evidence plus the organism name).
**Comparison:** macro AUPRC on the same 96 genomes against PB2 (released pathway model), PB1 and PE, with a paired
genus-clustered bootstrap.
**Predictions.** L-P1: P_H_genome is worse than PB2 (CI excluding zero), as Claude-genome was for donors/acceptors
(−0.088 vs B2). L-P2: P_H_named beats P_H_genome, because the name allows literature recall instead of genome
reading. L-P3: the named arm's advantage is larger here than for donors/acceptors, because pathway labels derive
from FAPROTAX, which is itself taxonomy-driven and therefore recallable from the organism name.

## Outcome of the pathway LLM baseline (revised 2026-09-21 after audit; supersedes the first version)
All 1 920 answers returned for both arms (96 genomes x 20 traits). Macro AUPRC **on the 96-genome subset**, which
is enriched for rare traits and therefore reads far above the full-set numbers (PB2 0.659 here vs 0.500 on all
3 378 genomes; the taxonomic null 0.531 here vs 0.340). Comparisons are only valid within this column:

| Arm | subset macro AUPRC | full-set macro AUPRC |
|---|---|---|
| ET | 0.681 | 0.506 |
| **PB3 (the released pathway model)** | **0.675** | 0.501 |
| RF | 0.667 | 0.498 |
| **P_H_named** (LLM, organism name given) | **0.663** | — |
| PB2 | 0.659 | 0.500 |
| PE | 0.638 | 0.486 |
| PB1 (852-KO panel, the evidence the LLM saw) | 0.609 | 0.467 |
| **null_taxNN (taxonomic null)** | **0.531** | 0.340 |
| **P_H_genome** (LLM, evidence only) | **0.471** | — |
| PA (KEGG module rule) | 0.317 | 0.213 |

Paired genus-clustered bootstrap (results/pathway_llm_contrasts.csv):
* P_H_genome − PB3 = **−0.203 [−0.282, −0.123]**; − PB1 (matched evidence) = −0.138 [−0.231, −0.068].
  **L-P1 confirmed**: reading gene evidence, the LLM is far worse than the KO model.
* P_H_genome − null_taxNN = **−0.059 [−0.150, +0.043]**. **The evidence-only LLM does not beat the taxonomic
  null on this task** — the point estimate is below it. This is the pre-registered baseline 2, it was missing from
  the first version of this outcome section, and it *inverts* the donor/acceptor result, where the same arm cleared
  the null (+0.093 [+0.011, +0.150]). Any statement that the pathway result simply mirrors donors/acceptors is
  wrong: there, the LLM read genomes better than taxonomy; here it does not.
* P_H_named − P_H_genome = **+0.191 [+0.131, +0.252]**: **L-P2 confirmed**, the organism name is worth more than
  the gene evidence. Against the released arm, P_H_named − PB3 = −0.012 [−0.077, +0.059], a tie.
* **L-P3 underpowered rather than refuted**: +0.191 here against +0.199 for donors/acceptors, with CIs of about
  ±0.06 each and no test of the difference of differences.
**Corrections recorded by the audit.** (i) The released pathway model is **PB3**, not PB2 (results/models/pathways/
RELEASED.json); the first version of this section compared against PB2 and thereby reported a parity claim that
reverses against the released arm. (ii) The claim that the energy panel "contains no photosystem KOs" is false —
the panel holds 55 photosystem/light-harvesting KOs; what is true is that it holds no bacteriochlorophyll or puf
KOs and no arsenic KOs, and that every photo-positive genome in this subset is anoxygenic. (iii) The subset's
rare-trait enrichment (5-7x) was not quantified for the pathway task as it was for donors/acceptors
(results/qc/llm_subset_seed_sensitivity.json); the full-set column above is given instead.

## Provenance limitation of these addenda (added 2026-09-21 after audit)
These pathway addenda were written in the working tree and hashed into a local .sha256 chain, but the repository
had no commit between 2026-09-18 and 2026-09-21, so **the claim that each addendum predates its experiment cannot
be verified from the artefacts**, and the hash chain never covered the one committed version. The predictions and
their outcomes are reported exactly as recorded, but they should be read as *stated in advance and recorded
locally*, not as externally timestamped pre-registration. From 2026-09-21 the directory is committed, and any
future addendum is committed before the run it describes.

## Training-budget asymmetry in the fine-tuned arms (added 2026-09-21 after audit)
The Kaggle fine-tuning jobs (FTP_BFL2, FTP_BFL0) hold out one fold for testing **and a second fold for epoch
selection**, so they train on 3/5 of the genomes, while PB1, PB2, PE_BFL and PE_MLP train on 4/5 under nested CV.
The fine-tuned arms are therefore handicapped by 25% less training data in FT-P1, FT-P2, FT-P3 and in FTP_FUSE.
The size of the effect is visible between the two frozen-vector models: PE_BFL (4/5 of the data) 0.401 versus
FTP_BFL0 (3/5) 0.345. The conclusion "a genome foundation model does not reach gene content" is still supported by
the frozen arms, which are budget-matched (PE_BFL 0.401 and PE_MLP 0.417 against PB2 0.500), but it is **not**
established at full strength by the fine-tuned arms alone.

## Addendum — external prospective test for the pathway task (written 2026-09-21, before any pathway prediction
## exists for these genomes)
**Why:** the pathway task has no external test at all. Its labels come entirely from the Madin compilation, and
eight of its twenty traits have 100% FAPROTAX-derived positives, so a test built from primary literature is the
single most valuable thing that can be added. The 167 organisms of the donor/acceptor prospective validation
(selection_v2.md) already carry blind, quoted growth tests from their description papers; several of those tests
map directly onto pathway traits. No pathway model has been run on any of these genomes at the time of writing.

**Mapping, fixed now.** A published donor/acceptor growth test is carried over to a pathway trait only where the
two denote the same process:

| Evidence category (task 1) | Pathway trait (task 2) | Exact? |
|---|---|---|
| acceptor_sulfate | sulfate_reduction | yes |
| acceptor_sulfur | sulfur_reduction | yes |
| donor_sulfur_oxidation | sulfur_oxidation | yes |
| acceptor_nitrate | nitrate_reduction | yes (the trait covers denitrification and DNRA together) |
| acceptor_iron | iron_reduction | yes |
| donor_iron_oxidation | iron_oxidation | yes |
| donor_hydrogen_oxidation | hydrogen_oxidation | yes |
| acceptor_oxygen | aerobic_chemoheterotrophy | **approximate** — an aerobic chemolithoautotroph satisfies the acceptor test without being a chemoheterotroph |
| acceptor_CO2_methanogen_acetogen | methanogenesis or acetogenesis | **approximate** — the evidence does not separate the two; scored against max(P_methanogenesis, P_acetogenesis) |

Expected yield from the current evidence: 166 organisms, 439 calls, of which 337 fall under the seven exact
mappings. Traits with no mapping are not tested here.

**Model:** the released pathway model **PB3** (full KO profile + KEGG module completeness), trained on all 3 378
genomes, unchanged. Binary calls at its stored cut-offs.
**Scoring:** identical to the donor/acceptor prospective test. Primary analysis = blind evidence only and exact
mappings only; the two approximate mappings are reported separately; "deployed mode" uses PB3's own 11 reliable
traits. Baselines on the same calls: a constant rule at training prevalence, and an oracle per-category majority.
Confidence intervals are cluster bootstraps over description papers.
**Predictions.** X-P1: PB3 beats the training-prevalence constant rule on MCC by a clear margin (> 0.15).
X-P2: the per-trait ordering broadly follows internal cross-validation — sulfate reduction and sulfur oxidation
near the top, iron reduction near the bottom. X-P3: iron reduction is close to blind (sensitivity < 0.25), as it
was for the donor/acceptor model, because KEGG has no orthology for dissimilatory metal reduction.

## Outcome — external prospective test for the pathway task (written 2026-09-22 after scoring)
Balanced 99 organisms (one per genus, labels-only selection over the seven exact mappings), 203 scored calls from
99 description papers; PB3 trained on all 3 378 genomes, unchanged. Results (validation_new_organisms/
pathway_summary.csv, pathway_per_trait.csv, novel100_pathway_comparison.csv):
* PB3 accuracy 59.6% [51.1, 67.0], **MCC 0.218**, sensitivity 46.3% (50/108), specificity 74.7% (71/95);
  training-prevalence constant rule 34.0% (MCC −0.323); oracle per-trait majority 72.4%.
* **X-P1 confirmed**: PB3 beats the deployable constant rule by 0.54 MCC (> 0.15).
* **X-P2 partly confirmed**: sulfate reduction is the best trait (accuracy 90%, MCC 0.76) and iron reduction is weak
  (MCC 0.17), as predicted; but sulfur oxidation (MCC 0.00, 12 calls) and iron oxidation (0 of 5 positives found)
  are worse than internal CV suggested, and hydrogen oxidation finds only 4 of 17 positives.
* **X-P3 refuted in the favourable direction**: iron-reduction sensitivity is 0.476 (10/21), not < 0.25 — the
  pathway model, unlike the donor/acceptor model, finds about half of the published iron reducers.
* Head-to-head of every saved pathway model on the same 203 calls (AUROC / AUPRC averaged over traits):
  **PE (PB2 + Bacformer-large) 0.701 / 0.806**, **PB3 0.673 / 0.812**, PB2 0.650 / 0.807, PB1 0.628 / 0.781.
  PE − PB3 = +0.028 AUROC [−0.056, +0.128]: a tie, PE ahead on ranking, PB3 ahead on AUPRC and on the released
  yes/no calls (MCC 0.193 vs 0.179). The 100-organism scale cannot separate them.
The approximate mappings (aerobic chemoheterotrophy; methanogenesis-or-acetogenesis), reported apart: 50 calls,
accuracy 64.0%, MCC 0.299.
