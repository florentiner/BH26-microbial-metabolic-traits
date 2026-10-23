# Addendum A to `preregistration_additional_experiments.md` — METABOLIC and DRAM as X1 comparators

Written 2026-10-14 ~23:10 MSK, before METABOLIC or DRAM was run on any genome. Same evaluation set, metrics and
scoring as X1 (100 newly described species, 255 blind growth tests; accuracy, MCC, sensitivity, specificity;
paired cluster bootstrap over papers for B2 − tool). Both tools give presence/absence only, so no AUROC.

## METABOLIC-G v4 (github.com/AnantharamanLab/METABOLIC @9723633)
Run on the RefSeq protein FASTA of each genome (`-in`), default thresholds and motif checks, KOfam profiles from
the same KOfam release used by our models (2026-08-25; only the ~2 800 profiles METABOLIC references are pressed).
A category is predicted positive when METABOLIC reports **any** of the listed functions as present:

| Our category | METABOLIC Category \| Function |
|---|---|
| donor_sulfur_oxidation | Sulfur cycling \| Sulfide oxidation; Sulfur oxidation; Thiosulfate oxidation |
| donor_ammonia_oxidation | Nitrogen cycling \| Ammonia oxidation |
| donor_nitrite_oxidation | Nitrogen cycling \| Nitrite oxidation |
| donor_hydrogen_oxidation | Hydrogenases \| Ni-Fe Hydrogenase |
| donor_CO_oxidation | C1 metabolism \| Aerobic CO oxidation |
| donor_iron_oxidation | Iron cycling \| Iron oxidation |
| acceptor_oxygen | Oxygen metabolism \| cytochrome c oxidase caa3- or cbb3-type, or quinol oxidase aa3- or bo-type (bd-type excluded: common in anaerobes) |
| acceptor_nitrate | Nitrogen cycling \| Nitrate reduction |
| acceptor_sulfate | Sulfur cycling \| Sulfate reduction |
| acceptor_sulfur | Sulfur cycling \| Sulfur reduction; Sulfite reduction |
| acceptor_CO2_methanogen_acetogen | Methane metabolism \| Methane production; Carbon fixation \| Wood Ljungdahl pathway |
| acceptor_iron | Iron cycling \| Iron reduction; Metal (Iron/Manganese) reduction |

## DRAM (KOfam mode, distill step)
DRAM's database build needs ~128 GB RAM, which is not available locally. Its trait calls come from its `distill`
step applied to annotations; we supply KOfam annotations made with the same KOfam release and the per-profile
thresholds that DRAM's KOfam mode uses, restricted to the KOs referenced by DRAM's distill forms, and run DRAM's own
distill code. The category mapping onto DRAM's distillate (function names in DRAM's `genome_summary_form` /
`product` output) will be fixed in Addendum B after reading DRAM's forms and before running distill. If DRAM's
distill cannot be run faithfully, it is reported as not run.

## Predictions
* A-P1: B2's MCC exceeds METABOLIC's on the 255 calls (direction only; no threshold).
* A-P2: METABOLIC's sensitivity exceeds its specificity less than B2's does (rule-based tools over-call
  capabilities from single marker genes).
