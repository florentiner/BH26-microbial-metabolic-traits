# Predicting microbial energy metabolism from genomes — BH26 preprint

Paper, released models, data and code for the BioHackrXiv preprint from DBCLS BioHackathon 2026
(Kawashima, Ezhov, Kinjo, Koblitz, Nakazato, Okabbepu, Otsuka). This folder holds the prediction part
(donor and acceptor classifiers, pathway classifiers, their evaluation); the RDF resource, TogoMCP integration,
literature collection and viewer are described in the paper and maintained by the co-authors.

## Contents

| Folder | What it holds |
|---|---|
| `paper/` | `BH26_preprint.pdf`, a Markdown version (`BH26_preprint.md`), the Overleaf project (`overleaf/`, compile with LuaLaTeX) and the pandoc source (`source/`) |
| `models/` | The five released models with metadata and `SHA256SUMS` (Supplementary Table S19) |
| `data/` | Training labels, features and frozen folds, the external genomes, the newly described species with quoted evidence, the blind tests, MAG lists, comparator-tool calls, LLM inputs and outputs, and all preregistration documents |
| `results/tables/` | The metric tables behind the tables and figures of the paper |
| `code/` | Every script of the project (`src/`), helpers (`scripts/`), scoring of new species (`validation_scoring/`) and the build of the paper (`paper_build/`), with pinned package versions (`requirements.lock`) |

## Models

| Name in the paper | Folder | Predicts | Input |
|---|---|---|---|
| RedoxKO (final model) | `models/RedoxKO` | 14 donor and acceptor categories, calls in 6, confidence labels | 852-KO panel, random forest, CPU |
| RedoxKO-LX | `models/RedoxKO-LX` | same, first final model, for comparison | 852-KO panel, logistic regression or XGBoost |
| RedoxKO-ens-BF | `models/RedoxKO-ens-BF` | same | KO panel and Bacformer-large genome vector (GPU) |
| EnergyKO-M | `models/EnergyKO-M` | 20 pathway traits, calls in 11 | full KO profile and KEGG module completeness |
| EnergyKO-BF | `models/EnergyKO-BF` | 20 pathway traits | EnergyKO and Bacformer-large vector |

Internal code names used in the scripts and result tables: RedoxKO = `RF`, RedoxKO-LX = `B2`, RedoxKO-ens =
`B2XL`, RedoxKO-ens-BF = `B2XLE_large`, EnergyKO-M = `PB3`, EnergyKO-BF = `PE`, extra-trees pathway model = `ET`.
Each `model.joblib` is a dict from category to {model, threshold, calibrator}. Check files with
`cd models && shasum -a 256 -c SHA256SUMS`.

## Predicting a genome

```bash
python -m venv venv && venv/bin/pip install -r code/requirements.lock
sh code/setup_layout.sh          # links models/ and data/ into code/ in the layout the scripts expect
# put the KOfam profiles/ and ko_list.gz (release 25 August 2026, https://www.genome.jp/ftp/db/kofam/)
# into code/data/raw/kofam/, then:
venv/bin/python code/src/predict_with_labels.py GCF_000006985.1 --out out.csv
venv/bin/python code/src/predict_with_labels.py --fasta my_genome.fna --out out.csv --pathways none
```
The default model is RedoxKO; `--arm B2` or `--arm B2XLE_large` selects another. `predict_with_labels.py` runs
`predict.py` and adds the confidence labels. It takes about 15 s per genome from an accession on a CPU.

## Data

| Folder | Paper section |
|---|---|
| `data/donor_acceptor/` | Methods, donor and acceptor classifiers: 358 labeled genomes x 14 categories, 852-KO panel, KO matrix, Bacformer-large vectors, GTDB taxonomy, frozen genus/family folds with SHA-256 and the other fold allocations (Table S5) |
| `data/external_genomes/` | 562 genomes of unseen genera labeled from the larger dataset (Madin et al. 2020) |
| `data/pathways/` | 3,378 genomes x 20 traits with folds, the 71 KEGG energy modules, module completeness, full KO matrix |
| `data/new_species/` | Tests of newly described species, each with its source and quoted string; the primary 100, pathway 99 and autotroph-like 99 sets; protocol fixed before prediction |
| `data/blind_tests/` | The two blind tests of the frozen model (Table S20): extracted tests, audit logs, frozen labels with hashes, predictions and scores |
| `data/mags/` | MAG-isolate pairs and the 300 habitat MAGs (GTDB release 232) |
| `data/comparator_tools/` | Calls of METABOLIC, DRAM, FAPROTAX, Microbe Decoder and the BacDive aerobic model |
| `data/llm_comparison/` | Prompts (input batches), subset keys and raw outputs of the LLM comparison |
| `data/preregistration/` | Every preregistration document with its SHA-256, and the log of deviations (Supplementary Note 1) |

## Not included

Genome FASTA and protein files (downloaded from NCBI by `code/src/01_download_genomes.py`), the KOfam profiles
(8.7 GB), the weights of ESM-2, ESM-C and Bacformer (fetched from Hugging Face on first use) and the phenotype
collection itself, which belongs to the co-authors' resource.
