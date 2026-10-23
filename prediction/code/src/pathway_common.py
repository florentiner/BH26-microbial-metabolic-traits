"""Shared loading and feature construction for the pathway-prediction task (docs/preregistration_pathways.md)."""
import importlib
import os
from pathlib import Path

import numpy as np
import pandas as pd

from redox_common import PROC, sha256

ROOT = Path(__file__).resolve().parents[1]
PATH = Path(os.environ.get("PATHWAY_DATA", ROOT / "data/pathways"))  # overridable for smoke tests
RES = ROOT / "results"
# taxonomic ranks for the group-mean null, coarse -> fine (Madin columns; superkingdom renamed to domain)
TAX = ["domain", "phylum", "class", "order", "family", "genus"]


def load_pathways():
    """Labelled genomes with frozen folds, trait list and label matrix Y."""
    lab = pd.read_csv(PATH / "labels.tsv", sep="\t")
    expected = (PATH / "labels.sha256").read_text().split()[0]
    assert sha256(PATH / "labels.tsv") == expected, "pathway labels.tsv changed after freezing"
    meta = {"species_tax_id", "species", "genus", "family", "order", "class", "phylum", "superkingdom", "genome_id",
            "gtdb_genus", "gtdb_family", "faprotax_source", "data_source", "pathways", "group_genus", "group_family",
            "fold_genus", "fold_family"}
    traits = [c for c in lab.columns if c not in meta]
    lab = lab.rename(columns={"superkingdom": "domain"})
    return lab, traits, lab[traits].to_numpy(dtype=float)


def ko_matrix(genome_ids):
    m = pd.read_csv(PATH / "ko_matrix.tsv", sep="\t", index_col=0)
    return m.reindex(genome_ids).fillna(0).astype(np.float32)


def module_completeness(genome_ids):
    """KEGG energy-module completeness for the pathway genomes (computed from their KO profiles)."""
    f = PATH / "module_completeness.tsv"
    if not f.exists():
        pm = importlib.import_module("29_pathway_modules")
        mods = pd.read_csv(PROC / "kegg_modules_energy.tsv", sep="\t")
        ko = pd.read_csv(PATH / "ko_matrix.tsv", sep="\t", index_col=0)
        cols = np.array(ko.columns)
        out = pd.DataFrame(index=ko.index, columns=mods.module, dtype=float)
        steps_by_mod = {m.module: pm.parse_steps(m.definition) for m in mods.itertuples()}
        for acc, row in zip(ko.index, ko.to_numpy()):
            kos = set(cols[row > 0])
            for mid, steps in steps_by_mod.items():
                out.at[acc, mid] = np.mean([pm.step_present(st, kos) for st in steps]) if steps else np.nan
        out.index.name = "genome_id"
        out.to_csv(f, sep="\t", float_format="%.4f")
    return pd.read_csv(f, sep="\t", index_col=0).reindex(genome_ids).fillna(0).astype(np.float32)


def energy_panel_subset(genome_ids):
    """The 852-KO energy-metabolism panel, as used by the donor/acceptor arm B2."""
    panel = pd.read_csv(PROC / "ko_panel.tsv", sep="\t").KO
    m = ko_matrix(genome_ids)
    keep = [k for k in panel if k in m.columns]
    return m[keep]


# pathway trait -> KEGG modules whose completeness is the rule-based baseline (arm PA), fixed before scoring
TRAIT_MODULES = {
    "denitrification": ["M00529"],
    "dnra": ["M00530"],
    "nitrate_reduction": ["M00529", "M00530"],
    "nitrogen_fixation": ["M00175"],
    "nitrification": ["M00528", "M00804"],
    "methanogenesis": ["M00567", "M00357", "M00563", "M00356", "M00358"],
    "acetogenesis": ["M00377"],
    "sulfate_reduction": ["M00596"],
    "sulfur_reduction": ["M00596", "M00986"],
    "sulfur_oxidation": ["M00595"],
    "iron_reduction": [],
    "iron_oxidation": [],
    "manganese_reduction": [],
    "arsenate_reduction": [],
    "fumarate_reduction": [],
    "hydrogen_oxidation": [],
    "fermentation": ["M00579", "M00357"],
    "aerobic_chemoheterotrophy": ["M00149", "M00151", "M00155"],
    "photoautotrophy": ["M00161", "M00163", "M00165"],
    "photoheterotrophy": ["M00161", "M00597"],
    "methylotrophy": ["M00174", "M00346", "M00345"],
    "methanotrophy": ["M00174"],
}


def rule_scores(lab, traits):
    """Arm PA: per trait, the maximum completeness over its mapped modules (0 where no module maps)."""
    mc = module_completeness(lab.genome_id)
    S = np.zeros((len(lab), len(traits)), dtype=np.float32)
    for j, t in enumerate(traits):
        mods = [m for m in TRAIT_MODULES.get(t, []) if m in mc.columns]
        if mods:
            S[:, j] = mc[mods].to_numpy().max(axis=1)
    return S


def completeness_from_matrix(X, columns, modules=None):
    """KEGG module completeness from an arbitrary KO presence matrix (used for degraded/MAG-like proteomes)."""
    pm = importlib.import_module("29_pathway_modules")
    mods = modules if modules is not None else pd.read_csv(PROC / "kegg_modules_energy.tsv", sep="\t")
    cols = np.asarray(columns)
    steps_by_mod = [(m.module, pm.parse_steps(m.definition)) for m in mods.itertuples()]
    out = np.zeros((X.shape[0], len(steps_by_mod)), dtype=np.float32)
    for i in range(X.shape[0]):
        kos = set(cols[X[i] > 0])
        for j, (_, steps) in enumerate(steps_by_mod):
            out[i, j] = np.mean([pm.step_present(st, kos) for st in steps]) if steps else 0.0
    return out
