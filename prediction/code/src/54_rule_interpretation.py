"""Biological interpretation of every trusted rule (high probability and informative), one rule at a time.

For each rule: what each condition means, why the combination implies the trait (pathway genes vs aerobic/anaerobic
lifestyle markers, checked against the lifestyle the trait requires), which genera it captures, where it misfires,
and a trust verdict. Genera come from the data; mechanisms from src/gene_mechanisms.py. Interpretation only - no
rule or number is changed here.

usage: python src/54_rule_interpretation.py --task {da,ox,pw}
Output: results/rulebook/<task>_trusted_interpretation.csv
"""
import argparse
import importlib
import re

import numpy as np
import pandas as pd

from gene_mechanisms import (ABSENCE_MEANING, AEROBIC_MARKERS, ANAEROBIC_MARKERS, MECH, PATHWAY, TRAIT_LIFESTYLE,
                             condition_support, register_task)

tx = importlib.import_module("47_trees_xai")
rb = importlib.import_module("51_rulebook")
OUT = tx.RES / "rulebook"


def parse(rule):
    return [(k, s == "present") for k, s in re.findall(r"\((K\d{5})\) (present|absent)", rule)]


def condition_text(ko, present, trait):
    name, mech = MECH.get(ko, (ko, ""))
    if present:
        what = mech or f"{rb.ko_name(ko)} — no known link to this trait"
        return f"{name} present: {what}"
    meaning = ABSENCE_MEANING.get(ko)
    if meaning:
        return f"{name} absent: {meaning}"
    if ko in ANAEROBIC_MARKERS:
        return f"{name} absent: no sign of {mech.split(': ')[0] if mech else 'this anaerobic enzyme'}, an enzyme typical of anaerobes"
    if ko in AEROBIC_MARKERS:
        return (f"{name} absent: no {mech.split(': ')[0] if mech else 'aerobic enzyme'} (on its own this says little "
                f"about oxygen — most organisms lack it)")
    return f"{name} absent: {mech or rb.ko_name(ko)}" + ("" if mech else " — no known link to this trait")


def lifestyle_signal(conds, path=frozenset()):
    """What the lifestyle markers in a rule say about oxygen; genes already on the trait's pathway are not counted."""
    sig = []
    for ko, present in conds:
        if ko in path:
            continue
        if ko in ANAEROBIC_MARKERS:
            sig.append(("anaerobic" if present else "no anaerobic marker", MECH.get(ko, (ko,))[0]))
        elif ko in AEROBIC_MARKERS and present:          # an absent aerobic marker is uninformative
            sig.append(("aerobic", MECH.get(ko, (ko,))[0]))
    return sig


def logic_text(conds, cls, trait):
    path = PATHWAY.get(trait, set())
    core = [MECH.get(k, (k,))[0] for k, p in conds if k in path and p]
    missing_core = [MECH.get(k, (k,))[0] for k, p in conds if k in path and not p]
    life = lifestyle_signal(conds, path)
    need = TRAIT_LIFESTYLE.get(trait, "")
    unknown = [MECH.get(k, (k,))[0] if k in MECH else k for k, p in conds
               if condition_support(k, p, trait, cls)[0] == "none"]
    parts = []
    if cls == "has trait":
        if core:
            parts.append(f"The rule requires the pathway itself ({', '.join(core)}).")
        if missing_core:
            parts.append(f"It requires {', '.join(missing_core)} to be absent, so it covers only the organisms that run "
                         f"this physiology without {', '.join(missing_core)} — one branch of it, not all.")
    else:
        if missing_core:
            parts.append(f"The rule fires when pathway gene(s) ({', '.join(missing_core)}) are missing. Whether that "
                         f"rules the trait out depends on whether another enzyme can do the same job — the measured "
                         f"precision, not the biology alone, says how often it does.")
        if core:
            parts.append(f"It also requires {', '.join(core)}, which on its own belongs to a different use of the "
                         f"pathway genes.")
    if life:
        said = "; ".join(f"{g} → {s}" for s, g in life)
        kinds = [condition_support(k, p, trait, cls)[0] for k, p in conds
                 if k not in path and k in ANAEROBIC_MARKERS | AEROBIC_MARKERS]
        consistent = bool(kinds) and all(k in ("lifestyle", "neutral") for k in kinds) and "lifestyle" in kinds
        if need in ("either", "facultative", ""):
            parts.append(f"Its oxygen-lifestyle markers ({said}) carry no information about a trait found in both "
                         f"aerobes and anaerobes; this part of the rule most likely reflects lineage.")
        elif cls == "lacks trait" and consistent:
            parts.append(f"Its lifestyle markers ({said}) point away from the {need} lifestyle the trait requires.")
        elif cls == "lacks trait":
            parts.append(f"Its lifestyle markers ({said}) do not clearly point away from the {need} lifestyle the "
                         f"trait requires — treat this part with care.")
        elif consistent:
            parts.append(f"Its lifestyle markers ({said}) agree with the {need or 'expected'} lifestyle the trait "
                         f"requires" + (", so the rule reads the pathway and the conditions it runs under together."
                                        if core else "; it contains no pathway gene, so it reads the lifestyle, not "
                                                     "the pathway itself."))
        else:
            parts.append(f"Its lifestyle markers ({said}) do not obviously fit the {need or 'expected'} lifestyle the "
                         f"trait requires — treat this part with care.")
    if unknown:
        parts.append(f"{', '.join(unknown)} has no known role here and most likely marks a lineage rather than a "
                     f"mechanism — this part of the rule may not hold outside the training lineages.")
    return " ".join(parts) or "The rule’s genes have no curated role; treat it as a statistical association."


def grounding(conds, cls, trait):
    """Is every condition connected to the trait by known biology? Computed, not judged (gene_mechanisms.
    condition_support; direction-aware since audit C5): a pathway gene counts only in the state that fits the rule, an
    informative absence only when the rule rules out the trait it argues against, an oxygen-lifestyle marker only when
    its signal fits the trait's lifestyle (and never for traits found in both lifestyles).
    Verdict: grounded (every condition supports the rule), partly grounded (at least one supports it, and at least one
    is unsupported or uninformative), not grounded (none supports it). Since audit round 3 an uninformative (neutral)
    condition also prevents a full "grounded" verdict: the rule needs it, but biology does not explain why."""
    g, u, n, why = 0, 0, 0, []
    for ko, present in conds:
        kind, text = condition_support(ko, present, trait, cls)
        if kind in ("pathway", "absence", "lifestyle"):
            g += 1
        elif kind == "neutral":
            n += 1
        else:
            u += 1
        why.append(text)
    verdict = ("grounded" if g and not (u or n) else "partly grounded" if g else "not grounded")
    return verdict, g, u, "; ".join(why)


def covered_genera(lab, X, conds, y, cls, top=5):
    m = np.ones(len(lab), bool)
    for ko, present in conds:
        m &= (X[ko].to_numpy() > 0) == present
    g = lab.genus.fillna("?").to_numpy()
    hit = m & (y == (1 if cls == "has trait" else 0))
    miss = m & (y != (1 if cls == "has trait" else 0))
    tp = pd.Series(g[hit]).value_counts().head(top)
    fp = pd.Series(g[miss]).value_counts().head(3)
    fmt = lambda s: ", ".join(f"{k} ({v})" for k, v in s.items())   # noqa: E731
    phyla = lab.phylum.fillna("?").to_numpy()
    ph = pd.Series(phyla[hit]).value_counts()
    own = pd.Series(phyla[y == (1 if cls == "has trait" else 0)]).value_counts()   # the trait's own spread
    breadth = dict(n_phyla=int(len(ph)), top_phylum=str(ph.index[0]) if len(ph) else "",
                   top_phylum_share=float(ph.iloc[0] / ph.sum()) if len(ph) else np.nan,
                   trait_n_phyla=int(len(own)),
                   trait_share_of_that_phylum=float(own.get(ph.index[0], 0) / own.sum()) if len(ph) else np.nan)
    return fmt(tp), fmt(fp), int(pd.Series(g[hit]).nunique()), int(miss.sum()), breadth


def breadth_text(br, cls):
    """Compare the phylum spread of what the rule captures with the spread of the trait itself."""
    if not br["n_phyla"]:
        return ""
    who = "organisms with the trait" if cls == "has trait" else "organisms without the trait"
    excess = br["top_phylum_share"] - br["trait_share_of_that_phylum"]
    if excess <= 0.15:
        return (f"as broad as the trait: {br['n_phyla']} of {br['trait_n_phyla']} phyla; {br['top_phylum']} makes up "
                f"{br['top_phylum_share']:.0%} of hits and {br['trait_share_of_that_phylum']:.0%} of all {who} — the "
                f"rule follows the trait’s own distribution")
    return (f"narrower than the trait: {br['top_phylum_share']:.0%} of its hits are {br['top_phylum']}, against "
            f"{br['trait_share_of_that_phylum']:.0%} of all {who} — it captures mainly that lineage’s version of the "
            f"trait ({br['n_phyla']} of {br['trait_n_phyla']} phyla)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=tx.TASKS, required=True)
    task = ap.parse_args().task
    rules = pd.read_csv(OUT / f"{task}_rules.csv")
    trusted = rules.copy()   # every rule is interpreted; the trusted subset is written separately below
    lab, traits, Y, Xdf = tx.load_task(task)
    register_task(Xdf, traits, Y)
    Xdf.index = range(len(Xdf))
    ext = None
    if task == "da":
        ex, el = rb.load_external()
        ext = (ex.reset_index(drop=True), el.reset_index())
    rows = []
    for r in trusted.itertuples():
        conds = parse(r.rule)
        j = traits.index(r.trait)
        y = Y[:, j].astype(int)
        tp, fp, n_gen, n_fp, br = covered_genera(lab, Xdf, conds, y, r.cls)
        ext_fp = ""
        if ext is not None:
            ex, el = ext
            col = rb.EXTERNAL_LABEL.get(r.trait, r.trait)
            if col in el.columns and all(k in ex.columns for k, _ in conds):
                ye = el[col].to_numpy()
                ok = ~pd.isna(ye)
                sub_lab, sub_X, sub_y = el[ok].reset_index(drop=True), ex[ok].reset_index(drop=True), ye[ok].astype(int)
                _, ext_fp, _, _, _ = covered_genera(sub_lab, sub_X, conds, sub_y, r.cls)
        gv, n_g, n_u, gwhy = grounding(conds, r.cls, r.trait)
        n_all = len(lab)
        fires = np.ones(n_all, bool)
        for ko, present in conds:
            fires &= (Xdf[ko].to_numpy() > 0) == present
        rows.append(dict(task=task, trait=r.trait, cls=r.cls, rule=r.rule, trusted=bool(r.high_probability_and_informative),
                         grounding=gv, n_grounded_conditions=n_g, n_ungrounded_conditions=n_u, grounding_reasons=gwhy,
                         fires_on_share=float(fires.mean()), fires_on_n=int(fires.sum()),
                         probability_ci_low=(r.external_ci_low if "external_ci_low" in rules.columns and
                                             pd.notna(r.external_ci_low) else r.heldout_ci_low),
                         probability_ci_high=(r.external_ci_high if "external_ci_high" in rules.columns and
                                              pd.notna(r.external_ci_high) else r.heldout_ci_high),
                         probability_source=("external genomes" if "external_precision" in rules.columns and
                                             pd.notna(r.external_precision) else "held-out genera"),
                         base_rate=(r.external_base_rate if "external_precision" in rules.columns and   # audit M7:
                                    pd.notna(r.external_precision) else r.base_rate), probability=(   # same genomes
            r.external_precision if "external_precision" in rules.columns and pd.notna(r.external_precision)
            else r.heldout_precision), coverage=r.coverage, interpretation_level=r.interpretation_level,
            external_n_pos=getattr(r, "external_n_pos", np.nan),
            conditions=" | ".join(condition_text(k, p, r.trait) for k, p in conds),
            logic=logic_text(conds, r.cls, r.trait),
            captures=tp, n_genera_captured=n_gen, misfires_training=fp, n_misfires_training=n_fp,
            misfires_external=ext_fp, **br,
            breadth=breadth_text(br, r.cls)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / f"{task}_all_interpretation.csv", index=False, float_format="%.4f")
    out[out.trusted].drop(columns="trusted").to_csv(OUT / f"{task}_trusted_interpretation.csv", index=False,
                                                    float_format="%.4f")
    print(f"[{task}] {len(out)} rules interpreted ({int(out.trusted.sum())} trusted); grounding: "
          f"{out.grounding.value_counts().to_dict()}")
    for r in out.head(3).itertuples():
        print(f"\n{r.rule}\n  logic: {r.logic}\n  captures: {r.captures}\n  misfires: {r.misfires_training} | ext: {r.misfires_external}")


if __name__ == "__main__":
    main()
