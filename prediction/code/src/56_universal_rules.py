"""Universal rules: gene rules that hold across phyla, built only from genes with a known function
(docs/preregistration_trees_xai.md, Addendum J).

Candidates: the trait's curated pathway genes + aerobic/anaerobic lifestyle markers. A rule is eligible only if its
hits span >= 3 phyla (>= 3 genomes each); sequential covering picks the rule with the highest worst-phylum precision.
Validation: leave one phylum out - learn on all other phyla, score on the unseen phylum.

usage: python src/56_universal_rules.py --task {da,ox,pw}
Output: results/rulebook/universal_<task>.csv, universal_<task>_lopo.csv
"""
import argparse
import importlib
from itertools import combinations

import numpy as np
import pandas as pd

from gene_mechanisms import AEROBIC_MARKERS, ANAEROBIC_MARKERS, PATHWAY, TRAIT_LIFESTYLE
from pathway_common import ko_matrix

tx = importlib.import_module("47_trees_xai")
rb = importlib.import_module("51_rulebook")
ri = importlib.import_module("54_rule_interpretation")
OUT = tx.RES / "rulebook"
MIN_PHYLA, MIN_PER_PHYLUM, MIN_LB = 3, 3, 0.60
N_RULES = 3
LOPO_MIN_POS = 5
RANK_BY_COVERAGE = False   # set True by --posthoc (Addendum K b)


# enzyme / complex families: subunits of one enzyme are interchangeable for rule identity (Addendum K, post hoc)
ENZYME_FAMILY = {
    "K00370": "nar", "K00371": "nar", "K00373": "nar", "K00374": "nar", "K02567": "nap", "K02568": "nap",
    "K02274": "aa3", "K02275": "aa3", "K02276": "aa3", "K00404": "cbb3", "K00405": "cbb3", "K00406": "cbb3",
    "K00425": "bd", "K00426": "bd", "K02297": "bo3", "K02298": "bo3", "K02257": "heme-o/a", "K02259": "heme-o/a",
    "K00410": "bc1", "K00413": "bc1", "K11180": "dsr", "K11181": "dsr", "K27187": "dsrMK", "K27188": "dsrMK",
    "K27190": "dsrMK", "K00394": "apr", "K00395": "apr", "K16950": "asr", "K16951": "asr", "K16952": "asr",
    "K00399": "mcr", "K00401": "mcr", "K00402": "mcr", "K03388": "hdr", "K03389": "hdr", "K03390": "hdr",
    "K00174": "kor", "K00176": "kor", "K00177": "kor", "K00192": "cdh", "K00193": "cdh", "K00194": "cdh",
    "K00197": "cdh", "K02588": "nif", "K02586": "nif", "K02591": "nif", "K17222": "sox", "K17223": "sox",
    "K17224": "sox", "K17225": "sox", "K17226": "sox", "K17227": "sox", "K10944": "amo", "K10945": "amo",
    "K10946": "amo", "K04561": "nor", "K02305": "nor", "K08928": "puf", "K08929": "puf", "K02689": "psa",
    "K02690": "psa", "K02703": "psb", "K02706": "psb", "K06281": "hya/hyb", "K05922": "hya/hyb", "K06282": "hya/hyb",
    "K00436": "hox", "K18007": "hox", "K14028": "mdh", "K14029": "mdh", "K23995": "mdh", "K00016": "ldh",
    "K03778": "ldh", "K03518": "cox-CO", "K03519": "cox-CO", "K03520": "cox-CO",
}


def signature(rule):
    return tuple(sorted({(ENZYME_FAMILY.get(g, g), v) for g, v in rule}))


def lineage_unit(lab, y):
    """Phylum if the trait occurs in >= 3 phyla (>= 3 genomes each), else class, else order (Addendum K)."""
    for col in ("phylum", "class", "order"):
        if col not in lab.columns:
            continue
        u = lab[col].fillna("?").to_numpy()
        n = pd.Series(u[y == 1]).value_counts()
        if (n >= MIN_PER_PHYLUM).sum() >= 3:
            return col, u
    col = "order" if "order" in lab.columns else "phylum"
    return col, lab[col].fillna("?").to_numpy()


def features(task, lab):
    if task == "da":
        return tx.load_task("da")[3].reset_index(drop=True)
    X = ko_matrix(lab.genome_id)            # full profile: rare markers must not be filtered away
    return X.reset_index(drop=True)


MAX_GENES = 20


def candidates(trait, cols):
    """All pathway genes, then the lifestyle markers matching the trait's oxygen lifestyle (both kinds when the
    trait is facultative or either), capped at MAX_GENES."""
    need = TRAIT_LIFESTYLE.get(trait, "either")
    life = (ANAEROBIC_MARKERS if need == "anaerobic" else AEROBIC_MARKERS if need.startswith("aerobic") or
            need == "mostly aerobic" else ANAEROBIC_MARKERS | AEROBIC_MARKERS)
    genes = [k for k in sorted(PATHWAY.get(trait, set())) if k in cols]
    genes += [k for k in sorted(life) if k in cols and k not in genes]
    return genes[:MAX_GENES]


def enumerate_rules(genes):
    conds = [(g, v) for g in genes for v in (1, 0)]
    return [tuple(sorted(c)) for L in (1, 2, 3) for c in combinations(conds, L) if len({g for g, _ in c}) == L]


def rule_masks(X, rules):
    B = {g: X[g].to_numpy() > 0 for g in {g for r in rules for g, _ in r}}
    M = np.ones((len(rules), len(X)), dtype=bool)
    for i, r in enumerate(rules):
        for g, v in r:
            M[i] &= B[g] if v else ~B[g]
    return M


def learn(M, y, phy_onehot, cls, idx, min_units=MIN_PHYLA):
    """Sequential covering on the genomes in idx; returns chosen rule indices and their per-phylum stats."""
    target = (y[idx] == cls).astype(np.float32)
    Msub = M[:, idx].astype(np.float32)
    P = phy_onehot[idx]
    hits_p = Msub @ P                                   # rules x phyla: genomes the rule fires on
    corr_p = Msub @ (P * target[:, None])               # ... of which correct
    n, k = hits_p.sum(1), corr_p.sum(1)
    base = target.mean()
    ok_ph = hits_p >= MIN_PER_PHYLUM
    n_ph = ok_ph.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        prec_p = np.where(ok_ph, corr_p / np.maximum(hits_p, 1), np.nan)
        worst = np.nanmin(np.where(ok_ph, prec_p, np.inf), axis=1)
    lb = np.array([rb.wilson_lb(int(a), int(b)) for a, b in zip(k, n)])
    eligible = (n_ph >= min_units) & (lb >= MIN_LB) & (lb > base)
    uncovered = target.astype(bool).copy()
    chosen = []
    for _ in range(N_RULES):
        new = (Msub[:, uncovered].sum(1) > 0) if uncovered.any() else np.zeros(len(M), bool)
        cand = np.where(eligible & new)[0]
        if not len(cand):
            break
        good = cand[worst[cand] >= 0.8]
        if RANK_BY_COVERAGE and len(good):                             # Addendum K (b): broadest reliable rule
            best = good[np.argmax(Msub[good][:, target.astype(bool)].sum(1))]
        else:
            best = cand[np.lexsort((-n[cand], -worst[cand]))[0]]      # pre-registered: worst-unit precision first
        chosen.append(best)
        uncovered &= ~Msub[best].astype(bool)
    stats = {i: dict(n_phyla=int(n_ph[i]), worst_phylum_precision=float(worst[i]), overall_precision=float(k[i] / n[i]),
                     support=int(n[i]), lb=float(lb[i])) for i in chosen}
    return chosen, stats, prec_p, ok_ph


def posthoc(task, lab, traits, Y, X, ext):
    """Addendum K: lineage unit adapts to the trait; enzyme-level rule identity; pooled unseen-unit precision."""
    rows = []
    for j, t in enumerate(traits):
        y = Y[:, j].astype(int)
        genes = candidates(t, set(X.columns))
        if not genes:
            continue
        level, units = lineage_unit(lab, y)
        uu = np.unique(units)
        onehot = (units[:, None] == uu[None, :]).astype(np.float32)
        n_trait_units = int((pd.Series(units[y == 1]).value_counts() >= MIN_PER_PHYLUM).sum())
        min_units = max(2, min(3, n_trait_units))
        rules = enumerate_rules(genes)
        M = rule_masks(X, rules)
        test_units = [u for u in uu if ((units == u) & (y == 1)).sum() >= LOPO_MIN_POS]
        for cls, label in ((1, "has trait"), (0, "lacks trait")):
            learned = {}
            for u in test_units:
                tr, te = np.where(units != u)[0], np.where(units == u)[0]
                remaining = int((pd.Series(units[tr][y[tr] == 1]).value_counts() >= MIN_PER_PHYLUM).sum())
                chosen, _, _, _ = learn(M, y, onehot, cls, tr, max(2, min(3, remaining)))
                learned[u] = (te, chosen)
            chosen, stats, prec_u, ok_u = learn(M, y, onehot, cls, np.arange(len(y)), min_units)
            for rank, i in enumerate(chosen, 1):
                r, sig = rules[i], signature(rules[i])
                # audit C3: rules of one enzyme signature learned in the same run fire on overlapping genomes, so
                # their hits are merged (union) per held-out unit - each unseen genome counts once
                kh = nh = runs = n_unseen_good = 0
                fired_units = []
                for u, (te, lst) in learned.items():
                    same = [ii for ii in lst if signature(rules[ii]) == sig]
                    if same:
                        runs += 1
                        hit = M[same][:, te].any(axis=0)
                        c_u, h_u = int((hit & (y[te] == cls)).sum()), int(hit.sum())
                        kh += c_u
                        nh += h_u
                        n_unseen_good += int(h_u >= MIN_PER_PHYLUM and c_u / h_u >= 0.8)
                        if h_u:
                            fired_units.append(f"{u} {c_u}/{h_u}")
                # audit round 2 N5: large lineages (>= 50 in-sample hits) in which the rule was never tested unseen
                tested = {f.rsplit(" ", 1)[0] for f in fired_units}
                big = [(uu[q], int((M[i] & (units == uu[q])).sum())) for q in range(len(uu))]
                untested_big = [f"{u_} ({n_} hits)" for u_, n_ in sorted(big, key=lambda x: -x[1])
                                if n_ >= 50 and u_ not in tested]
                lo, hi = rb.wilson(kh, nh) if nh else (np.nan, np.nan)
                per_unit = sorted(((uu[q], prec_u[i, q], int((M[i] & (units == uu[q])).sum()))
                                   for q in np.where(ok_u[i])[0]), key=lambda x: -x[2])
                n_good = sum(1 for _, pr, _ in per_unit if pr >= 0.8)
                stab = runs / len(learned) if learned else np.nan
                gv, _, _, gwhy = ri.grounding(list(r), label, t)
                base = float((y == cls).mean())
                row = dict(task=task, trait=t, cls=label, rank=rank, lineage_unit=level, n_trait_units=n_trait_units,
                           rule=rb.rule_text(r, cls, t), genes=", ".join(g for g, _ in r),
                           enzyme_signature=" + ".join(f"{f} {'present' if v else 'absent'}" for f, v in sig),
                           base_rate=base, overall_precision=stats[i]["overall_precision"], support=stats[i]["support"],
                           coverage=float((M[i] & (y == cls)).sum() / max((y == cls).sum(), 1)),
                           n_units=stats[i]["n_phyla"], n_units_precision_ge_80=n_good,
                           worst_unit_precision=stats[i]["worst_phylum_precision"],
                           per_unit="; ".join(f"{p} {pr:.0%} (n={n})" for p, pr, n in per_unit[:8]),
                           lopo_runs=len(learned), lopo_stability_enzyme=stab, unseen_hits=nh, unseen_correct=kh,
                           unseen_units_tested=runs, unseen_units_precision_ge_80=n_unseen_good,
                           unseen_units_fired=len(fired_units), unseen_by_unit="; ".join(fired_units),
                           large_units_never_tested="; ".join(untested_big),
                           unseen_precision=kh / nh if nh else np.nan, unseen_ci_low=lo, unseen_ci_high=hi,
                           grounding=gv, grounding_reasons=gwhy,
                           conditions=" | ".join(ri.condition_text(g, bool(v), t) for g, v in r),
                           logic=ri.logic_text([(g, bool(v)) for g, v in r], label, t))
                row["universal"] = bool(nh >= MIN_PER_PHYLUM and gv != "not grounded" and n_good >= min_units
                                        and kh / nh >= 0.8 and lo >= 0.6 and lo > base and stab >= 0.5)
                # audit round 3: a universal claim needs unseen evidence in >= 2 lineages where the rule fired and
                # >= 20 distinct unseen genomes
                row["enough_unseen_evidence"] = bool(len(fired_units) >= 2 and nh >= 20)
                row["universal"] = row["universal"] and row["enough_unseen_evidence"]
                if ext is not None:
                    ex, el = ext
                    col = rb.EXTERNAL_LABEL.get(t, t)
                    if col in el.columns and all(g in ex.columns for g, _ in r):
                        le = el[col]
                        okm = le.notna().to_numpy()
                        me = np.ones(len(ex), bool)
                        for g, v in r:
                            me &= (ex[g].to_numpy() > 0) == bool(v)
                        me &= okm
                        ke, ne = int((me & (le.fillna(-1).to_numpy() == cls)).sum()), int(me.sum())
                        elo, _ = rb.wilson(ke, ne) if ne else (np.nan, np.nan)
                        row.update(external_precision=ke / ne if ne else np.nan, external_support=ne,
                                   external_ci_low=elo, external_base_rate=float((le[okm] == cls).mean()))
                        # audit M9a: a rule the independent external test contradicts is not universal; the test
                        # needs >= 5 external genomes of the predicted class to be informative (round 2 n1)
                        n_cls_ext = int((le[okm] == cls).sum())
                        row["external_n_class"] = n_cls_ext
                        if n_cls_ext < 5 or (le[okm] != cls).sum() < 5:
                            row["external_uninformative"] = True
                        elif ne >= MIN_PER_PHYLUM and not (elo > row["external_base_rate"] and elo >= 0.6):   # round 3
                            row["universal"] = False
                            row["failed_external"] = True
                rows.append(row)
        print(f"[{task}/post hoc] {t}: unit={level} ({n_trait_units} units with the trait), "
              f"{sum(1 for r in rows if r['trait'] == t and r['universal'] and r['cls'] == 'has trait')} universal "
              f"positive rules", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / f"universal_{task}_posthoc.csv", index=False, float_format="%.4f")
    u = d[d.universal & (d.cls == "has trait")]
    print(f"\n[{task}/post hoc] universal rules predicting a trait: {len(u)}")
    print(u[["trait", "enzyme_signature", "lineage_unit", "n_units", "worst_unit_precision", "unseen_precision",
             "unseen_hits", "lopo_stability_enzyme", "coverage"]].round(2).to_string(index=False))


# the core marker whose presence should imply each trait (complete enzyme = all listed genes present)
CORE_MARKER = {
    "methanogenesis": (["K00399", "K00401", "K00402"], "methyl-CoM reductase mcrABG"),
    "nitrogen_fixation": (["K02588", "K02586", "K02591"], "nitrogenase nifHDK"),
    "sulfate_reduction": (["K11180", "K11181", "K00394"], "dsrAB + aprA"),
    "acceptor_sulfate": (["K11180", "K11181", "K00394"], "dsrAB + aprA"),
    "denitrification": (["K04561", "K00376"], "norB + nosZ"),
    "photoautotrophy": (["K02703", "K02706"], "photosystem II psbA + psbD"),
    "photoheterotrophy": (["K08928", "K08929"], "reaction centre pufL + pufM"),
    "donor_ammonia_oxidation": (["K10944", "K10945", "K10946"], "amoABC"),
    "methylotrophy": (["K14028", "K14029"], "methanol dehydrogenase mxaFI"),
    "acetogenesis": (["K14138", "K00194", "K00197"], "Wood-Ljungdahl acsB + cdhDE"),
    "uses_O2": (["K02274", "K02275", "K02276"], "aa3 oxidase coxABC"),
    "acceptor_oxygen": (["K02274", "K02275", "K02276"], "aa3 oxidase coxABC"),
}


def label_completeness(task, lab, traits, Y, X):
    """For each trait's complete core enzyme: carriers, share labelled with the trait, share of trait carriers."""
    rows = []
    for t, (kos, name) in CORE_MARKER.items():
        if t not in traits or not all(k in X.columns for k in kos):
            continue
        y = Y[:, traits.index(t)]
        m = np.all([X[k].to_numpy() > 0 for k in kos], axis=0)
        rows.append(dict(task=task, trait=t, core_marker=name, carriers=int(m.sum()), carriers_labelled=int(y[m].sum()),
                         share_of_carriers_labelled=float(y[m].mean()) if m.any() else np.nan,
                         share_of_trait_with_marker=float(y[m].sum() / max(y.sum(), 1)),
                         n_phyla_with_marker=int(pd.Series(lab.phylum.fillna("?").to_numpy()[m]).nunique())))
    pd.DataFrame(rows).to_csv(OUT / f"label_completeness_{task}.csv", index=False, float_format="%.4f")
    if rows:
        print(pd.DataFrame(rows).round(2).to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=tx.TASKS, required=True)
    ap.add_argument("--posthoc", action="store_true", help="Addendum K revision (labelled post hoc)")
    a = ap.parse_args()
    task = a.task
    lab, traits, Y, _ = tx.load_task(task)
    lab = lab.reset_index(drop=True)
    X = features(task, lab)
    rb.register_task(X, traits, Y)
    phyla = lab.phylum.fillna("?").to_numpy()
    uph = np.unique(phyla)
    onehot = (phyla[:, None] == uph[None, :]).astype(np.float32)
    ext = None
    if task == "da":
        ex, el = rb.load_external()
        ext = (ex, el)
    if a.posthoc:
        global RANK_BY_COVERAGE
        RANK_BY_COVERAGE = True
        posthoc(task, lab, traits, Y, X, ext)
        label_completeness(task, lab, traits, Y, X)
        return
    rows, lopo_rows, proc_rows = [], [], []
    all_idx = np.arange(len(lab))
    for j, t in enumerate(traits):
        y = Y[:, j].astype(int)
        genes = candidates(t, set(X.columns))
        if not genes:
            continue
        rules = enumerate_rules(genes)
        M = rule_masks(X, rules)
        test_phyla = [p for p in uph if ((phyla == p) & (y == 1)).sum() >= LOPO_MIN_POS]
        for cls, label in ((1, "has trait"), (0, "lacks trait")):
            # ---- leave one phylum out
            lopo_sets, pooled = {}, {}

            for p in test_phyla:
                tr = np.where(phyla != p)[0]
                te = np.where(phyla == p)[0]
                chosen, _, _, _ = learn(M, y, onehot, cls, tr)
                lopo_sets[p] = set(chosen)
                for i in chosen:
                    m = M[i, te]
                    a = pooled.setdefault(i, [0, 0, 0])
                    a[0] += int((m & (y[te] == cls)).sum())
                    a[1] += int(m.sum())
                    a[2] += int(m.sum() > 0)
                    lopo_rows.append(dict(task=task, trait=t, cls=label, held_out_phylum=p,
                                          rule=rb.rule_text(rules[i], cls, t), hits=int(m.sum()),
                                          correct=int((m & (y[te] == cls)).sum())))
                union = M[chosen][:, te].any(axis=0) if chosen else np.zeros(len(te), bool)   # audit C3
                proc_rows.append(dict(trait=t, cls=label, held_out_phylum=p, hits=int(union.sum()),
                                      correct=int((union & (y[te] == cls)).sum())))
            # ---- final rules on all phyla
            chosen, stats, prec_p, ok_ph = learn(M, y, onehot, cls, all_idx)
            for rank, i in enumerate(chosen, 1):
                r = rules[i]
                kh, nh, nph = pooled.get(i, [0, 0, 0])
                lo, hi = rb.wilson(kh, nh) if nh else (np.nan, np.nan)
                per_phylum = sorted(((uph[q], prec_p[i, q], int((M[i] & (phyla == uph[q])).sum()))
                                     for q in np.where(ok_ph[i])[0]), key=lambda x: -x[2])
                n_good = sum(1 for _, pr, _ in per_phylum if pr >= 0.8)
                stab = sum(i in s for s in lopo_sets.values()) / len(lopo_sets) if lopo_sets else np.nan
                gv, _, _, gwhy = ri.grounding(list(r), label, t)
                row = dict(task=task, trait=t, cls=label, rank=rank, rule=rb.rule_text(r, cls, t),
                           genes=", ".join(g for g, _ in r), base_rate=float((y == cls).mean()),
                           overall_precision=stats[i]["overall_precision"], support=stats[i]["support"],
                           coverage=float((M[i] & (y == cls)).sum() / max((y == cls).sum(), 1)),
                           n_phyla=stats[i]["n_phyla"], n_phyla_precision_ge_80=n_good,
                           worst_phylum_precision=stats[i]["worst_phylum_precision"],
                           per_phylum="; ".join(f"{p} {pr:.0%} (n={n})" for p, pr, n in per_phylum[:8]),
                           lopo_runs=len(lopo_sets), lopo_stability=stab, lopo_precision=kh / nh if nh else np.nan,
                           lopo_ci_low=lo, lopo_ci_high=hi, lopo_phyla_fired=nph, grounding=gv, grounding_reasons=gwhy,
                           conditions=" | ".join(ri.condition_text(g, bool(v), t) for g, v in r),
                           logic=ri.logic_text([(g, bool(v)) for g, v in r], label, t))
                if ext is not None:
                    ex, el = ext
                    col = rb.EXTERNAL_LABEL.get(t, t)
                    if col in el.columns and all(g in ex.columns for g, _ in r):
                        le = el[col]
                        okm = le.notna().to_numpy()
                        me = np.ones(len(ex), bool)
                        for g, v in r:
                            me &= (ex[g].to_numpy() > 0) == bool(v)
                        me &= okm
                        ke, ne = int((me & (le.fillna(-1).to_numpy() == cls)).sum()), int(me.sum())
                        elo, ehi = rb.wilson(ke, ne) if ne else (np.nan, np.nan)
                        row.update(external_precision=ke / ne if ne else np.nan, external_support=ne,
                                   external_ci_low=elo, external_ci_high=ehi,
                                   external_base_rate=float((le[okm] == cls).mean()))
                # criterion (ii) of Addendum J needs precision in phyla the rule never saw: a rule that was never
                # tested in an unseen phylum cannot be called universal (no fallback to in-sample precision)
                row["tested_in_unseen_phyla"] = bool(nh >= MIN_PER_PHYLUM)
                row["universal"] = bool(row["tested_in_unseen_phyla"] and gv != "not grounded" and n_good >= 3
                                        and row["lopo_precision"] >= 0.8 and lo >= 0.6 and stab >= 0.5
                                        and lo > row["base_rate"])
                rows.append(row)
        n_u = sum(1 for r in rows if r["trait"] == t and r["universal"])
        print(f"[{task}] {t}: {len(genes)} candidate genes, {len(test_phyla)} test phyla, {n_u} universal rules",
              flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / f"universal_{task}.csv", index=False, float_format="%.4f")
    lp = pd.DataFrame(lopo_rows)
    lp.to_csv(OUT / f"universal_{task}_lopo.csv", index=False)
    if proc_rows:   # procedure level: rule list learned without a phylum, scored on it (each genome once)
        pr = pd.DataFrame(proc_rows).groupby(["trait", "cls"]).agg(
            hits=("hits", "sum"), correct=("correct", "sum"), phyla=("held_out_phylum", "nunique")).reset_index()
        pr["precision_in_unseen_phyla"] = pr.correct / pr.hits.clip(lower=1)
        pr.to_csv(OUT / f"universal_{task}_procedure.csv", index=False, float_format="%.4f")
    u = d[d.universal]
    print(f"\n[{task}] universal rules: {len(u)} ({int((u.cls == 'has trait').sum())} predicting a trait)")
    print(u[u.cls == "has trait"][["trait", "rule", "n_phyla", "worst_phylum_precision", "lopo_precision",
                                   "lopo_stability", "coverage"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
