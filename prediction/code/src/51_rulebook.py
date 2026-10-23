"""Rulebook: high-probability IF-THEN rules built from the XAI results (docs/preregistration_trees_xai.md, Addendum G).

usage: python src/51_rulebook.py --task {da,ox,pw}
Outputs (results/rulebook/): <task>_rules.csv (final rules + all validation), <task>_nested.csv (per fold),
                             <task>_procedure.csv (held-out performance of the whole rule list per trait)
"""
import argparse
import importlib
import re
from itertools import combinations

import numpy as np
import pandas as pd

from gene_mechanisms import ABSENCE_MEANING, LIFESTYLE, MECH, PATHWAY, condition_support, register_task

tx = importlib.import_module("47_trees_xai")
OUT = tx.RES / "rulebook"
MAX_GENES, MAX_LEN = 16, 3
MIN_GENERA, MIN_GENOMES, MIN_LB = 5, 8, 0.60
N_POS_RULES, N_NEG_RULES = 3, 2
HIGH_PREC, HIGH_LB = 0.80, 0.60
EXTERNAL_LABEL = {"donor_ammonia_oxidation": "nitrification"}   # external set names this trait more broadly


def load_external():
    """External genomes (data/external), restricted to lineage-novel genomes: deviation #13 removes the 20 whose
    genus is in training (11 of them are among the genomes the rules were learned on)."""
    root = tx.RES.parent / "data/external"
    ex = pd.read_csv(root / "ko_matrix_panel.tsv", sep="\t", index_col=0)
    novel = pd.read_csv(root / "lineage.tsv", sep="\t").set_index("genome_id").novel
    ex = ex[novel.reindex(ex.index).fillna(False).astype(bool).to_numpy()]
    el = pd.read_csv(root / "labels.tsv", sep="\t").set_index("genome_id").reindex(ex.index)
    return ex, el


def wilson_lb(k, n, z=1.96):
    if n == 0:
        return 0.0
    p = k / n
    return (p + z * z / (2 * n) - z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / (1 + z * z / n)


def wilson(k, n):
    lo = wilson_lb(k, n)
    hi = 1 - wilson_lb(n - k, n) if n else 1.0
    return lo, hi


def candidate_genes(task, trait, cols):
    d = tx.OUT / task
    order = []

    def add(kos):
        for k in kos:
            if k in cols and k not in order:
                order.append(k)
    if (d / "consensus.csv").exists():
        c = pd.read_csv(d / "consensus.csv")
        c = c[c.trait == trait]
        add(c[c.status == "robust"].sort_values("mean_rank").KO)
    if (d / "xai_stable_rf.csv").exists():
        s = pd.read_csv(d / "xai_stable_rf.csv")
        add(s[s.trait == trait].sort_values("mean_rank").KO)
    if (d / "surrogate_rules.csv").exists():
        r = pd.read_csv(d / "surrogate_rules.csv")
        add(re.findall(r"K\d{5}", " ".join(r[r.trait == trait].rule.astype(str))))
    cmh_f = d / "alt_cmh_size_adjusted.csv"
    known = tx.known_markers(task, trait, tx.reference_groups())
    if cmh_f.exists():
        cm = pd.read_csv(cmh_f)
        cm = cm[(cm.trait == trait) & (cm.q_bh < 0.05) & cm.KO.isin(known)].sort_values("p")
        add(cm.KO)
    return order[:MAX_GENES]


def candidate_genes_infold(Xtr, ytr, genera_tr, min_genera=MIN_GENERA):
    """Strict nesting (audit M1): candidate genes chosen from the training folds only. Genomes are collapsed to genera
    (mean presence, mean label) so one large genus cannot dominate; genes are ranked by |correlation| between genus
    presence and genus label across genera, among genes present in at least `min_genera` genera."""
    g = pd.Series(genera_tr)
    P = pd.DataFrame(Xtr.to_numpy() > 0, columns=Xtr.columns).groupby(g.to_numpy()).mean()
    yl = pd.Series(ytr).groupby(g.to_numpy()).mean().reindex(P.index).to_numpy()
    keep = (P > 0).sum(axis=0).to_numpy() >= min_genera
    A = P.to_numpy()[:, keep]
    if yl.std() == 0 or not keep.any():
        return []
    a = (A - A.mean(0)) / np.where(A.std(0) > 0, A.std(0), np.inf)
    r = (a * ((yl - yl.mean()) / yl.std())[:, None]).mean(0)
    order = np.argsort(-np.abs(r))[:MAX_GENES]
    return list(np.asarray(P.columns)[keep][order])


def all_rules(genes):
    conds = [(g, v) for g in genes for v in (1, 0)]
    out = []
    for L in range(1, MAX_LEN + 1):
        for combo in combinations(conds, L):
            if len({g for g, _ in combo}) == L:
                out.append(tuple(sorted(combo)))
    return out


def masks(Xdf, rules):
    B = {g: Xdf[g].to_numpy() > 0 for g in {g for r in rules for g, _ in r}}
    M = np.ones((len(rules), len(Xdf)), dtype=bool)
    for i, r in enumerate(rules):
        for g, v in r:
            M[i] &= B[g] if v == 1 else ~B[g]
    return M


def learn(Xdf, y, genera, rules, M, cls):
    """Sequential covering for one class (1 = trait present, 0 = absent)."""
    target = (y == cls)
    uncovered = target.copy()
    chosen = []
    n_max = N_POS_RULES if cls == 1 else N_NEG_RULES
    gcodes = pd.factorize(genera)[0]
    for _ in range(n_max):
        best, best_lb = None, MIN_LB
        for i in range(len(rules)):
            m = M[i]
            n = int(m.sum())
            if n < MIN_GENOMES or not (m & uncovered).any():
                continue
            if len(np.unique(gcodes[m])) < MIN_GENERA:
                continue
            lb = wilson_lb(int((m & target).sum()), n)
            if lb > best_lb:
                best, best_lb = i, lb
        if best is None:
            break
        chosen.append(best)
        uncovered &= ~M[best]
    return chosen


def rule_text(rule, cls, trait):
    parts = [f"{MECH.get(g, (g, ''))[0]} ({g}) {'present' if v else 'absent'}" for g, v in rule]
    return "IF " + " AND ".join(parts) + f" THEN {'has' if cls == 1 else 'lacks'} {trait}"


KO_NAME = {}


def ko_name(g):
    if not KO_NAME:
        import gzip
        with gzip.open(tx.RES.parent / "data/raw/kofam/ko_list.gz", "rt") as f:
            d = pd.read_csv(f, sep="\t", usecols=["knum", "definition"])
        cut = lambda x: x if len(x) <= 60 else x[:60].rsplit(" ", 1)[0].rstrip(",;[(") + " …"   # noqa: E731
        KO_NAME.update(zip(d.knum, d.definition.map(cut)))
    return KO_NAME.get(g, g)


def interpret(rule, cls, trait):
    """Why the rule should work, from gene function (Addendum H levels; direction-aware since audit C5):
    mechanistic = every informative condition is a pathway gene in the fitting state or an absence that rules the trait
    out; physiological = those plus oxygen-lifestyle markers that fit the trait; contains unexplained conditions = at least one
    condition with no known link, or in the state that contradicts the biology."""
    cl = "has trait" if cls == 1 else "lacks trait"
    notes, proxies, lifestyle, core = [], [], [], 0
    for g, v in rule:
        name, mech = MECH.get(g, (g, ""))
        if v == 1:
            notes.append(f"{name}: {mech}" if mech else f"{g}: {ko_name(g)} (no known link to {trait})")
        else:
            notes.append(f"{name} absent: {ABSENCE_MEANING.get(g, mech or ko_name(g))}")
        kind, _ = condition_support(g, bool(v), trait, cl)
        if kind in ("pathway", "absence"):
            core += 1
        elif kind == "lifestyle":
            lifestyle.append(name)
        elif kind in ("none", "opposite", "neutral"):     # neutral: required by the rule, not explained (round 3)
            proxies.append(name if name != g else g)
    if proxies or not (core or lifestyle):
        level = "contains unexplained conditions"
        verdict = (f"contains conditions that biology does not explain — no known link, uninformative, or in a state "
                   f"that contradicts the pathway ({', '.join(proxies) or 'no informative gene'}): they may reflect "
                   f"lineage rather than mechanism")
    elif lifestyle:
        level = "physiological"
        verdict = (f"pathway genes plus lifestyle markers ({', '.join(lifestyle)}): consistent physiology" if core else
                   f"oxygen-lifestyle markers only ({', '.join(lifestyle)}), no pathway gene: consistent with the "
                   f"lifestyle the trait implies")
    else:
        level = "mechanistic"
        verdict = "every informative gene belongs to the physiology the rule predicts, in the fitting state"
    return " | ".join(notes), verdict, level


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=tx.TASKS, required=True)
    task = ap.parse_args().task
    OUT.mkdir(parents=True, exist_ok=True)
    lab, traits, Y, Xdf = tx.load_task(task)
    register_task(Xdf, traits, Y)          # carrier shares for direction-aware grounding (audit round 2 N3)
    genera = lab.group_genus.to_numpy()
    fold = lab.fold_genus.to_numpy()
    mt = pd.read_csv(tx.OUT / task / "metrics_per_trait.csv")
    best = mt[mt.split_level == "genus"].groupby("arm").auprc.mean()
    best = best[~best.index.str.startswith("null")].idxmax()
    predictable = set(mt[(mt.arm == best) & (mt.split_level == "family") & mt.predictable].trait)
    ext = None
    if task == "da":
        ext = load_external()
    final_rows, nested_rows, proc_rows = [], [], []
    for j, t in enumerate(traits):
        y = Y[:, j].astype(int)
        genes = candidate_genes(task, t, set(Xdf.columns))
        if not genes:
            continue
        rules = all_rules(genes)
        M = masks(Xdf[genes], rules)
        # ---- nested validation: learn on 4 genus folds, score on the 5th
        fold_sets = {}
        for k in sorted(np.unique(fold)):
            tr, te = fold != k, fold == k
            for cls in (1, 0):
                ids = learn(Xdf[tr], y[tr], genera[tr], rules, M[:, tr], cls)
                fold_sets.setdefault(cls, {})[k] = [rules[i] for i in ids]
                for i in ids:
                    m = M[i, te]
                    nested_rows.append(dict(task=task, trait=t, fold=int(k), cls=cls, rule=rule_text(rules[i], cls, t),
                                            heldout_support=int(m.sum()),
                                            heldout_correct=int((m & (y[te] == cls)).sum())))
                if cls == 1:  # whole positive rule list as a classifier on the held-out fold (no rule: fires on none)
                    fire = M[ids][:, te].any(axis=0) if ids else np.zeros(int(te.sum()), bool)
                    yt = y[te]
                    proc_rows.append(dict(task=task, trait=t, fold=int(k), n_rules=len(ids), fired=int(fire.sum()),
                                          tp=int((fire & (yt == 1)).sum()), n_pos=int(yt.sum()), selection="main"))
        # ---- strict nesting (audit M1): the candidate genes too are chosen inside each training fold
        strict = {}
        for k in sorted(np.unique(fold)):
            tr, te = fold != k, fold == k
            g_k = candidate_genes_infold(Xdf[tr], y[tr], genera[tr])
            if not g_k:
                continue
            rules_k = all_rules(g_k)
            M_k = masks(Xdf[g_k], rules_k)
            for cls in (1, 0):
                ids = learn(Xdf[tr], y[tr], genera[tr], rules_k, M_k[:, tr], cls)
                for i in ids:
                    m = M_k[i, te]
                    key = (cls, rule_text(rules_k[i], cls, t))
                    a_, b_ = strict.get(key, (0, 0))
                    strict[key] = (a_ + int((m & (y[te] == cls)).sum()), b_ + int(m.sum()))
                if cls == 1:
                    fire = M_k[ids][:, te].any(axis=0) if ids else np.zeros(int(te.sum()), bool)
                    yt = y[te]
                    proc_rows.append(dict(task=task, trait=t, fold=int(k), n_rules=len(ids), fired=int(fire.sum()),
                                          tp=int((fire & (yt == 1)).sum()), n_pos=int(yt.sum()),
                                          selection="strict"))
        # ---- final rules on all genomes
        for cls in (1, 0):
            ids = learn(Xdf, y, genera, rules, M, cls)
            for rank, i in enumerate(ids, 1):
                r, m = rules[i], M[i]
                k_all, n_all = int((m & (y == cls)).sum()), int(m.sum())
                stab = sum(r in fold_sets.get(cls, {}).get(k, []) for k in np.unique(fold))
                nd = pd.DataFrame([x for x in nested_rows if x["trait"] == t and x["cls"] == cls
                                   and x["rule"] == rule_text(r, cls, t)])
                kh, nh = (int(nd.heldout_correct.sum()), int(nd.heldout_support.sum())) if len(nd) else (0, 0)
                mech, verdict, level = interpret(r, cls, t)
                row = dict(task=task, trait=t, predictable=t in predictable, cls="has trait" if cls else "lacks trait",
                           rank=rank, rule=rule_text(r, cls, t), n_conditions=len(r),
                           in_sample_precision=k_all / n_all, in_sample_support=n_all,
                           coverage=k_all / max(int((y == cls).sum()), 1), genera_covered=int(len(set(genera[m]))),
                           cv_stability=f"{stab}/{len(np.unique(fold))}",
                           heldout_precision=kh / nh if nh else np.nan, heldout_support=nh,
                           heldout_ci_low=wilson(kh, nh)[0] if nh else np.nan,
                           heldout_ci_high=wilson(kh, nh)[1] if nh else np.nan,
                           mechanism=mech, interpretation=verdict, interpretation_level=level,
                           base_rate=float((y == cls).mean()))
                ks, ns = strict.get((cls, rule_text(r, cls, t)), (0, 0))
                row.update(strict_heldout_precision=ks / ns if ns else np.nan, strict_heldout_support=ns,
                           strict_heldout_ci_low=wilson(ks, ns)[0] if ns else np.nan)
                row["lift_heldout"] = row["heldout_precision"] / row["base_rate"] if nh else np.nan
                row["informative_heldout"] = bool(nh and row["heldout_ci_low"] > row["base_rate"])
                if ext is not None:
                    ex, el = ext
                    col = EXTERNAL_LABEL.get(t, t)
                    if col in el.columns and all(g in ex.columns for g, _ in r):
                        lab_e = el[col].reindex(ex.index)
                        ok = lab_e.notna().to_numpy()
                        me = np.ones(len(ex), bool)
                        for g, v in r:
                            me &= (ex[g].to_numpy() > 0) == bool(v)
                        me &= ok
                        ke, ne = int((me & (lab_e.fillna(-1).to_numpy() == cls)).sum()), int(me.sum())
                        base_e = float((lab_e[ok] == cls).mean())
                        row.update(external_label=col, external_support=ne,
                                   external_n_pos=int((lab_e[ok] == 1).sum()),
                                   external_precision=ke / ne if ne else np.nan,
                                   external_ci_low=wilson(ke, ne)[0] if ne else np.nan,
                                   external_ci_high=wilson(ke, ne)[1] if ne else np.nan,
                                   external_base_rate=base_e,
                                   external_lift=(ke / ne) / base_e if ne and base_e else np.nan,
                                   informative_external=bool(ne and wilson(ke, ne)[0] > base_e))
                prec = row.get("external_precision") if not np.isnan(row.get("external_precision", np.nan)) \
                    else row["heldout_precision"]
                lb = row.get("external_ci_low") if not np.isnan(row.get("external_ci_low", np.nan)) \
                    else row["heldout_ci_low"]
                row["high_probability"] = bool(prec is not None and not np.isnan(prec) and prec >= HIGH_PREC
                                               and lb >= HIGH_LB)
                inf = row.get("informative_external") if "informative_external" in row else row["informative_heldout"]
                row["high_probability_and_informative"] = bool(row["high_probability"] and inf)
                final_rows.append(row)
        print(f"[{task}] {t}: {len(genes)} candidate genes, "
              f"{sum(1 for r in final_rows if r['trait'] == t and r['high_probability'])} high-probability rules")
    fr = pd.DataFrame(final_rows)
    fr.to_csv(OUT / f"{task}_rules.csv", index=False, float_format="%.4f")
    pd.DataFrame(nested_rows).to_csv(OUT / f"{task}_nested.csv", index=False)
    pr = pd.DataFrame(proc_rows)
    if len(pr):
        agg = pr.groupby(["trait", "selection"])[["fired", "tp", "n_pos"]].sum()
        agg["heldout_precision"] = agg.tp / agg.fired.clip(lower=1)
        agg["heldout_recall"] = agg.tp / agg.n_pos.clip(lower=1)
        agg.reset_index().to_csv(OUT / f"{task}_procedure.csv", index=False, float_format="%.4f")
    hp = fr[fr.high_probability_and_informative]
    print(f"\n[{task}] high-probability rules: {int(fr.high_probability.sum())} of {len(fr)} as pre-registered; "
          f"{len(hp)} also informative (beat the base rate)")
    if "external_precision" in fr.columns:
        pos = fr[(fr.cls == "has trait") & fr.external_precision.notna() & fr.heldout_precision.notna()]
        tr = pos.groupby("interpretation_level").agg(n=("rule", "size"), heldout=("heldout_precision", "mean"),
                                                    external=("external_precision", "mean"))
        tr["drop"] = tr.heldout - tr.external
        tr.to_csv(OUT / f"{task}_transfer_by_level.csv", float_format="%.4f")
        print("\ntransfer of positive rules to the external genomes, by interpretation level:")
        print(tr.round(3).to_string())
    print(hp[["trait", "cls", "rule", "heldout_precision", "cv_stability"] +
             (["external_precision"] if "external_precision" in hp.columns else [])].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
