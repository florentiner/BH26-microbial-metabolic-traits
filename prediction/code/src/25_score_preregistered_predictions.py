"""Score the pre-registered predictions P1-P6 (docs/preregistration.md) against the results (audit round 3, m7).

Output: results/preregistered_predictions_scored.csv
"""
import numpy as np
import pandas as pd

from redox_common import RES

# every pre-registered embedding arm that was run, including the protocol's preferred large encoders (arm C = pooled
# ESM-C, arm E = Bacformer-large) and the attention arm D
EMB = ["C", "D", "E", "E_MAG", "C_large", "E_large", "E_large_MAG"]
MARKER_RICH = ["donor_ammonia_oxidation", "acceptor_sulfate", "acceptor_CO2_methanogen_acetogen", "acceptor_oxygen"]
WEAK = ["donor_iron_oxidation", "acceptor_iron", "donor_other", "acceptor_other"]


def main():
    m = pd.read_csv(RES / "metrics_macro.csv")
    pt = pd.read_csv(RES / "metrics_per_trait.csv")
    pc = pd.read_csv(RES / "paired_comparisons.csv")
    g = m[m.split_level == "genus"].set_index("arm").macro_auprc
    f = m[m.split_level == "family"].set_index("arm").macro_auprc
    rows = []

    emb_present = [a for a in EMB if a in g.index]
    emb_max = g[emb_present].max()
    rows.append(dict(prediction="P1", statement="Arm B (KO+model) beats or ties every embedding arm on macro AUPRC "
                                                "(held-out genus)",
                     evidence=f"B2 {g['B2']:.4f} vs best embedding arm {g[emb_present].idxmax()} {emb_max:.4f} "
                              f"(arms compared: {emb_present})",
                     outcome="SUPPORTED" if g["B2"] >= emb_max else
                     "NOT SUPPORTED (statistical tie: the protocol's preferred encoder arm edges ahead)"))

    d = pc[(pc.split_level == "genus") & (pc.contrast == "B2 - B1r") & (pc.trait == "MACRO")]
    b1, b2 = g.get("B1", np.nan), g["B2"]
    rows.append(dict(prediction="P2", statement="Arm B2 (852-KO energy panel) >= B1 (84 provided markers)",
                     evidence=f"B2 {b2:.3f} vs B1 {b1:.3f} (delta {b2 - b1:+.4f}); B2-B1r macro delta "
                              f"{d.delta_auprc.iloc[0]:+.4f} [{d.ci_low.iloc[0]:.4f}, {d.ci_high.iloc[0]:.4f}]",
                     outcome="NOT SUPPORTED (statistical tie, point estimate below B1)"))

    b = pt[(pt.arm == "B2") & (pt.split_level == "genus")].set_index("trait").auprc.sort_values()
    worst4 = set(b.index[:4])
    rows.append(dict(prediction="P3", statement="Iron oxidation, iron reduction and the two 'other' categories are "
                                                "the worst traits for every arm",
                     evidence=f"four lowest B2 AUPRC traits: {sorted(worst4)}",
                     outcome="PARTLY SUPPORTED" if len(worst4 & set(WEAK)) >= 3 else "NOT SUPPORTED"))

    ec = pc[(pc.split_level == "genus") & (pc.contrast == "E - C") & (pc.trait == "MACRO")]
    tie = bool(ec.ci_low.iloc[0] <= 0 <= ec.ci_high.iloc[0])
    rows.append(dict(prediction="P4", statement="Arm E (Bacformer) ~ arm C (pooled ESM-2): synteny adds little",
                     evidence=f"E-C macro delta {ec.delta_auprc.iloc[0]:+.4f} [{ec.ci_low.iloc[0]:.4f}, "
                              f"{ec.ci_high.iloc[0]:.4f}]",
                     outcome="SUPPORTED" if tie else "NOT SUPPORTED"))

    fam = pt[(pt.arm == "B2") & (pt.split_level == "family")].set_index("trait")
    rich_ok = [t for t in MARKER_RICH if fam.delta_vs_null_q_bh[t] < 0.05 and fam.delta_vs_null[t] > 0]
    weak_ok = [t for t in WEAK if fam.delta_vs_null_q_bh[t] < 0.05 and fam.delta_vs_null[t] > 0]
    rows.append(dict(prediction="P5", statement="On held-out family the model beats the taxonomic null on "
                                                "marker-rich traits but not on iron / 'other' traits",
                     evidence=f"marker-rich beating null (BH q<0.05): {rich_ok} of {MARKER_RICH}; "
                              f"weak traits beating null: {weak_ok}",
                     outcome=f"PARTLY SUPPORTED ({len(rich_ok)}/{len(MARKER_RICH)} marker-rich traits beat the null; "
                             f"but so do {len(weak_ok)} of the {len(WEAK)} traits predicted not to: {weak_ok})"
                     if weak_ok else f"SUPPORTED ({len(rich_ok)}/{len(MARKER_RICH)} marker-rich, 0 weak)"))

    arms = [a for a in g.index if not a.startswith("null")]
    drops = {a: g[a] - f[a] for a in arms if a in f.index}
    n_drop = sum(v > 0 for v in drops.values())
    rows.append(dict(prediction="P6", statement="Performance drops from genus- to family-level hold-out for every "
                                                "learned arm",
                     evidence=f"{n_drop}/{len(drops)} arms drop; exceptions: "
                              f"{[a for a, v in drops.items() if v <= 0]}",
                     outcome="MOSTLY SUPPORTED" if n_drop >= len(drops) - 2 else "NOT SUPPORTED"))

    out = pd.DataFrame(rows)
    out.to_csv(RES / "preregistered_predictions_scored.csv", index=False)
    print(out.to_string())


if __name__ == "__main__":
    main()
