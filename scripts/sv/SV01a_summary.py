"""SV01a: write results/sv/SV01a/SUMMARY.md (rule 1.7) from the saved SV01a files only.
Line 1 (not done / substituted) is NOT_DONE below.
Run: .venv/bin/python -I scripts/sv/SV01a_summary.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV01a"
NOT_DONE = (
    "Not done / substituted: (1) docs/sv/briefs/SV01a.md carried rules 1-3 but not the SV01 text, so the step text was "
    "read from docs/sv/RUN_SV1028.md section 4 (SV01); (2) the run file gives no rule for missing predictor values: "
    "pre-declared median imputation within the cluster plus a 'not reported' indicator in block X (pH missing for "
    "{ph7} cl7 and {ph4} cl4 chains), so every model of a cluster uses the same rows; (3) 'mean B' = the entry mean B of "
    "Q05/K05 (the QC-tier and Y14 value), not the chain C-alpha mean B, which averages over the outcome residues; "
    "(4) entry level = one chain per entry (the entry's first matrix chain in the cluster, Y14 convention); "
    "(5) the step's drop rule was extended, before fitting, to exact collinearity and to joint separation in the full "
    "model (drop the predictor with the largest |log OR| while the fit separates); (6) the [DEFAULT] verdict rule names "
    "no level, so it was applied at chain and at entry level and the step verdict is 'artifact'/'biology' only when "
    "both levels agree; (7) chain-level errors: logit MLE with entry-clustered sandwich errors (= GEE, independence "
    "working correlation); the fraction outcome uses the fractional-logit QMLE (statsmodels Logit). Lattice contacts "
    "and the mobility outcome belong to SV01b and were not attempted here.")


def rd(name):
    return pd.read_csv(OUT / name, low_memory=False)


def size(p):
    s = p.stat().st_size
    return f"{s / 1e3:.1f} kB" if s < 1e6 else f"{s / 1e6:.2f} MB"


def main():
    nums = json.loads((OUT / "numbers.json").read_text())
    sc = rd("h12_status_counts.csv")
    pe = rd("partial_ends.csv")
    prof = rd("completeness_profile.csv")
    cs = rd("concordance_summary.csv")
    mb = rd("model_blocks.csv")
    orr = rd("odds_ratios.csv")
    vt = rd("verdict_table.csv")
    verdict_md = (OUT / "verdict.md").read_text().strip()

    scl = rd("predictor_scaling.csv")
    phm = scl[(scl.predictor == "ph") & scl.n_missing.notna()].set_index("cluster").n_missing
    L = [NOT_DONE.format(ph7=int(phm[7]), ph4=int(phm[4])), "", "# SV01a - H12 completeness, first pass", ""]
    step = nums["step_verdict"]
    lv = nums["level_verdicts"]
    L += [f"**Verdict ([DEFAULT] rule, H12 complete on 538-547): {step['complete']}** - chain level "
          f"{lv['complete|chain']}, entry level {lv['complete|entry']}. Sensitivity (fraction of 538-547 observed): "
          f"{step['fraction_observed']} (chain {lv['fraction_observed|chain']}, entry {lv['fraction_observed|entry']}).", ""]
    L += ["## Inputs", "",
          f"- results/sv/SV00/chain_table.csv and observed_mask.npy (tracker card SV00 out; SV00 DONE in tracker/esr1.log); "
          f"{nums['n_chains']} chains, {nums['n_entries']} entries, remediated chains excluded: 0.",
          f"- Outcome from the mask (all 10 C-alpha of 538-547 observed) agrees with the S16 h12_label for "
          f"{nums['mask_vs_s16_label_agree']}/{nums['n_chains']} chains.", ""]

    L += ["## 1. Distribution", "", "H12 status (h12_status_counts.csv; entry level = first matrix chain of each entry in the group):", "",
          "| level | group | n | complete | partial | not modelled | fraction complete |", "|---|---|---|---|---|---|---|"]
    for _, r in sc[sc.group_type.isin(["cluster", "all", "mutation"])].iterrows():
        L.append(f"| {r.level} | {r.group} | {r.n} | {r.complete} | {r.partial} | {r.not_modelled} | {r.frac_complete:.3f} |")
    L.append("")
    cm = sc[(sc.group_type == "cluster_mutation") & (sc.level == "chain")]
    L.append("Within cluster by construct (chain level, fraction complete): " + "; ".join(
        f"{r.group} {r.complete}/{r.n}" for r in cm.itertuples()) + ".")
    L.append("")
    side = pe.side.value_counts()
    bycl = pd.crosstab(pe.cluster_group, pe.side)
    L += [f"Partial H12 ({len(pe)} chains, partial_ends.csv): " + ", ".join(f"{k} {v}" for k, v in side.items()) +
          "; by cluster " + "; ".join(f"{c}: " + ", ".join(f"{k} {int(v)}" for k, v in row.items() if v)
                                      for c, row in bycl.iterrows()) +
          f". No partial chain lacks only N-side residues. Missing at the C end: median "
          f"{pe.n_missing_C_end.median():.0f}, max {pe.n_missing_C_end.max()} residues; internal gaps in "
          f"{int((pe.n_missing_internal > 0).sum())} chain(s). Last observed residue of partial chains: " +
          ", ".join(f"{int(k)}: {v}" for k, v in pe.last_observed_residue.value_counts().sort_index().items()) + ".", ""]
    pc = prof[(prof.group_type == "cluster") & prof.group.isin(["cl7", "cl4"])].set_index(["group", "residue"])
    w = prof[(prof.group_type == "cluster") & prof.group.isin(["cl7", "cl4"]) & prof.residue.between(520, 550)]
    mins = w.loc[w.groupby("group").frac_observed.idxmin()]
    loop = w[w.residue.between(525, 537)]
    lmin = loop.loc[loop.groupby("group").frac_observed.idxmin()]
    L += ["Per-residue observation (completeness_profile.csv, chain fraction): " +
          "; ".join(f"{g}: 538 {pc.loc[(g, 538), 'frac_observed']:.2f}, 545 {pc.loc[(g, 545), 'frac_observed']:.2f}, "
                    f"546 {pc.loc[(g, 546), 'frac_observed']:.2f}, 547 {pc.loc[(g, 547), 'frac_observed']:.2f}, "
                    f"550 {pc.loc[(g, 550), 'frac_observed']:.2f}" for g in ("cl7", "cl4")) +
          ". H11-H12 loop minimum (525-537): " + "; ".join(
              f"{r.group} {r.frac_observed:.2f} at {int(r.residue)}" for r in lmin.itertuples()) + ".", ""]

    L += ["## 2. Within-entry concordance (concordance.csv, concordance_summary.csv)", "",
          "| measure | n | k | fraction |", "|---|---|---|---|"]
    for r in cs.itertuples():
        L.append(f"| {r.measure} | {r.n} | {r.k} | {r.fraction:.3f} |")
    L.append("")

    L += ["## 3. Block models (model_blocks.csv, odds_ratios.csv)", "",
          "McFadden R2 gain, outcome H12 complete (Q quality, X experiment, B biology; Q+X combined):", "",
          "| cluster | level | n | incomplete | block | first | last | test (entered last) | q (BH) |",
          "|---|---|---|---|---|---|---|---|---|"]
    m = mb[mb.outcome == "complete"]
    for (c, lvl, blk), g in m.groupby(["cluster", "level", "block"], sort=False):
        f = g[g.order == "first"].iloc[0]
        la = g[g.order == "last"].iloc[0]
        L.append(f"| cl{c} | {lvl} | {f.n} | {f.n_incomplete} | {blk} | {f.mcfadden_gain:.3f} | {la.mcfadden_gain:.3f} | "
                 f"{la.test.split(' (')[0]} {la.statistic:.1f} (df {la.df:.0f}), p {la.p_value:.2g} | {la.q_value:.2g} |")
    L.append("")
    L.append(f"BH: one family over every p-value of the step ({nums['n_tests_bh']} tests).")
    L.append("")
    L += ["Sensitivity outcome (fraction observed), gain entered first, Q+X vs B: " + "; ".join(
        f"cl{r.cluster} {r.level} {r.gain_QX:.3f} vs {r.gain_B:.3f}"
        for r in vt[(vt.outcome == "fraction_observed") & (vt.order == "first")].itertuples()) + ".", ""]
    sig = orr[(orr.outcome == "complete") & (orr.q_value < 0.05)].sort_values(["cluster", "level", "q_value"])
    L += ["Odds ratios of the full model with BH q < 0.05 (complete outcome): " + ("; ".join(
        f"cl{r.cluster} {r.level} {r.label} ({r.block}, {r.unit}) OR {r.odds_ratio:.3g} [{r.ci_low:.3g}, {r.ci_high:.3g}], "
        f"q {r.q_value:.2g}" for r in sig.itertuples()) or "none") + ".", ""]
    dr = orr[(orr.outcome == "complete") & (orr.status != "fitted")]
    L += ["Dropped predictors (complete outcome), cluster/level: " + "; ".join(
        f"cl{c} {lvl}: " + ", ".join(f"{r.predictor} [{r.block}] {r.status.replace('dropped: ', '')}" for r in g.itertuples())
        for (c, lvl), g in dr.groupby(["cluster", "level"], sort=False)) + ".", ""]
    L += ["## Verdict (verdict.md)", "", verdict_md, "",
          "Reading: the chain-level cl7 comparison is a near tie (B exceeds Q+X by "
          f"{float(vt[(vt.outcome == 'complete') & (vt.level == 'chain') & (vt.cluster == 7) & (vt.order == 'first')].eval('gain_B - gain_QX').iloc[0]):.4f}"
          "), cl7 has few incomplete chains for the number of predictors, and the drop rule removed several B predictors "
          "that separate the outcome (their information is not credited to B). SV01b refits with lattice contacts.", ""]

    L += ["## Files", ""]
    files = sorted(OUT.glob("*")) + sorted((ROOT / "figures/sv/SV01a").glob("*")) + \
        sorted((ROOT / "figures/panels/sv/SV01a").glob("*")) + \
        [ROOT / "scripts/sv/SV01a_completeness.py", ROOT / "scripts/sv/SV01a_figures.py",
         ROOT / "scripts/sv/SV01a_summary.py", ROOT / "scripts/sv/common/sv_fig.py"]
    for p in files:
        if p.name == "SUMMARY.md":
            continue
        L.append(f"- {p.relative_to(ROOT)} ({size(p)})")
    L.append("- results/sv/SV01a/SUMMARY.md (this file)")
    (OUT / "SUMMARY.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
