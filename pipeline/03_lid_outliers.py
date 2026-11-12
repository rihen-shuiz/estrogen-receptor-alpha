"""Stage 03 - outliers of the lid-to-H11 distance (H12 centroid to the H11 axis, core-v3 frame) and what they are.

Declared before looking at individual chains:
  * Outlier rule per group (cl7, cl4): Tukey fences, below Q1 - 1.5 IQR or above Q3 + 1.5 IQR.
  * Minor clusters ("grey"): every chain is examined (14 chains; 10 have an H12 axis).
  * cl4 tails: chains with |robust z| >= 2 (z = (x - median) / (1.4826 MAD)) are examined as "cl4 extremes" whether
    or not they pass the Tukey rule, because cl4 is the group the meeting singled out.
  * An outlier "exceeds ordinary variation" when |x - group median| is larger than both the within-group SD and the
    95th percentile of copy-to-copy differences in that group (common.py).
  * What drives the spread inside cl7 and cl4: entry-level (mean over the entry's chains in the group) Mann-Whitney
    tests for binary covariates and Spearman for numeric ones, BH over all tests of the step; an effect is reported
    as real only if q < 0.05 AND the difference of entry medians exceeds the ordinary-variation threshold.
  * Action class vs outlier status: Fisher exact per group, chain and entry level (uncurated ChEMBL-mechanism classes only).
Outputs (results/03_lid_outliers/): lid_distance_by_group.csv, lid_outliers.csv, lid_distance_per_chain.csv, lid_spread_drivers.csv,
lid_outlier_class_table.csv, numbers_lid_outliers.json.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

OUT = c.out_dir("03_lid_outliers")

COL = "dist_H12c_H11axis_A"


def why_unusual(r, grp: pd.DataFrame, t: pd.DataFrame, med: float, thr: float) -> str:
    """Plain-language reasons: attributes that are rare in the chain's own group, plus copy agreement."""
    out = []
    if pd.notna(r[COL]):
        dev = r[COL] - med
        out.append(f"lid {r[COL]:.2f} A, {dev:+.2f} A from the group median "
                   f"({'beyond' if abs(dev) > thr else 'within'} ordinary variation {thr:.2f} A)")
    else:
        out.append("no H12 axis (fewer than 6 helical H12 C-alpha)")
    share_pep = grp.peptide_bound.mean()
    if bool(r.peptide_bound) != (share_pep >= 0.5):
        out.append(f"{'peptide' if r.peptide_bound else 'no peptide'} (group: {share_pep:.0%} with peptide)")
    m = grp.mutation_status.value_counts(normalize=True)
    if m.get(r.mutation_status, 0) < 0.10:
        out.append(f"construct {r.construct_mutations} ({m.get(r.mutation_status, 0):.0%} of group are "
                   f"{r.mutation_status})")
    sg = grp.space_group.value_counts(normalize=True)
    if sg.get(r.space_group, 0) < 0.05:
        out.append(f"space group {r.space_group} ({sg.get(r.space_group, 0):.0%} of group)")
    if r.h12_label != "complete":
        out.append(f"H12 {r.h12_label} ({int(r.n_h12_ca)} of 10 C-alpha 538-547)")
    if pd.notna(r.resolution_A) and r.resolution_A > grp.resolution_A.quantile(0.95):
        out.append(f"resolution {r.resolution_A} A (> group 95th pct)")
    if pd.notna(r.ligand_heavy_atoms):
        q05, q95 = grp.ligand_heavy_atoms.quantile([0.05, 0.95])
        if r.ligand_heavy_atoms < q05 or r.ligand_heavy_atoms > q95:
            out.append(f"ligand {int(r.ligand_heavy_atoms)} heavy atoms (outside group 5-95th pct)")
    if bool(r.consensus_outlier):
        out.append(f"also a consensus RMSD outlier ({r.outlier_verdict}: {r.outlier_reasons})")
    cop = t[(t.entry_id == r.entry_id) & (t.chain_key != r.chain_key)]
    if len(cop):
        vals = ", ".join(f"{k.split('|')[1]} {v:.2f}" if pd.notna(v) else f"{k.split('|')[1]} n/a"
                         for k, v in zip(cop.chain_key, cop[COL]))
        out.append(f"other copies in the crystal: {vals}")
    return "; ".join(out)


def lid_position(x, med: dict, thr: dict) -> str:
    """Where a lid distance sits relative to the two main groups (each group's median +- its ordinary threshold)."""
    if pd.isna(x):
        return "no H12 axis"
    if x < med["cl7"] - thr["cl7"]:
        return "below the cl7 range"
    if x <= med["cl7"] + thr["cl7"]:
        return "cl7-like"
    if x < med["cl4"] - thr["cl4"]:
        return "between cl7 and cl4"
    if x <= med["cl4"] + thr["cl4"]:
        return "cl4-like"
    return "beyond the cl4 range"


def main():
    t = c.chains()
    rows, flagged = [], []
    main_med = {g: float(t[t.group == g][COL].median()) for g in ("cl7", "cl4")}
    main_thr = {g: c.ordinary_variation(t, COL, g)["threshold"] for g in ("cl7", "cl4")}
    for g in ("cl7", "cl4", "minor"):
        d = t[t.group == g]
        x = d[COL].dropna()
        q1, q3 = x.quantile([0.25, 0.75])
        iqr = q3 - q1
        lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
        med = float(x.median())
        mad = float(1.4826 * np.median(np.abs(x - med)))
        ov = c.ordinary_variation(t, COL, g)
        tuk = d[(d[COL] < lo) | (d[COL] > hi)]
        rows.append(dict(group=g, n_chains=len(d), n_with_value=len(x), n_entries=d.entry_id.nunique(), median=med,
                         q1=q1, q3=q3, iqr=iqr, fence_low=lo, fence_high=hi, sd=float(x.std()), mad_scaled=mad,
                         min=float(x.min()), max=float(x.max()), copy_p95=ov["copy_p95"],
                         n_copy_pairs=ov["n_copy_pairs"], ordinary_threshold=ov["threshold"], n_tukey=len(tuk)))
        if g == "minor":
            sel = d.assign(reason_flagged="minor cluster (grey)")
        else:
            z = (d[COL] - med) / mad
            sel = d[(d[COL] < lo) | (d[COL] > hi) | (z.abs() >= 2 if g == "cl4" else False)].copy()
            sel["reason_flagged"] = np.where((sel[COL] < lo) | (sel[COL] > hi), "Tukey outlier", "cl4 extreme |z|>=2")
        for _, r in sel.iterrows():
            dev = r[COL] - med if pd.notna(r[COL]) else np.nan
            near, near_dev, near_thr = g, dev, ov["threshold"]
            if g == "minor" and pd.notna(r[COL]):
                # a grey chain is compared with the nearer main group, not with the pooled minor median
                cand = [(abs(r[COL] - main_med[m]), m) for m in ("cl7", "cl4")]
                near = min(cand)[1]
                near_dev, near_thr = r[COL] - main_med[near], main_thr[near]
            lid_state = lid_position(r[COL], main_med, main_thr)
            flagged.append(dict(
                chain_key=r.chain_key, entry_id=r.entry_id, group=g, cluster_1p5=int(r.cluster_1p5),
                reason_flagged=r.reason_flagged, lid_to_h11_A=r[COL], deviation_A=dev,
                robust_z=(dev / mad) if pd.notna(dev) else np.nan,
                beyond_ordinary=bool(abs(near_dev) > near_thr) if pd.notna(near_dev) else None,
                compared_with=near, deviation_from_compared_A=near_dev, lid_state=lid_state,
                ligand=r.ligand_id_primary, ligand_name=r.ligand_name, ligand_heavy_atoms=r.ligand_heavy_atoms,
                basic_amine=r.ligand_basic_amine, action_class=r.action_class, peptide_bound=r.peptide_bound,
                peptide=r.peptide_class_entry, mutations=r.construct_mutations, mutation_status=r.mutation_status,
                additives=r.additive_codes, resolution_A=r.resolution_A, r_free=r.r_free, space_group=r.space_group,
                h12_label=r.h12_label, n_h12_ca=r.n_h12_ca, h12_mobility_z=r.mobility_h12_z,
                lattice_sym_atoms_528_550=r.sym_atoms, consensus_outlier=r.consensus_outlier,
                rmsd_outlier_reason=r.outlier_reasons, why_unusual=why_unusual(r, t[t.group == near] if g == "minor" and near != "minor" else d, t,
                                        main_med.get(near, med), near_thr)))
    # compact per-chain table behind the lid figures
    t[["chain_key", "entry_id", "group", COL]].rename(columns={COL: "lid_to_h11_A"}).assign(
        examined=t.chain_key.isin({r["chain_key"] for r in flagged})).to_csv(OUT / "lid_distance_per_chain.csv",
                                                                            index=False)
    summ = pd.DataFrame(rows)
    summ.to_csv(OUT / "lid_distance_by_group.csv", index=False)
    fl = pd.DataFrame(flagged)
    fl.to_csv(OUT / "lid_outliers.csv", index=False)

    # ---------------- what drives the spread inside each main group (entry level)
    drivers, pv = [], []
    for g in ("cl7", "cl4"):
        d = t[(t.group == g) & t[COL].notna()].copy()
        d["L536S"] = d.mut_L536S if "mut_L536S" in d else d.mutation_status.eq("L536S")
        d["Y537S"] = d.mutation_status.eq("Y537S")
        d["agonist_vs_SERM"] = d.action_class.map({"agonist": True, "SERM": False})
        d["h12_complete"] = d.h12_label.eq("complete")
        thr = summ.set_index("group").loc[g, "ordinary_threshold"]
        ent = d.groupby("entry_id").agg(y=(COL, "mean"), peptide_bound=("peptide_bound", "max"),
                                        L536S=("L536S", "max"), Y537S=("Y537S", "max"),
                                        basic_amine=("ligand_basic_amine", "max"),
                                        agonist_vs_SERM=("agonist_vs_SERM", "first"),
                                        h12_complete=("h12_complete", "mean"),
                                        heavy_atoms=("ligand_heavy_atoms", "mean"),
                                        resolution_A=("resolution_A", "mean"), sym_atoms=("sym_atoms", "mean"),
                                        sub_1p0=("cut_1.0", lambda s: s.mode().iloc[0]))
        for cov in ("peptide_bound", "L536S", "Y537S", "basic_amine", "agonist_vs_SERM", "h12_complete"):
            e = ent.dropna(subset=[cov])
            if cov == "h12_complete":
                e = e.assign(h12_complete=e.h12_complete >= 0.5)
            a, b = e[e[cov].astype(bool)].y, e[~e[cov].astype(bool)].y
            if len(a) >= 3 and len(b) >= 3:
                pval = mannwhitneyu(a, b).pvalue
                diff = float(a.median() - b.median())
            else:
                pval, diff = np.nan, (float(a.median() - b.median()) if len(a) and len(b) else np.nan)
            drivers.append(dict(group=g, covariate=cov, kind="binary (true minus false, entry medians)",
                                n_true=len(a), n_false=len(b), effect=diff, p=pval, threshold=thr))
            pv.append(pval)
        # largest 1.0 A sub-cluster vs the rest of the group
        top = ent.sub_1p0.value_counts().index[0]
        a, b = ent[ent.sub_1p0 == top].y, ent[ent.sub_1p0 != top].y
        pval = mannwhitneyu(a, b).pvalue if len(a) >= 3 and len(b) >= 3 else np.nan
        drivers.append(dict(group=g, covariate=f"1.0 A sub-cluster {top} vs rest", kind="binary",
                            n_true=len(a), n_false=len(b), effect=float(a.median() - b.median()), p=pval,
                            threshold=thr))
        pv.append(pval)
        for cov in ("heavy_atoms", "resolution_A", "sym_atoms"):
            e = ent.dropna(subset=[cov])
            rho, pval = spearmanr(e[cov], e.y)
            # effect in A: difference of entry medians between the top and bottom quartile of the covariate
            q1, q3 = e[cov].quantile([0.25, 0.75])
            eff = float(e[e[cov] >= q3].y.median() - e[e[cov] <= q1].y.median())
            drivers.append(dict(group=g, covariate=cov, kind=f"numeric (rho {rho:.2f}; top minus bottom quartile)",
                                n_true=len(e), n_false=0, effect=eff, p=pval, threshold=thr))
            pv.append(pval)
    dr = pd.DataFrame(drivers)
    dr["q"] = c.bh(pv)
    dr["real"] = (dr.q < 0.05) & (dr.effect.abs() > dr.threshold)
    dr.to_csv(OUT / "lid_spread_drivers.csv", index=False)

    # ---------------- action class vs outlier status
    out_keys = set(fl[fl.reason_flagged == "Tukey outlier"].chain_key)
    ct_rows = []
    for g in ("cl7", "cl4", "minor"):
        d = t[(t.group == g) & t[COL].notna()].copy()
        d["outlier"] = d.chain_key.isin(out_keys) if g != "minor" else True
        for cls in ("agonist", "SERM", "SERD", "unlabelled"):
            ct_rows.append(dict(group=g, action_class=cls, level="chain",
                                outlier=int(((d.action_class == cls) & d.outlier).sum()),
                                not_outlier=int(((d.action_class == cls) & ~d.outlier).sum())))
            e = d.groupby("entry_id").agg(cls=("action_class", lambda s: (s == cls).any()), o=("outlier", "max"))
            ct_rows.append(dict(group=g, action_class=cls, level="entry", outlier=int((e.cls & e.o).sum()),
                                not_outlier=int((e.cls & ~e.o).sum())))
    ctab = pd.DataFrame(ct_rows)
    tests = []
    for (g, lvl), d in ctab[ctab.group != "minor"].groupby(["group", "level"]):
        d = d.set_index("action_class")
        lab = d.loc[["agonist", "SERM", "SERD"]].sum()
        un = d.loc["unlabelled"]
        tab = [[lab.outlier, lab.not_outlier], [un.outlier, un.not_outlier]]
        tests.append(dict(group=g, level=lvl, test="labelled (any class) vs unlabelled",
                          table=str(tab), p=fisher_exact(tab)[1] if np.array(tab).sum() else np.nan))
    tt = pd.DataFrame(tests)
    tt["q"] = c.bh(tt.p)
    ctab.to_csv(OUT / "lid_outlier_class_table.csv", index=False)
    tt.to_csv(OUT / "lid_outlier_class_tests.csv", index=False)

    s = summ.set_index("group")
    tuk = fl[fl.reason_flagged == "Tukey outlier"]
    c.write_numbers(OUT, {
        "median_A": s["median"].to_dict(), "sd_A": s["sd"].to_dict(), "copy_p95_A": s["copy_p95"].to_dict(),
        "threshold_A": s["ordinary_threshold"].to_dict(), "n_with_value": s["n_with_value"].to_dict(),
        "n_tukey": s["n_tukey"].to_dict(),
        "tukey_cl7": tuk[tuk.group == "cl7"].chain_key.tolist(),
        "tukey_cl7_beyond_ordinary": tuk[(tuk.group == "cl7") & tuk.beyond_ordinary.astype(bool)].chain_key.tolist(),
        "cl4_extremes": fl[fl.reason_flagged == "cl4 extreme |z|>=2"].chain_key.tolist(),
        "minor_with_value": int(fl[(fl.group == "minor")].lid_to_h11_A.notna().sum()),
        "minor_lid_state": fl[fl.group == "minor"].lid_state.value_counts().to_dict(),
        "minor_intermediate": fl[(fl.group == "minor") & (fl.lid_state == "intermediate")].chain_key.tolist(),
        "drivers_real": dr[dr.real][["group", "covariate", "effect", "q"]].to_dict("records"),
        "class_tests_min_q": float(tt.q.min()) if len(tt) else None,
        "h12_h3_angle_median": t.groupby("group").angle_H12_H3_deg.median().to_dict(),
    })
    print(summ.to_string())
    print(fl[["chain_key", "group", "reason_flagged", "lid_to_h11_A", "deviation_A", "beyond_ordinary",
              "action_class", "peptide_bound", "mutation_status", "ligand"]].to_string())
    print(dr.to_string())
    print(tt.to_string())


if __name__ == "__main__":
    main()
