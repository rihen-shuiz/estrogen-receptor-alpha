"""SV07: write results/sv/SV07/SUMMARY.md (rule 1.7) and env.txt from the saved SV07 files only.
Run after SV07_lattice.py, SV07_spikes.py and SV07_figures.py: .venv/bin/python -B scripts/sv/SV07_summary.py
"""
import json
import platform
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV07"
FIG = [ROOT / "figures/sv/SV07/SV07_spike_profile.png", ROOT / "figures/sv/SV07/SV07_spike_profile.svg",
       ROOT / "figures/panels/sv/SV07/SV07_spike_profile__profile.png",
       ROOT / "figures/panels/sv/SV07/SV07_spike_profile__observation.png"]
SCRIPTS = ["SV07_spikes.py", "SV07_lattice.py", "SV07_figures.py", "SV07_summary.py"]


def size(p):
    s = p.stat().st_size
    return f"{s / 1e3:.1f} kB" if s < 1e6 else f"{s / 1e6:.2f} MB"


def f2(v, nd=2):
    return "-" if pd.isna(v) else f"{v:.{nd}f}"


def env():
    import gemmi
    import matplotlib
    import numpy
    import scipy
    import statsmodels
    now = datetime.now(timezone(timedelta(hours=5))).strftime("%Y-%m-%dT%H:%M+0500")
    lines = [f"# SV07 environment, {now}",
             f"python {platform.python_version()} at {sys.executable}",
             f"platform {platform.platform()}",
             f"numpy {numpy.__version__}", f"pandas {pd.__version__}", f"scipy {scipy.__version__}",
             f"statsmodels {statsmodels.__version__}", f"gemmi {gemmi.__version__}",
             f"matplotlib {matplotlib.__version__}",
             "code reused unchanged: scripts/sv/SV05_core_v4.y10_profile (Y10 formulas), "
             "scripts/sv/common/core_frame.load_ca_core, scripts/sv/common/lattice_contacts.py (SV01b module)",
             "nothing installed, upgraded or downgraded by SV07; no network used"]
    (OUT / "env.txt").write_text("\n".join(lines) + "\n")


def main():
    env()
    n = json.loads((OUT / "numbers.json").read_text())
    sp = pd.read_csv(OUT / "spikes.csv")
    vd = pd.read_csv(OUT / "verdict.csv")
    top = pd.read_csv(OUT / "top_entries.csv")
    cov = pd.read_csv(OUT / "covariate_tests.csv")
    prof = pd.read_csv(OUT / "restricted_profile.csv")
    lat = json.loads((OUT / "lattice_checks.json").read_text())
    rv = n["restricted_vs_y10"]
    vc = vd.verdict.value_counts()
    gen = vd[vd.verdict == "genuine loop difference"].position.tolist()
    few = vd[vd.verdict == "driven by <= 3 entries"]
    med_only = few[few.exceeds == "medoid"].position.tolist()
    few_mean = few[few.exceeds != "medoid"].position.tolist()
    m = sp.set_index("position")
    n_sparse_medoid = int(((vd.exceeds.str.contains("medoid"))).sum())

    L = [("Not done / substituted: (1) PARTLY BLOCKED - the 'sparse-sample effect (gone when restricted)' rule is "
          "degenerate for cluster-mean spikes: Y10 already averages, at each position, only the chains that observe it "
          "(nan-mean), so the cluster means recomputed on the chains observing the position in both clusters equal Y10 "
          f"at all 246 positions (max |diff| {rv['max_abs_diff_d_mean_A']:.1g} A) and no cluster-mean spike can be 'gone "
          "when restricted'; the rule could still act through the medoid measure (complete-H12 medoids re-picked among "
          f"the chains observing the position; {n_sparse_medoid} spikes exceed on the medoid). Tried: the other reading "
          f"of 'in both clusters' (entries holding chains of both clusters) - only {n['n_entries_with_cl7_and_cl4_chains']} "
          f"of {n['n_entries_cl7_or_cl4']} cl7/cl4 entries do, too few to recompute a profile; not used. "
          "(2) The run file gives no operational rule for 'driven by <= 3 entries', 'numbering error', 'persists across "
          "many entries' or the covariate tests; declared in scripts/sv/SV07_spikes.py before any diagnostic was "
          "computed: greedy removal of up to 3 entries (each time the entry whose removal lowers d_mean most); numbering "
          "flag = residue differs from UniProt with no deposited explanation (residue_status unexplained_difference) or "
          "UniProt-author offset not constant over the LBD / not from SIFTS (SV00); covariate tests = OLS of the signed "
          "C-alpha coordinate along the cl7->cl4 displacement on cluster + covariate, entry-clustered Wald chi2 (chain) "
          f"and nested F (entry x cluster units), one BH family over {n['n_covariate_tests']} tests, 'associated' = "
          "q < 0.05 at both levels; primary verdict at chain level (the Y10 unit), entry level reported. (3) SV01b saved "
          "B z only for 528-547; the SV01b formula was applied to all positions and checked against that file "
          f"({n['b_z_check_vs_SV01b']['n_compared']} residues, max |diff| {n['b_z_check_vs_SV01b']['max_abs_diff']:.1g}). "
          "(4) SSE context is SV02's pydssp fallback ('H' = any helix). (5) Supporting files not named in the run file: "
          "top_entries.csv, restricted_profile.csv, covariate_tests.csv, lattice_contacts_spikes_per_residue.csv, "
          "lattice_contacts_spikes_per_chain.csv, lattice_checks.json, numbers.json.")]
    L += ["", "# SV07 - Spike audit of the Y10 displacement profile (slide 5)", ""]
    L.append(f"**Verdict:** {len(vd)} positions outside 528-550 meet the [DEFAULT] spike rule (Y10 cluster-mean or "
             f"medoid-pair C-alpha distance > 1.0 A). **{vc.get('genuine loop difference', 0)} are genuine loop "
             f"differences** ({', '.join(map(str, gen))}): the cl7-cl4 distance stays > 1.0 A after removing the 3 "
             f"most influential entries, at chain and entry level. **{len(few)} are driven by <= 3 entries**: "
             f"{', '.join(map(str, med_only))} exceed only on the medoid pair (cluster means "
             f"{f2(m.loc[med_only].d_mean_y10_A.min())}-{f2(m.loc[med_only].d_mean_y10_A.max())} A), and "
             f"{', '.join(map(str, few_mean))} drops to "
             f"{f2(m.loc[few_mean[0]].d_mean_minus_1_entries_A) if few_mean else '-'} A without one entry "
             f"({m.loc[few_mean[0]].entries_removed_greedy.split(';')[0] if few_mean else '-'}). "
             f"No spike is a numbering error (0 unexplained residue differences at the spike positions; removing the "
             f"{n['n_chains_offset_flag']} chains with a non-constant or non-SIFTS offset changes d_mean by <= "
             f"{f2((sp.d_mean_without_flagged_A - sp.d_mean_restricted_A).abs().max(), 2)} A) and none is a "
             f"sparse-sample effect under the rule as written (see line 1). Chain- and entry-level verdicts agree "
             f"for {int(vd.levels_agree.sum())} of {len(vd)} spikes.")
    # covariate pattern among the genuine spikes (from covariate_tests.csv / verdict.csv)
    assoc = {c: [p for p in gen if c in str(vd[vd.position == p].covariates_associated.iloc[0])]
             for c in ("ligand_chemotype", "dimer_partner_source", "peptide_bound")}
    none = [p for p in gen if str(vd[vd.position == p].covariates_associated.iloc[0]).startswith("none")]
    ch = cov[cov.level == "chain"].set_index(["position", "covariate"])
    c335 = ch.loc[(335, "ligand_chemotype")] if (335, "ligand_chemotype") in ch.index else None
    c467 = ch.loc[(467, "peptide_bound")] if (467, "peptide_bound") in ch.index else None
    txt = (f" Among the genuine spikes, ligand chemotype is associated (q < 0.05 at chain and entry level) at "
           f"{', '.join(map(str, assoc['ligand_chemotype'])) or 'none'}, dimer-partner source at "
           f"{', '.join(map(str, assoc['dimer_partner_source'])) or 'none'} and peptide_bound at "
           f"{', '.join(map(str, assoc['peptide_bound'])) or 'none'}; {', '.join(map(str, none)) or 'none'} show no "
           "association.")
    if c335 is not None:
        txt += (f" At 335 the cl4-cl7 coefficient falls from {c335.b_cl4_without_cov_A:.2f} to "
                f"{c335.b_cl4_with_cov_A:.2f} A with chemotype in the model (mean signed coordinate "
                f"{c335.chain_mean_s_by_cluster_level_A}).")
    if c467 is not None:
        txt += (f" At 467 it falls from {c467.b_cl4_without_cov_A:.2f} to {c467.b_cl4_with_cov_A:.2f} A with "
                f"peptide_bound (mean signed coordinate {c467.chain_mean_s_by_cluster_level_A}; chains "
                f"{c467.chain_counts_by_cluster}). These covariates are unevenly spread over the clusters, so the "
                "adjusted coefficients rest on the few chains whose value is rare in their cluster.")
    L[-1] += txt
    L += ["", "## Inputs (rule 1.10)", "",
          "- Y10 card out results/ext/cluster_displacement_profile.csv (spike rule applied to it).",
          f"- Core v3 C-alpha: SV00 gate fit data/sv/gate/aligned_coordinates.npz = S03 results/alignment/aligned_coordinates.npz "
          f"(max |diff| {n['gate_vs_s03_coords']['max_abs_diff_A']:.1g} A, same gaps {n['gate_vs_s03_coords']['same_gaps']}); "
          f"Y10 reproduced with SV05's y10_profile: max |diff| {n['y10_reproduction_max_abs_diff']:.1g}, same NaN "
          f"{n['y10_reproduction_same_nan']}, medoids {n['y10_medoids'][0]} / {n['y10_medoids'][1]}.",
          "- SV00 observed_mask.npy, chain_table.csv (cluster, peptide_bound, ligand_chemotype, dimer_partner_source, "
          f"offsets), chain_offsets.csv; residue identity data/sv/lbd_extract/residue_mapping.csv "
          f"({n['residue_mapping_rows_vs_mask']['mapped_and_masked']} of {n['residue_mapping_rows_vs_mask']['mask_observed']} "
          "observed C-alpha mapped).",
          "- SV02 sse_per_position.csv, sse_codes_per_chain.csv, helix_boundaries.csv; SV01b B z formula; SV01b module "
          f"lattice_contacts.py ({lat['residues_found']} spike residues found = {lat['residues_ca_observed_mask']} "
          "observed in the SV00 mask, 0 mismatches).", ""]

    L += ["## Spikes (spikes.csv, verdict.csv)", "",
          "d = cl7-cl4 C-alpha distance (A). Entry = entry-level means. -3 = after greedy removal of 3 entries. "
          "Medoid restr. = complete-H12 medoids re-picked among chains observing the position. SSE = SV02 consensus "
          "helix cl7 / cl4 and H share among observing chains. B z = median per cluster. Lattice = share of observing "
          "chains with a symmetry-copy contact; partner = ASU copy of another chain.", "",
          "| pos | exceeds | obs cl7 / cl4 | d mean Y10 | d entry | d -3 (chain / entry) | d medoid Y10 / restr. | SSE cl7 / cl4 | B z cl7 / cl4 | lattice cl7 / cl4 | partner cl7 / cl4 | verdict |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in sp.itertuples():
        v = vd[vd.position == r.position].iloc[0]
        L.append(f"| {r.position} | {r.exceeds} | {r.obs_cl7:.0%} ({r.n_obs_cl7}) / {r.obs_cl4:.0%} ({r.n_obs_cl4}) | "
                 f"{f2(r.d_mean_y10_A)} | {f2(r.d_mean_entry_A)} | {f2(r.d_mean_minus_3_entries_A)} / "
                 f"{f2(r.d_mean_entry_minus_3_entries_A)} | {f2(r.d_medoid_y10_A)} / {f2(r.d_medoid_restricted_A)} | "
                 f"{r.consensus_helix_cl7} {r.frac_H_cl7:.2f} / {r.consensus_helix_cl4} {r.frac_H_cl4:.2f} | "
                 f"{f2(r.b_z_median_cl7, 1)} / {f2(r.b_z_median_cl4, 1)} | {r.frac_lattice_sym_cl7:.2f} / "
                 f"{r.frac_lattice_sym_cl4:.2f} | {r.frac_partner_asu_cl7:.2f} / {r.frac_partner_asu_cl4:.2f} | "
                 f"{v.verdict} |")
    L.append("")
    L.append(f"- Restricted recomputation: cluster means equal Y10 (max |diff| {rv['max_abs_diff_d_mean_A']:.1g} A). "
             f"Re-picked medoids differ from Y10 at {rv['n_positions_d_medoid_changed']} position where Y10 has a value "
             f"and give a value at {rv['n_positions_d_medoid_defined_only_restricted']} positions where the Y10 medoids "
             f"do not observe the residue; there they exceed 1.0 A at "
             f"{', '.join(map(str, rv['positions_restricted_medoid_gt_1_outside_528_550_where_y10_nan']))} "
             "(not spikes: the step-1 rule uses the Y10 values).")
    L.append(f"- Entry-level cluster means exceed 1.0 A outside 528-550 at "
             f"{', '.join(map(str, rv['positions_d_mean_entry_gt_1_outside_528_550']))} (max "
             f"{rv['max_d_mean_entry_outside_528_550_A']:.2f} A); 465 exceeds at entry level only "
             f"({f2(prof.set_index('position').loc[465].d_mean_entry_A)} A vs Y10 "
             f"{f2(prof.set_index('position').loc[465].d_mean_y10_A)} A)." if 465 in
             rv['positions_d_mean_entry_gt_1_outside_528_550'] and 465 not in sp.position.tolist() else
             f"- Entry-level cluster means exceed 1.0 A outside 528-550 at "
             f"{', '.join(map(str, rv['positions_d_mean_entry_gt_1_outside_528_550']))}.")
    t1 = top[top["rank"] == 1].set_index("position")
    L.append("- Largest single-entry influence per spike (top_entries.csv; 5 per spike with offsets, residue identity, "
             "SSE code, B z, peptide, ligand): " + "; ".join(
                 f"{p} {t1.loc[p].entry_id} {t1.loc[p].influence_A:+.2f} A" for p in sp.position) + ".")
    tc = []
    for r in top.itertuples():
        for ck, off, oc, src in zip(r.chains.split(";"), str(r.uniprot_author_offset).split(";"),
                                    str(r.offset_constant).split(";"), str(r.offset_source).split(";")):
            if float(off) != 0 or oc != "True" or src != "sifts_pdbe_api":
                tc.append(f"{ck} at {r.position} (offset {float(off):g}, constant {oc}, {src})")
    stat = sorted(set(";".join(top.residue_status.astype(str)).split(";")))
    L.append(f"- Numbering of the top entries (top_entries.csv): residue status at the spike position {stat} for every "
             f"listed chain; chains with a non-zero, non-constant or non-SIFTS offset among them: "
             f"{'; '.join(sorted(set(tc))) if tc else 'none'}. A constant non-zero offset (e.g. -1000 for author "
             "numbering 1305-1550) is a SIFTS-mapped constant, not a flag. Flagged chains overall (offset not constant "
             f"or not SIFTS): {', '.join(n['chains_offset_flag'])}.")
    L.append("")

    L += ["## Genuine spikes against peptide, ligand chemotype and dimer source (covariate_tests.csv)", "",
          "q = BH over all covariate tests (chain: OLS, entry-clustered Wald chi2; entry: entry x cluster units, nested "
          "F). Coef = cl4-minus-cl7 difference of the signed coordinate along the displacement (A) without -> with the "
          "covariate.", "",
          "| pos | peptide_bound q chain / entry | coef | ligand_chemotype q chain / entry | coef | dimer source q chain / entry | coef | associated (q < 0.05 both) |",
          "|---|---|---|---|---|---|---|---|"]
    for p in gen:
        g = cov[cov.position == p]
        cells = []
        for c in ("peptide_bound", "ligand_chemotype", "dimer_partner_source"):
            ch = g[(g.covariate == c) & (g.level == "chain")].iloc[0]
            en = g[(g.covariate == c) & (g.level == "entry")].iloc[0]
            cells.append(f"{ch.q_bh:.2g} / {en.q_bh:.2g} | {ch.b_cl4_without_cov_A:.2f} -> {ch.b_cl4_with_cov_A:.2f}")
        a = vd[vd.position == p].covariates_associated.iloc[0]
        a = ", ".join(x.split(" (")[0] for x in a.split("), ")) if not str(a).startswith("none") else "none"
        L.append(f"| {p} | " + " | ".join(cells) + f" | {a} |")
    L.append("")
    L.append("")
    L += ["## Figure", "",
          "- figures/sv/SV07/SV07_spike_profile.png/.svg: (a) Y10 profile, all chains and chains observing the "
          "position in both clusters (means coincide; medoids re-picked), 1.0 A threshold, spikes labelled with their "
          "verdict; (b) share of cl7 / cl4 chains observing each position. Panels in figures/panels/sv/SV07/.", ""]
    L += ["## Files", ""]
    files = [OUT / f for f in ("spikes.csv", "verdict.csv", "top_entries.csv", "restricted_profile.csv",
                               "covariate_tests.csv", "lattice_contacts_spikes_per_residue.csv",
                               "lattice_contacts_spikes_per_chain.csv", "lattice_checks.json", "numbers.json",
                               "env.txt")] + FIG + [ROOT / "scripts/sv" / s for s in SCRIPTS]
    for f in files:
        L.append(f"- {f.relative_to(ROOT)} ({size(f)})")
    (OUT / "SUMMARY.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
