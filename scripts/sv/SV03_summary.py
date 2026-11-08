"""SV03 (RUN_SV1028): write results/sv/SV03/env.txt and SUMMARY.md (rule 1.7) from the saved SV03 outputs only.
Run after SV03_geometry.py, SV03_tests.py and SV03_figures.py:
.venv/bin/python -B -I scripts/sv/SV03_summary.py
"""
import json
import platform
import sys
from datetime import datetime
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV03"
SMALL_ENTRIES = 5   # flag (not a rule): a DEFAULT hit whose smaller group has <= 5 entries
LABEL = {"ulen_H1": "length H1", "ulen_H3": "length H3", "ulen_H45": "length H4/5", "ulen_H7": "length H7",
         "ulen_H8": "length H8", "ulen_H9": "length H9", "ulen_H1011": "length H10/11", "ulen_H12": "length H12",
         "angle_H12_H3_deg": "angle H12-H3", "angle_H12_H11_deg": "angle H12-H11", "angle_H12_H45_deg": "angle H12-H4/5",
         "angle_H11_H3_deg": "angle H11-H3", "dist_H12c_H3axis_A": "H12 centroid-H3 axis",
         "dist_H12c_H11axis_A": "H12 centroid-H11 axis", "kink_H11_deg": "H11 kink", "hinge_H11C_H12_deg": "H11-H12 hinge"}
UNIT = {"ulen": "res", "angle": "deg", "dist": "A", "kink": "deg", "hinge": "deg"}
COV = {"ligand_heavy_atoms": "heavy atoms", "ligand_basic_amine": "basic amine",
       "ligand_carboxylic_acid": "carboxylic acid", "peptide_bound": "peptide", "mut_Y537S": "Y537S",
       "mut_L536S": "L536S", "mut_D538G": "D538G", "mut_other": "other mutation"}


def size(path):
    n = path.stat().st_size
    return f"{n / 1e6:.2f} MB" if n >= 1e6 else f"{n / 1e3:.1f} kB"


def fq(q):
    return "n.t." if pd.isna(q) else (f"{q:.1e}" if q < 1e-3 else f"{q:.3f}")


def fv(x, d=1):
    return "NA" if pd.isna(x) else f"{x:.{d}f}"


def unit(d):
    return UNIT[d.split("_")[0]]


def env():
    pk = ["numpy", "pandas", "scipy", "statsmodels", "matplotlib"]
    lines = [f"# SV03 environment, {datetime.now().astimezone().strftime('%Y-%m-%dT%H:%M%z')}",
             f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}"]
    lines += [f"{p} {metadata.version(p)}" for p in pk]
    lines += ["secondary structure: SV02 outputs (pydssp 0.9.1 fallback for mkdssp, run by SV02; not re-run here)",
              "nothing installed, upgraded or downgraded by SV03", "network: none"]
    (OUT / "env.txt").write_text("\n".join(lines) + "\n")


def main():
    env()
    fc = json.loads((OUT / "frame_check.json").read_text())
    ng = json.loads((OUT / "numbers_geometry.json").read_text())
    nt = json.loads((OUT / "numbers_tests.json").read_text())
    hl = pd.read_csv(OUT / "helix_length_summary.csv")
    aw = pd.read_csv(OUT / "cl7_vs_cl4_wide.csv")
    sc = pd.read_csv(OUT / "sensitivity_check.csv")
    t = pd.read_csv(OUT / "tests.csv")
    hs = pd.read_csv(OUT / "h12_window_sensitivity.csv")
    w = pd.read_csv(OUT / "geometry_per_chain.csv", keep_default_na=False, na_values=[""])
    am = pd.read_csv(OUT / "axis_method_check.csv")

    gaps = hl[(hl.start_gap + hl.end_gap) > 0]
    gap_txt = ", ".join(f"{r.helix} {int(r.start_gap)} start / {int(r.end_gap)} end" for r in gaps.itertuples()
                        if r.helix != "H12_w544")
    hits = sc[sc.meets_default].copy()
    # entries0 / entries1 in sensitivity_check.csv are the chain-level rows' entry counts per group
    hits["min_entries"] = np.where(hits.covariate == "ligand_heavy_atoms", hits.entries1,
                                   hits[["entries0", "entries1"]].min(axis=1))
    hits["small"] = hits.min_entries <= SMALL_ENTRIES
    nb = nt["n_tests"]

    L = []
    L.append(
        "Not done / substituted: (1) 'helical' and helix ends come from SV02's pydssp fallback (mkdssp not installed; "
        "pydssp 'H' = any helix, 3-, 4- and 5-turn together), not from DSSP alpha-H. (2) H12 range = 538-547, the "
        "primary H12 window of run file section 2, not the SV02 data-derived consensus H12 538-544 (a proposal); the "
        "H12 descriptors with 538-544 are kept as a descriptive sensitivity copy (geometry_per_chain.csv *_w544, "
        "h12_window_sensitivity.csv), not tested. (3) 'H11' = SV02 consensus helix H10/11 497-526 (SV02 keeps "
        "H10/H11 combined); its halves split that range at the midpoint (H11N 497-511, H11C 512-526); 'H11-H12 hinge "
        "angle' = angle H11C-H12 (the half next to the hinge), distinct from angle H12-H11 (whole H11) - an "
        "interpretation, fixed before the run. (4) A third terminus state the step does not name: 'gap' = next "
        "residue not helical but the second one helical (1-residue break before another H run): neither censored "
        f"nor a 'real end', so not an uncensored length ({gap_txt}). (5) Step-5 [DEFAULT] for the one continuous "
        "covariate (ligand heavy atoms): Spearman |rho| >= 0.3 used in place of |rank-biserial r|. (6) 'mutations' "
        "coded as in SV01a: Y537S, L536S, D538G, other construct mutation (each present vs absent). (7) Chain-level "
        "rank test = GEE on average ranks with entry clusters (the SV11 convention); entry level = Mann-Whitney U / "
        "Spearman on entry medians; BH within four families ((a) chain, (a) entry, (b) chain, (b) entry). Tests with "
        "< 2 entries in a group were not run; a group of 2-5 entries passes this rule, and cluster-robust GEE errors "
        f"are optimistic with so few clusters: {int(hits.small.sum())} of the {len(hits)} [DEFAULT] hits rest on such a "
        "group (flagged below). (8) Frame = SV00 gate transforms (S03 code path, equal to results/alignment/ within "
        "3.6e-14 A per SV00); frozen copy data/sv/SV03/core_v3_transforms.csv. (9) Supporting files not named in the "
        "run file: helix_per_chain.csv, helix_length_summary.csv, h12_window_sensitivity.csv, axis_method_check.csv, "
        "sensitivity_check.csv, "
        "cl7_vs_cl4_wide.csv, frame_check.json, numbers_geometry.json, numbers_tests.json, data/sv/SV03/, "
        "scripts/sv/common/{core_frame,helix_geometry,h12_map}.py.")
    L += ["", "# SV03 - Helix geometry", ""]
    n_hit_cl = hits.cluster.value_counts().to_dict()
    n_hit_h12 = int(hits.involves_H12.sum())
    L.append(f"**Verdict ([DEFAULT] step 5):** {len(hits)} within-cluster comparisons meet |r| >= 0.3 and q < 0.05 at "
             f"both chain and entry level ({n_hit_cl.get('cl7', 0)} in cl7, {n_hit_cl.get('cl4', 0)} in cl4; "
             f"{n_hit_h12} use H12, {len(hits) - n_hit_h12} only other helices), so by the rule the geometric "
             f"descriptors are **more sensitive than RMSD** for those subgroups; {int(hits.small.sum())} of them rest on "
             f"a group of <= {SMALL_ENTRIES} entries. cl7 vs cl4 (step 4a): "
             f"{len(nt['a_separates_both_levels'])} of 16 descriptors differ at q < 0.05 at both levels, all 7 H12 "
             f"descriptors among them; not separated: {', '.join(LABEL[d] for d in nt['a_not_separating'])}.")
    a_ = aw.set_index("descriptor")
    bh = hits[hits.covariate != "ligand_heavy_atoms"].copy()
    bh["dmed"] = (bh.median1_chain - bh.median0_chain).abs()
    mag = {k: bh[bh.descriptor.str.startswith(k)].dmed for k in ("angle", "hinge", "kink", "dist", "ulen")}
    L.append("")
    L.append("Size of the within-cluster differences (binary hits, |median yes - median no| at chain level): angles "
             f"incl. kink and hinge {pd.concat([mag['angle'], mag['hinge'], mag['kink']]).min():.1f}-"
             f"{pd.concat([mag['angle'], mag['hinge'], mag['kink']]).max():.1f} deg, distances "
             f"{mag['dist'].min():.2f}-{mag['dist'].max():.2f} A, lengths {mag['ulen'].min():.0f}-{mag['ulen'].max():.0f} "
             "residues; for comparison the cl7 - cl4 chain medians differ by "
             f"{abs(a_.loc['angle_H12_H3_deg'].median1_chain - a_.loc['angle_H12_H3_deg'].median0_chain):.1f} deg "
             "(angle H12-H3) and "
             f"{abs(a_.loc['dist_H12c_H11axis_A'].median1_chain - a_.loc['dist_H12c_H11axis_A'].median0_chain):.1f} A "
             "(H12 centroid-H11 axis) (step 4a table).")
    L += ["", "## Frame (step inputs)", "",
          f"- Core v3 frame of the SV00 gate (reference 7NFB|B|1, {fc['n_chains']} chains, status {fc['status_counts']}); "
          f"convention `{fc['convention']}`; deposited backbone C-alpha transformed with the frozen transforms vs "
          f"data/sv/gate/aligned_coordinates.npz: max |d| = {fc['max_abs_diff_backbone_CA_vs_aligned_A']:.1e} A, "
          f"same observed mask {fc['observed_mask_identical']}.",
          "- Reuse (SV04, SV13): data/sv/SV03/core_v3_transforms.csv + scripts/sv/common/core_frame.py "
          "(load_transforms, to_core, load_ca_core); axis fit scripts/sv/common/helix_geometry.py (fit_axis); H12 map "
          "scripts/sv/common/h12_map.py (plot_h12_map, add_legends; marker column and hollow markers are parameters).",
          "- Run ends and censoring reproduce SV02 sse_per_chain.csv exactly for every helix with the same range "
          "(status, start, end, start/end censored: " + "; ".join(
              f"{k} {v['start_end_equal']}/{v['n_helix']}" for k, v in ng["sv02_crosscheck"].items()) + ")."]
    L += ["", "## Step 1 - helix lengths and censoring (helix_length_summary.csv)", "",
          "Length = first to last pydssp-H residue of the chain's runs overlapping the range; censored = a terminus "
          "within 2 residues of an unobserved residue or chain break; uncensored = both termini real ends.", "",
          "| helix | range | chains with a run | start cens. | end cens. | gap ends | censoring fraction (all / cl7 / cl4) "
          "| uncensored | median [IQR] uncensored length cl7 | cl4 |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in hl.itertuples():
        if r.helix == "H12_w544":
            continue
        L.append(f"| {r.helix} | {r.range} | {r.n_with_run} | {r.start_censored} | {r.end_censored} | "
                 f"{r.start_gap + r.end_gap} | {r.censoring_fraction:.3f} / {r.censoring_fraction_cl7:.3f} / "
                 f"{r.censoring_fraction_cl4:.3f} | {r.n_uncensored} | {fv(r.median_uncensored_length_cl7, 0)} "
                 f"[{fv(r.q1_cl7, 0)}-{fv(r.q3_cl7, 0)}] (n={r.n_uncensored_cl7}) | "
                 f"{fv(r.median_uncensored_length_cl4, 0)} [{fv(r.q1_cl4, 0)}-{fv(r.q3_cl4, 0)}] (n={r.n_uncensored_cl4}) |")
    h12 = hl[hl.helix == "H12"].iloc[0]
    L += ["", f"- H12 538-547: {int(h12.n_not_observed)} chains observe none of it, {int(h12.n_not_helical)} have no "
              f"H residue in it; {int(h12.n_with_run)} have a run."]
    L += ["", "## Step 2-3 - axes and descriptors (geometry_per_chain.csv, helix_per_chain.csv)", "",
          "- Axis [DEFAULT]: >= 6 observed helical C-alpha in the range, centroids of consecutive 4-residue windows, "
          "principal line, oriented N -> C; helix centroid = mean window centroid (core v3 frame, "
          "centroid_<helix>_x/y/z).",
          "- Axes per helix: " + "; ".join(f"{h} {v.get('ok', 0)}/{sum(v.values())}" for h, v in ng["axis_status"].items()
                                          if h != "H12_w544") +
          f". H12 without an axis: {788 - ng['n_h12_axis']} chains (< 6 observed helical C-alpha in 538-547).",
          "- Axis rule on ideal alpha-helices (axis_method_check.csv): angle error " + ", ".join(
              f"{int(r.n_ca)} CA {r.angle_error_deg:.1f} deg" for r in am.itertuples()) + "; H12 axes use "
          + "; ".join(f"{g}: " + ", ".join(f"{k} CA {v}" for k, v in d.items())
                      for g, d in ng["h12_axis_n_ca_by_group"].items())
          + " chains - cl4 H12 axes rest mostly on 7 C-alpha, cl7 on 8, so H12 angle differences of a few degrees "
          "inside a cluster are within the method's range for short H12 runs.",
          f"- H12 map points (chains with an H12 axis): {ng['n_h12_map_points']} "
          f"({', '.join(f'{g} {int(((w.group == g) & w.axis_H12).sum())}' for g in ('cl7', 'cl4', 'minor'))})."]
    L += ["", "## Step 4a - cl7 vs cl4 (tests.csv part a, cl7_vs_cl4_wide.csv)", "",
          f"Chain level {int(aw.n0_chain.max())} cl7 / {int(aw.n1_chain.max())} cl4 chains at most; entry level "
          f"entry groups by the SV10 rule ({nt['entry_groups']}). r = rank-biserial, positive = higher in cl4.", "",
          "| descriptor | uses H12 | median cl7 | median cl4 | r chain | q chain | r entry | q entry | q < 0.05 both |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in aw.itertuples():
        d = r.descriptor
        L.append(f"| {LABEL[d]} ({unit(d)}) | {'yes' if 'H12' in d else 'no'} | {fv(r.median0_chain, 2)} | "
                 f"{fv(r.median1_chain, 2)} | {r.effect_chain:+.2f} | {fq(r.q_chain)} | {r.effect_entry:+.2f} | "
                 f"{fq(r.q_entry)} | {'yes' if r.separates_both_levels else 'no'} |")
    L += ["", "## Step 4b-5 - within cl7 and within cl4 (tests.csv part b, sensitivity_check.csv)", "",
          f"- Comparisons: 16 descriptors x 8 covariates x 2 clusters = {nt['b_comparisons']}; tested at both levels "
          f"{nt['b_tested_both_levels']} (chain-level tests {nb['(b) chain']['sum']}/{nb['(b) chain']['size']}, "
          f"entry-level {nb['(b) entry']['sum']}/{nb['(b) entry']['size']}); not tested (< 2 chains or entries in a "
          "group, per descriptor): " + "; ".join(
              f"{cl} {COV[cv]} {n} descriptors" for (cl, cv), n in
              t[(t.part == "b") & (t.level == "chain") & ~t.tested.astype(bool)].groupby(["cluster", "covariate"])
              .size().items()) + ".",
          f"- q < 0.05 at both levels: {nt['b_q_lt_005_both']}; of these, |effect| >= 0.3 at both levels "
          f"([DEFAULT] met): {nt['b_meets_default']}; sign equal at both levels for "
          f"{int(hits.same_sign.sum())}/{len(hits)}.",
          "- Hits by covariate: " + ", ".join(f"{COV[k]} {v}" for k, v in hits.covariate.value_counts().items()) + ".",
          "", "| cluster | descriptor | covariate | chains (no / yes) | entries (no / yes) | median no / yes (chain) "
          "| r chain | q chain | r entry | q entry | small group |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in hits.sort_values(["cluster", "involves_H12", "descriptor"], ascending=[False, True, True]).itertuples():
        d = r.descriptor
        if r.covariate == "ligand_heavy_atoms":
            chains, ents, med = f"{int(r.n1_chain)} (rho)", f"{int(r.entries1)}", "-"
        else:
            chains = f"{int(r.n0_chain)} / {int(r.n1_chain)}"
            ents = f"{int(r.entries0)} / {int(r.entries1)}"
            med = f"{fv(r.median0_chain, 2)} / {fv(r.median1_chain, 2)} {unit(d)}"
        L.append(f"| {r.cluster} | {LABEL[d]} | {COV[r.covariate]} | {chains} | {ents} | {med} | "
                 f"{r.effect_chain:+.2f} | {fq(r.q_chain)} | {r.effect_entry:+.2f} | {fq(r.q_entry)} | "
                 f"{'yes' if r.small else ''} |")
    L += ["", "r = rank-biserial (positive = higher with the covariate), rho = Spearman for heavy atoms; small group = "
          f"<= {SMALL_ENTRIES} entries in the smaller group (a flag for reading, not part of the rule)."]
    L += ["", "## H12 window sensitivity (h12_window_sensitivity.csv, descriptive)", "",
          "| descriptor | chains | Spearman rho 538-547 vs 538-544 | median abs diff | max abs diff |", "|---|---|---|---|---|"]
    for r in hs.itertuples():
        L.append(f"| {r.descriptor} | {r.n_both} | {r.spearman_rho:.3f} | {r.median_abs_diff:.2f} | {r.max_abs_diff:.2f} |")
    L += ["", "Run ends do not depend on the window here: every H12 run overlapping 538-547 also overlaps 538-544, so "
          "lengths are identical; axes change because 545-547 add C-alpha."]
    L += ["", "## Figures", "",
          "- figures/sv/SV03/SV03_h12_map.png/.svg: x = H12 centroid to H11 axis, y = angle H12-H3, colour = S08 "
          "cluster, marker = ligand_chemotype, one point per chain with an H12 axis.",
          "- figures/sv/SV03/SV03_descriptors_by_cluster.png/.svg: the 16 tested descriptors by cluster group with the "
          "step-4a r and q per panel."]
    files = sorted([p for p in OUT.iterdir() if p.is_file() and p.name != "SUMMARY.md"]) + \
        sorted((ROOT / "data/sv/SV03").iterdir()) + sorted((ROOT / "figures/sv/SV03").iterdir()) + \
        sorted((ROOT / "figures/panels/sv/SV03").iterdir()) + \
        [ROOT / "scripts/sv" / f for f in ("SV03_geometry.py", "SV03_tests.py", "SV03_figures.py", "SV03_summary.py")] + \
        [ROOT / "scripts/sv/common" / f for f in ("core_frame.py", "helix_geometry.py", "h12_map.py")]
    L += ["", "## Files", ""] + [f"- {p.relative_to(ROOT)} ({size(p)})" for p in files] + \
         ["- results/sv/SV03/SUMMARY.md (this file)"]
    (OUT / "SUMMARY.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
