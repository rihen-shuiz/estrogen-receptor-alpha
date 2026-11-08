"""SV04 (RUN_SV1028): write results/sv/SV04/env.txt and SUMMARY.md (rule 1.7) from the saved SV04 outputs only.
Run after SV04_midpoints.py, SV04_cluster.py and SV04_figures.py:
.venv/bin/python -B -I scripts/sv/SV04_summary.py
"""
import json
import platform
import sys
from datetime import datetime
from importlib import metadata
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV04"
HEL = ["H1", "H3", "H4/5", "H7", "H8", "H9", "H10/11", "H12"]


def size(path):
    n = path.stat().st_size
    return f"{n / 1e6:.2f} MB" if n >= 1e6 else f"{n / 1e3:.1f} kB"


def env():
    pk = ["numpy", "pandas", "scipy", "scikit-learn", "matplotlib"]
    lines = [f"# SV04 environment, {datetime.now().astimezone().strftime('%Y-%m-%dT%H:%M%z')}",
             f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}"]
    lines += [f"{p} {metadata.version(p)}" for p in pk]
    lines += ["secondary structure: SV02 outputs (pydssp 0.9.1 fallback for mkdssp, run by SV02; not re-run here)",
              "frame and axis code: SV03 (scripts/sv/common/core_frame.py, helix_geometry.py, h12_map.py; "
              "data/sv/SV03/core_v3_transforms.csv)",
              "nothing installed, upgraded or downgraded by SV04", "network: none"]
    (OUT / "env.txt").write_text("\n".join(lines) + "\n")


def f(x, d=2):
    return "NA" if pd.isna(x) else f"{x:.{d}f}"


def main():
    env()
    nm = json.loads((OUT / "numbers_midpoints.json").read_text())
    nc = json.loads((OUT / "numbers_cluster.json").read_text())
    ve = pd.read_csv(OUT / "validation_error.csv")
    ari = pd.read_csv(OUT / "ari.csv")
    ds = pd.read_csv(OUT / "distance_summary.csv")
    pa = pd.read_csv(OUT / "partial_h12_assignment.csv", keep_default_na=False, na_values=[""])
    mp = pd.read_csv(OUT / "midpoints.csv", keep_default_na=False, na_values=[""])
    ch = nm["complete_helices"]
    p = nc["partial"]
    off = nc["sv04_minus_sv03_map_x_by_status_group"]
    xc = nm["sv03_crosscheck"]
    ev = ve[ve.end == "either"]

    def a(dist, level, K):
        return ari[(ari.distance == dist) & (ari.level == level) & (ari.K == K)].iloc[0]

    L = []
    L.append(
        "Not done / substituted: (1) 'helical' comes from SV02's pydssp fallback (mkdssp not installed; 'H' = any "
        "helix). (2) Step 1 'observed C-alpha' read as SV03's 'observed helical C-alpha' (step 1: 'fit the axis as in "
        "SV03'): a residue observed but not pydssp-H counts as missing, so a helix is 'complete' (validation) only "
        "when every position of its consensus range is observed and helical, and 'observed length' = helical "
        "C-alpha used. (3) H12 = 538-547 (primary window, run file section 2), not the SV02 proposal 538-544: only "
        f"{ch['H12']['n']} chains have all of 538-547 helical (cl7 {ch['H12']['cl7']}, cl4 {ch['H12']['cl4']}), "
        "so the H12 error curve rests on them, and most chains' H12 'missing' residues are 546-547 observed but not "
        "helical, not unobserved. (4) [DEFAULT] error per observed length = median over N- and C-end deletions "
        "pooled (per-end curves are in validation_error.csv). (5) 'Report the axis angle separately' = RMS axis "
        "angle over the shared helices, clustered the same way with its own ARI rows, plus per-helix "
        "distance_summary.csv. (6) K = 2: linkage recomputed on the cl7 + cl4 sub-matrix; K = 7: all chains cut "
        "into 7 clusters vs the S08 1.5 A labels; entry level = the K07 chain per entry (results/entry_level_"
        "clusters.csv), sub-matrix re-linked; ARI is descriptive (no test, so no BH). (7) Step 5 'where they land' "
        "made operational three ways, fixed before the run: landing_h12 (primary; nearest of the cl7 / cl4 "
        "reference H12 midpoints = coordinate-wise median over complete-H12 chains), landing_k2 and landing_k7 "
        "(SV04 cluster at K = 2 / 7, labelled by the majority S08 group of its complete-H12 members); minor-cluster "
        "chains cannot agree with a cl7 / cl4 label, so their agreement is n/a. (8) Figure panel (b) is not in the "
        "step text: the SV04 H12 midpoint (axis point of 542.5) and the SV03 H12 centroid differ along the H12 axis, "
        f"which moves x by a cluster-dependent amount (complete cl4 median {off['complete|cl4']['median']:+.2f} A, "
        f"cl7 {off['complete|cl7']['median']:+.2f} A; y identical), so panel (a) as specified (complete chains at "
        "SV03 positions, partial at SV04 estimates) carries that method offset between filled and hollow markers; "
        "(b) draws every chain with the SV04 estimator. (9) Supporting files not named in the run file: "
        "numbers_midpoints.json, numbers_cluster.json, cluster_labels.csv, cluster_crosstab.csv, "
        "distance_summary.csv, h12_map_coords.csv, partial_h12_entry_level.csv, data/sv/SV04/, "
        "scripts/sv/common/helix_midpoint.py.")
    L.append("")
    L.append("# SV04 - Helix midpoints for incomplete helices (Ye)")
    L.append("")
    k2c, k7c, k2e, k7e = a("midpoint_rms", "chain", 2), a("midpoint_rms", "chain", 7), a("midpoint_rms", "entry", 2), \
        a("midpoint_rms", "entry", 7)
    h12d = ds[(ds.helix == "H12") & (ds.measure == "midpoint_distance_A")].iloc[0]
    oth = ds[(~ds.helix.isin(["H12", "RMS over shared helices"])) & (ds.measure == "midpoint_distance_A")]
    ag = p["agree_by_h12_eligible"]
    ok_e = ag["agree_h12|cl4|h12_eligible=True"][0] + ag["agree_h12|cl7|h12_eligible=True"][0]
    n_e = ag["agree_h12|cl4|h12_eligible=True"][1] + ag["agree_h12|cl7|h12_eligible=True"][1]
    k2_ne = [ag["agree_k2|cl4|h12_eligible=False"], ag["agree_k2|cl7|h12_eligible=False"]]
    dis = pa[(pa.group != "minor") & (pa.agree_k2.astype(str) == "False")]
    L.append(
        f"**Verdict:** clustering chains on helix midpoints reproduces S08: ARI {k2c.ari:.3f} at K = 2 (cl7 + cl4, "
        f"{int(k2c.n)} chains) and {k7c.ari:.3f} at K = 7 (all {int(k7c.n)} chains); entry level {k2e.ari:.3f} / "
        f"{k7e.ari:.3f}. All {nc['k2_chain_mismatch']['n']} chain-level K = 2 mismatches "
        f"({', '.join(f'{g} {v}' for g, v in nc['k2_chain_mismatch']['by_group'].items())}) are chains without a usable "
        f"H12 midpoint ({nc['k2_chain_mismatch']['n_without_usable_h12']} of "
        f"{nc['n_without_usable_h12_by_group'].get('cl7', 0) + nc['n_without_usable_h12_by_group'].get('cl4', 0)} "
        "cl7 / cl4 chains without one, clustered on the other 7 helices). The midpoint distance is dominated by H12 (cl7 vs cl4 H12 midpoints median "
        f"{h12d.cl7_vs_cl4_median:.1f} A apart; every other helix <= {oth.cl7_vs_cl4_median.max():.2f} A), so the "
        "agreement mostly restates the H12 split. Partial-H12 chains: "
        f"{p['h12_eligible']} of {p['n']} have an H12 midpoint usable by the [DEFAULT] rule; every cl7 / cl4 one "
        f"lands on its own S08 side by H12 position ({ok_e} of {n_e}). Partial chains without a usable H12 "
        f"midpoint are placed by the other helices only and agree less (K = 2: cl4 {k2_ne[0][0]} of {k2_ne[0][1]}, "
        f"cl7 {k2_ne[1][0]} of {k2_ne[1][1]}; disagreeing: " + ", ".join(
            f"{r.chain_key} ({r.group}, {r.n_obs_window_538_547} of 10 H12 C-alpha observed -> {r.landing_k2}-like)"
            for r in dis.itertuples()) + "). "
        "Agreement with S08 is not independent evidence: the S08 RMSD already used the same observed H12 C-alpha.")
    L.append("")

    # step 3
    L.append("## Steps 1-3 - midpoints and truncation validation (validation_error.csv, midpoints.csv)")
    L.append("")
    L.append(f"- Frame / axes: core v3 (SV03 frozen frame, load_ca_core); SV04 fit sets equal SV03's for "
             f"{xc['n_fit_ca_equal']} of {xc['n_rows']} chain x helix rows, axis status equal {xc['n_axis_status_equal']}"
             f"/{xc['n_rows']}, max axis angle difference {xc['max_axis_angle_diff_deg']:.1e} deg, every midpoint lies "
             f"on the SV03 axis (max offset {xc['max_midpoint_offset_from_sv03_axis_A']:.1e} A). H12 midpoint minus SV03 "
             f"H12 centroid along the axis: median {xc['midpoint_minus_sv03_centroid_along_axis_A']['H12']['median']:+.2f} A "
             f"(max |{xc['midpoint_minus_sv03_centroid_along_axis_A']['H12']['max_abs']:.2f}| A); other helices median "
             f"within +-{max(abs(v['median']) for k, v in xc['midpoint_minus_sv03_centroid_along_axis_A'].items() if k != 'H12'):.2f} A.")
    L.append("- Rise per residue on complete helices (median): " + ", ".join(
        f"{h} {nm['rise_per_residue_A_complete'][h]:.2f}" for h in HEL) + " A.")
    L.append("- Complete helices used for validation (all positions observed and helical): " + ", ".join(
        f"{h} {ch[h]['n']} (cl7 {ch[h]['cl7']}, cl4 {ch[h]['cl4']}, minor {ch[h]['minor']})" for h in HEL) + ".")
    L.append("")
    L.append("Median midpoint error (A) / median axis angle error (deg) after deleting k residues from either end "
             "(N and C pooled; n = truncations estimated); * = fails the [DEFAULT] (median > 1.0 A); - = fewer than "
             "6 C-alpha left:")
    L.append("")
    L.append("| helix | range | L | " + " | ".join(f"k={k}" for k in range(1, 7)) + " | eligible k |")
    L.append("|---|---|---|" + "---|" * 7)
    for h in HEL:
        x = ev[ev.helix == h].set_index("k_removed")
        cells = []
        for k in range(1, 7):
            r = x.loc[k]
            if r.n_estimated == 0:
                cells.append("-")
            else:
                star = "" if str(r.eligible_observed_length) == "True" else " *"
                cells.append(f"{r.midpoint_err_median_A:.2f} / {r.axis_err_median_deg:.1f}{star} (n={int(r.n_estimated)})")
        L.append(f"| {h} | {int(x.range_start.iloc[0])}-{int(x.range_end.iloc[0])} | {int(x.consensus_length.iloc[0])} | "
                 + " | ".join(cells) + f" | {','.join(map(str, nm['eligible_k'][h]))} |")
    L.append("")
    eb = nm["eligible_by_helix"]
    L.append("- [DEFAULT] applied: eligible chain x helix " + ", ".join(
        f"{h} {eb[h]['eligible']}/{nm['n_chains']}" for h in HEL) + ". Every chain x helix with an axis falls on a "
        "validated observed length that passes (H4/5 fails at k >= 5, but no chain misses more than 4 H4/5 residues; "
        "H12 k distribution " + ", ".join(f"k={k}: {v}" for k, v in nm['k_distribution']['H12'].items()) +
        "); the only exclusions are the "
        f"{nm['n_chains'] - eb['H12']['eligible']} H12 helices with < 6 observed helical C-alpha (no axis, as in SV03). "
        "Eligible helices per chain: " + ", ".join(f"{k} helices {v} chains" for k, v in
                                                   nm['eligible_helices_per_chain'].items()) + ".")
    hb = nm["h12_by_status"]
    L.append("- H12 by SV01a status: " + "; ".join(
        f"{st} {v['n']} chains, eligible {v['eligible']} (cl7 {v['eligible_by_group']['cl7']}/{v['n_by_group']['cl7']}, "
        f"cl4 {v['eligible_by_group']['cl4']}/{v['n_by_group']['cl4']}, minor {v['eligible_by_group']['minor']}/"
        f"{v['n_by_group']['minor']})" for st, v in hb.items()) + ".")
    L.append("")

    # step 4
    L.append("## Step 4 - clustering on midpoints (ari.csv, cluster_crosstab.csv, distance_summary.csv)")
    L.append("")
    L.append(f"- Distance = RMS of midpoint distances over shared eligible helices; pairs share "
             + ", ".join(f"{k} helices: {v}" for k, v in nc['n_shared_helices'].items())
             + f" pairs; no chain dropped ({len(nc['dropped_no_shared_helix'])}). Average linkage.")
    L.append("")
    L.append("| distance | level | K | chains | n | SV04 cluster sizes | ARI vs S08 |")
    L.append("|---|---|---|---|---|---|---|")
    for r in ari.itertuples():
        L.append(f"| {r.distance} | {r.level} | {r.K} | {r.chains} | {r.n} | {r.sv04_cluster_sizes} | {r.ari:.3f} |")
    L.append("")
    L.append("Per-helix pair distances (chain level, eligible helices; median [90th percentile]):")
    L.append("")
    L.append("| helix | midpoint within cl7 (A) | within cl4 | cl7 vs cl4 | axis angle within cl7 (deg) | within cl4 | cl7 vs cl4 |")
    L.append("|---|---|---|---|---|---|---|")
    for h in HEL + ["RMS over shared helices"]:
        m_ = ds[(ds.helix == h) & (ds.measure == "midpoint_distance_A")].iloc[0]
        a_ = ds[(ds.helix == h) & (ds.measure == "axis_angle_deg")].iloc[0]
        L.append(f"| {h} | " + " | ".join(f"{m_[f'{c}_median']:.2f} [{m_[f'{c}_q90']:.2f}]"
                                          for c in ("within_cl7", "within_cl4", "cl7_vs_cl4")) + " | "
                 + " | ".join(f"{a_[f'{c}_median']:.1f} [{a_[f'{c}_q90']:.1f}]"
                              for c in ("within_cl7", "within_cl4", "cl7_vs_cl4")) + " |")
    L.append("")

    # step 5
    L.append("## Step 5 - partial-H12 chains (partial_h12_assignment.csv, partial_h12_entry_level.csv)")
    L.append("")
    rf = nc["reference_h12_midpoints"]
    L.append(f"- Partial H12 (SV01a: 1-9 of the 10 C-alpha 538-547 observed): {p['n']} chains "
             f"({', '.join(f'{g} {n}' for g, n in p['by_group'].items())}); usable H12 midpoint {p['h12_eligible']}; "
             f"placed on the map (H12, H3, H10/11 usable) {p['placed_on_map']} "
             f"({', '.join(f'{g} {n}' for g, n in p['placed_by_group'].items())}).")
    L.append(f"- Reference H12 midpoints (complete-H12 chains with a usable H12): cl7 n = {nc['n_reference']['cl7']}, "
             f"cl4 n = {nc['n_reference']['cl4']}; {nc['reference_distance_cl7_cl4_A']:.1f} A apart. "
             "Margin d(cl4 ref) - d(cl7 ref) for placed partial chains: " + "; ".join(
                 f"{g} median {v['median']:+.1f} A (range {v['min']:+.1f} to {v['max']:+.1f})"
                 for g, v in p['margin_d4_minus_d7_A'].items()) + ".")
    L.append(f"- SV04 cluster labels (majority S08 group of complete-H12 members): K = 2 {nc['cluster_label_k2']}; "
             f"K = 7 {nc['cluster_label_k7']}.")
    L.append("")
    L.append("| S08 group | partial chains | landing_h12 (H12 midpoint) | agree | landing_k2 | agree | landing_k7 | agree |")
    L.append("|---|---|---|---|---|---|---|---|")
    for g in ("cl7", "cl4", "minor"):
        if g not in p["landing_h12"]:
            continue

        def cnt(col):
            d = {k: v for k, v in p[col][g].items() if k != "n"}
            return ", ".join(f"{k} {v}" for k, v in d.items())

        def agr(col):
            a2 = p[col][g]
            return "n/a" if a2[1] == 0 else f"{a2[0]}/{a2[1]}"
        L.append(f"| {g} | {p['landing_h12'][g]['n']} | {cnt('landing_h12')} | {agr('agree_h12')} | "
                 f"{cnt('landing_k2')} | {agr('agree_k2')} | {cnt('landing_k7')} | {agr('agree_k7')} |")
    L.append("")
    L.append("- Agreement split by whether the chain has a usable H12 midpoint (agree / chains): " + "; ".join(
        f"{k.replace('|h12_eligible=', ', H12 usable ')} {v[0]}/{v[1]}" for k, v in ag.items() if v[1] and
        "minor" not in k) + ".")
    mi = pa[(pa.group == "minor") & (pa.h12_eligible.astype(str) == "True")]
    L.append("- Minor-cluster partial chains with an H12 midpoint (distance to the cl7 / cl4 reference, A): " + "; ".join(
        f"{r.chain_key} (S08 cl{r.s08_cluster}) {r.d_ref_cl7_A:.1f} / {r.d_ref_cl4_A:.1f} -> {r.landing_h12}-like"
        for r in mi.itertuples()) + ".")
    el = nc["entry_level"]
    L.append("- Entry level (majority landing vs majority S08 group of an entry's landed partial chains): " + "; ".join(
        f"{c}: {v['agree']} agree, {v['disagree']} disagree, {v['tie']} tie, {v['n_a_minor']} n/a (minor) of "
        f"{v['n_entries']} entries" for c, v in el.items()) + ".")
    sm = p["sv04_minus_sv03_map"]
    L.append(f"- Map positions of the {sm['n']} placed partial chains, SV04 estimate minus SV03 position: x median "
             f"{sm['x_median']:+.2f} A (max |{sm['x_max_abs']:.2f}|), y max |{sm['y_max_abs']:.1e}| deg; every placed "
             f"partial chain already had an SV03 point ({sm['placed_without_sv03_point']} new). Same difference for "
             "complete-H12 chains: " + ", ".join(
                 f"{k.split('|')[1]} {v['median']:+.2f} A (n={v['n']})" for k, v in off.items() if k.startswith("complete"))
             + "; the x offset is a method difference (midpoint at residue 542.5 vs centroid of the observed helical "
             "windows), not a position of the partial chains.")
    L.append("")
    L.append("## Figure")
    L.append("")
    lay = nc["map_layers"]
    L.append(f"- figures/sv/SV04/h12_map_partial.png/.svg: SV03 H12 map axes (x = H12 point to H11 axis, y = angle "
             f"H12-H3; colour = S08 cluster, marker = ligand_chemotype). (a) as specified: "
             f"{lay['filled (SV03 position)']} complete-H12 chains filled at SV03 positions, "
             f"{lay['hollow (SV04 estimate)']} partial-H12 chains hollow at SV04 estimates; (b) the same chains, all at "
             f"SV04 estimates. Not drawn: {lay['not drawn']} chains (not_modelled H12, or no usable H12 / H3 / H10/11).")
    L.append("")
    L.append("## Files")
    L.append("")
    files = sorted(OUT.glob("*")) + sorted((ROOT / "data/sv/SV04").glob("*")) + \
        sorted((ROOT / "figures/sv/SV04").glob("*")) + sorted((ROOT / "figures/panels/sv/SV04").glob("*")) + \
        [ROOT / f"scripts/sv/{n}" for n in ("SV04_midpoints.py", "SV04_cluster.py", "SV04_figures.py",
                                            "SV04_summary.py", "common/helix_midpoint.py")]
    for fp in files:
        if fp.name == "SUMMARY.md":
            continue
        L.append(f"- {fp.relative_to(ROOT)} ({size(fp)})")
    L.append("- results/sv/SV04/SUMMARY.md (this file)")
    (OUT / "SUMMARY.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
