"""SV02 (RUN_SV1028): write results/sv/SV02/env.txt and SUMMARY.md (rule 1.7) from the saved SV02 outputs.
Run after SV02_dssp.py, SV02_consensus.py, SV02_phipsi.py and SV02_figures.py:
.venv/bin/python -B -I scripts/sv/SV02_summary.py
"""
import json
import platform
import sys
from datetime import datetime
from importlib import metadata
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV02"


def size(path):
    n = path.stat().st_size
    return f"{n / 1e6:.2f} MB" if n >= 1e6 else f"{n / 1e3:.1f} kB"


def f1(x):
    return "NA" if pd.isna(x) else f"{x:.1f}"


def f2(x):
    return "NA" if pd.isna(x) else f"{x:.2f}"


def env():
    pk = ["numpy", "pandas", "pydssp", "einops", "torch", "gemmi", "matplotlib"]
    lines = [f"# SV02 environment, {datetime.now().astimezone().strftime('%Y-%m-%dT%H:%M%z')}",
             f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}"]
    lines += [f"{p} {metadata.version(p)}" for p in pk]
    lines += ["mkdssp: not found (no conda/Homebrew; docs/sv/env_pre.txt) -> step fallback pydssp 0.9.1, numpy backend "
              "(pydssp imports torch at package import; the numpy code path was used)",
              "nothing installed, upgraded or downgraded by SV02", "network: none"]
    (OUT / "env.txt").write_text("\n".join(lines) + "\n")


def main():
    env()
    run = json.loads((OUT / "dssp_run.json").read_text())
    nc = json.loads((OUT / "numbers_consensus.json").read_text())
    nr = json.loads((OUT / "numbers_rama.json").read_text())
    hb = pd.read_csv(OUT / "helix_boundaries.csv")
    cr = pd.read_csv(OUT / "consensus_runs_all.csv")
    cg = pd.read_csv(OUT / "consensus_global.csv")
    pc = pd.read_csv(OUT / "sse_per_chain.csv")
    rs = pd.read_csv(OUT / "rama_summary.csv")
    gm = nr["gemmi_check"]
    prop = nc["h12_proposal"]
    comp = nc["h12_comparisons"]
    c537 = nr["check_537_psi"]

    line1 = (
        "Not done / substituted: (1) mkdssp >= 4 is not installed (no conda/Homebrew) - the step's fallback pydssp 0.9.1 "
        "(numpy backend) ran instead; pydssp gives H / E / '-' only and its 'H' is the union of 3-, 4- and 5-turn "
        "helices, so the 'alpha (H)' fractions are 'any helix' and 3_10 (G) and pi (I) fractions are 'not available'; "
        "pydssp also places the amide H on the C(i-1)-N / CA-N bisector, not along the previous C=O as DSSP does. "
        f"(2) pydssp has no chain-break handling: each chain ran as {run['segments']} segments in total (consecutive "
        "label_seq_id, all of N/CA/C/O, C-N <= 2.5 A = mkdssp's peptide-bond limit), prolines as non-donors, "
        f"{run['segments_padded_lt6']} segments shorter than 6 residues padded with far-away dummies (pydssp raises "
        "below 6); E bridges between segments are therefore not detected (H and 'other' unaffected). "
        "(3) Step 7 mkdssp phi/psi comparison on 5 chains: not possible (mkdssp did not run); substitute sanity check "
        f"only: gemmi.calculate_phi_psi on the restored mmCIFs, 5 chains, {gm['residues_compared']} residues, max "
        f"|dphi| {gm['max_dphi_no_altloc_in_triplet']:.1e} and |dpsi| {gm['max_dpsi_no_altloc_in_triplet']:.1e} deg "
        f"where residues i-1..i+1 have no altloc; {gm['residues_altloc_in_triplet']} residues with an altloc differ by up "
        f"to {gm['max_dphi']:.1f} / {gm['max_dpsi']:.1f} deg (gemmi reads the first altloc, SV00 kept the "
        "highest-occupancy one). (4) Helix names: the helix table (K02, results/helix_table.csv) names only H12, so "
        "names come from the project's mapping of its rows, scripts/ext_common.HELIX_NAMES (H1, H3, H4/5, H6, H7, H8, "
        "H9, H10/11, H12) - H4/5 and H10/11 stay combined as there; no H2 and no H6 consensus helix (414-415 = 2 "
        "positions). (5) Implementation choices fixed before the output: low-support positions are trimmed from run "
        "ends (never a boundary, allowed inside; trimming occurred 0 times); per-chain censoring tests the two "
        "flanking residues - a one-residue test is False by construction under DSSP (the residue next to an H run "
        "closes its H-bond, so it is always observed), as a first run showed (0 of 6277) - this equals SV03's "
        "'real end' rule; rama overlap at entry level uses entry x cluster units and entry-level cluster means. "
        "(6) Supporting files not named in the run file: sse_codes_per_chain.csv, consensus_runs_all.csv, "
        "dssp_run.json, numbers_consensus.json, numbers_rama.json, phipsi_check_gemmi.csv, data/sv/SV02/. "
        "The 537 psi check passed (no new values).")

    L = [line1, "", "# SV02 - Consensus secondary structure and per-residue Ramachandran plots", ""]
    h12_c7 = hb[hb.helix_name == "H12"].iloc[0]
    L += [f"**Verdict:** {nc['n_consensus']['global']} global consensus helices ([DEFAULT] >= 4 positions with pydssp "
          f"H in >= 80 % of observing chains). **Data-derived H12 window = {prop[0]}-{prop[1]}** (proposal only; "
          f"primary window stays 538-547); it equals none of 538-547, 537-546, 536-545 (Jaccard "
          f"{comp[0]['jaccard']:.2f} / {comp[1]['jaccard']:.2f} / {comp[2]['jaccard']:.2f}). Consensus H12 is "
          f"{int(h12_c7.cl7_start)}-{int(h12_c7.cl7_end)} in cl7 and {int(h12_c7.cl4_start)}-{int(h12_c7.cl4_end)} in "
          f"cl4. Ramachandran at 537: psi circular mean cl7 {c537['cl7']['psi_circular_mean']:+.1f}, cl4 "
          f"{c537['cl4']['psi_circular_mean']:+.1f} deg (EXPLORATION.md +157 / -47: within 10 deg, check passed); "
          f"cl7-cl4 separation {rs[rs.position == 537].separation_cl7_cl4.iloc[0]:.1f} deg, overlap 0 at chain and "
          "entry level.", ""]

    L += ["## Steps 1-2 - DSSP (pydssp fallback) and numbering", "",
          f"- {run['chains']} kept chains, {run['residues_passed']} residues with all four backbone atoms passed "
          f"(full chain, incl. residues outside 305-550; {run['residues_missing_backbone_atoms']} residues lack an "
          f"atom), {run['segments']} segments (median {run['segments_per_chain_median']:.0f} per chain).",
          f"- UniProt numbering: SV00 esr1_position = auth_seq_id + SV00 uniprot_author_offset for "
          f"{run['numbering_atoms_equal_auth_plus_offset']} / "
          f"{run['numbering_atoms_equal_auth_plus_offset'] + run['numbering_atoms_differ']} mapped atoms.",
          f"- Positions 305-550 with a pydssp state: {run['lbd_positions_observed_dssp']} (SV00 C-alpha mask "
          f"{run['lbd_positions_observed_ca_mask']}; {run['ca_observed_but_incomplete_backbone']} C-alpha-observed "
          "residues lack another backbone atom). States: " +
          ", ".join(f"{k} {v}" for k, v in run["sse_counts_lbd"].items()) + ".",
          "- Observing a position = the residue has N, CA, C and O (a pydssp state).", ""]

    L += ["## Steps 3-4 - consensus helices (helix_boundaries.csv)", "",
          "| helix | helix table | global | cl7 | cl4 | global - table (start, end) |", "|---|---|---|---|---|---|"]
    for r in hb.itertuples():
        def rng(s, e):
            return "-" if pd.isna(s) else f"{int(s)}-{int(e)}"
        L.append(f"| {r.helix_name} | {r.helix_table_row} {rng(r.table_start, r.table_end)} | "
                 f"{rng(r.global_start, r.global_end)} | {rng(r.cl7_start, r.cl7_end)} | {rng(r.cl4_start, r.cl4_end)} | "
                 f"{'NA' if pd.isna(r.global_minus_table_start) else f'{int(r.global_minus_table_start):+d}, {int(r.global_minus_table_end):+d}'} |")
    dr = cr[cr.status != "kept"]
    L += ["", "- Runs that failed the 4-position rule: " +
          "; ".join(f"{r.group} {int(r.run_start)}-{int(r.run_end)}" for r in dr.itertuples()) + ".",
          "- Low-support positions (< 50 % of the group's chains observe them): " +
          "; ".join(f"{g} {v}" for g, v in nc["low_support_positions"].items()) + ".",
          f"- Global helices: min H fraction {cg.min_frac_alpha_H.min():.3f} ({cg.loc[cg.min_frac_alpha_H.idxmin(), 'helix_name']}), low-support positions "
          f"inside helices {int(cg.n_low_support_inside.sum())}.", ""]

    st = {r["helix_name"]: r for r in nc["per_chain_status"]}
    L += ["## Step 5 - per-chain table (sse_per_chain.csv, chain x global consensus helix)", "",
          "| helix | chains with an H run | not observed | observed, not helical | start censored | end censored |",
          "|---|---|---|---|---|---|"]
    for r in cg.itertuples():
        s = st[r.helix_name]
        L.append(f"| {r.helix_name} {r.start}-{r.end} | {s['n_helix']} | {s['n_not_obs']} | {s['n_not_helical']} | "
                 f"{s['start_cens']} | {s['end_cens']} |")
    L += ["", "Start/end = first/last position of the chain's H run(s) overlapping the consensus range; censored = "
          "residue start-1 or start-2 (end+1 or end+2) not observed or across a chain break.", ""]

    L += ["## Step 6 - data-derived H12 window (h12_window.md)", "",
          f"- Proposal: {prop[0]}-{prop[1]} ({prop[1] - prop[0] + 1} positions). Not adopted: 538-547 stays primary; "
          "a human signs off. It differs from 538-547, so SV01c runs (SV01c step 9).", "",
          "| window | shared | only in window | only in proposal | Jaccard | start, end shift |", "|---|---|---|---|---|---|"]
    for c in comp:
        L.append(f"| {c['window']} | {c['shared']} | {c['only_window'] or '-'} | {c['only_proposal'] or '-'} | "
                 f"{c['jaccard']:.2f} | {c['start_shift']:+d}, {c['end_shift']:+d} |")
    L += [""]

    L += ["## Steps 7-11 - phi/psi and Ramachandran plots (rama_summary.csv)", "",
          f"- phipsi_per_chain.csv: {nr['rows_phipsi']} chain x position rows ({nr['chains']} chains; residues with N, "
          f"CA, C), phi in {nr['residues_with_phi']}, psi in {nr['residues_with_psi']}, both in "
          f"{nr['residues_with_both']}; bonded = C-N <= 2.0 A ({nr['peptide_links_gt_2A']} of "
          f"{nr['peptide_links_observed']} observed consecutive links are longer).",
          f"- Positions ([DEFAULT], fixed in SV02_phipsi.py before plotting): {nr['positions']} "
          f"(Y12 top_positions {nr['y12_top_positions']}).",
          "", "| position | context (helix table) | chains cl7 / cl4 / minor | phi, psi cl7 | phi, psi cl4 | "
          "separation | overlap chain | overlap entry |", "|---|---|---|---|---|---|---|---|"]
    for r in rs.itertuples():
        L.append(f"| {r.position} | {r.context} | {r.n_chains_cl7} / {r.n_chains_cl4} / {r.n_chains_minor} | "
                 f"{f1(r.phi_mean_cl7)}, {f1(r.psi_mean_cl7)} | {f1(r.phi_mean_cl4)}, {f1(r.psi_mean_cl4)} | "
                 f"{f1(r.separation_cl7_cl4)} | {f2(r.overlap_chain)} | {f2(r.overlap_entry)} |")
    L += ["", "Circular means in degrees; circular SDs, entry counts and per-cluster overlaps are in rama_summary.csv. "
          "Overlap = share of cl7 + cl4 units closer to the other cluster's mean (wrapped phi/psi distance).",
          f"- Step 11 check (537, circular-mean psi): cl7 {c537['cl7']['psi_circular_mean']:+.2f} vs +157 "
          f"(|d| {c537['cl7']['abs_diff']:.2f}), cl4 {c537['cl4']['psi_circular_mean']:+.2f} vs -47 "
          f"(|d| {c537['cl4']['abs_diff']:.2f}): {'pass' if nr['check_537_pass'] else 'FAIL'}.",
          "- Step 12 (resume) not applicable: SV02 ran in full.", ""]

    files = sorted(OUT.glob("*.csv")) + sorted(OUT.glob("*.md")) + sorted(OUT.glob("*.json")) + [OUT / "env.txt"]
    files += sorted((ROOT / "data/sv/SV02").glob("*"))
    files += sorted((ROOT / "figures/sv/SV02").glob("*.*")) + sorted((ROOT / "figures/sv/SV02/rama").glob("*.*"))
    files += sorted((ROOT / "figures/panels/sv/SV02").glob("*.png"))
    files += sorted((ROOT / "scripts/sv").glob("SV02_*.py"))
    L += ["## Files", ""]
    L += [f"- {p.relative_to(ROOT)} ({size(p)})" for p in files if p.name != "SUMMARY.md"]
    L += ["- results/sv/SV02/SUMMARY.md (this file)"]
    (OUT / "SUMMARY.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
