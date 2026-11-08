"""SV11 (RUN_SV1028): write results/sv/SV11/verdict.md, env.txt and SUMMARY.md (rule 1.7) from the saved SV11 outputs.
Run after SV11_contacts.py, SV11_sites.py, SV11_steroid.py and SV11_figures.py:
.venv/bin/python -B -I scripts/sv/SV11_summary.py
"""
import json
import platform
import sys
from datetime import datetime
from importlib import metadata
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV11"
SITES = ("ligand pocket", "AF-2 groove", "dimer interface")
OUTCOMES = ("share with pocket", "best rank", "best probability")


def size(path):
    n = path.stat().st_size
    return f"{n / 1e6:.2f} MB" if n >= 1e6 else f"{n / 1e3:.1f} kB"


def g(x, d=2):
    return "NA" if pd.isna(x) else f"{x:.{d}g}"


def f2(x):
    return "NA" if pd.isna(x) else f"{x:.2f}"


def main():
    chk = json.loads((OUT / "p2rank_input_check.json").read_text())
    n = json.loads((OUT / "numbers_partA.json").read_text())
    sr = pd.read_csv(OUT / "site_residues.csv", keep_default_na=False, na_values=[""])
    t = pd.read_csv(OUT / "site_tests.csv", keep_default_na=False, na_values=[""])
    sb = pd.read_csv(OUT / "site_by_chain.csv", keep_default_na=False, na_values=[""])
    pl = pd.read_csv(OUT / "pocket_labels.csv", keep_default_na=False, na_values=[""])
    cs = pd.read_csv(OUT / "chain_site_contacts.csv", keep_default_na=False, na_values=[""])
    sc = pd.read_csv(OUT / "steroid_smarts_check.csv", keep_default_na=False, na_values=[""])
    ls = pd.read_csv(OUT / "ligand_steroid.csv", keep_default_na=False, na_values=[""], dtype={"ligand_id": str})
    st = pd.read_csv(OUT / "steroid_tests.csv", keep_default_na=False, na_values=[""]).set_index("level")

    def row(site, outcome, level):
        return t[(t.site == site) & (t.outcome == outcome) & (t.level == level)].iloc[0]

    def differs(site, outcome):
        a, b = row(site, outcome, "chain"), row(site, outcome, "entry")
        return bool(a.tested and b.tested and a.q < 0.05 and b.q < 0.05)

    tested = t[t.tested.astype(str) == "True"]
    c = cs.copy()
    contacts = {s: (int(x.has_partner.astype(str).eq("True").sum()), int(x.sv00_match.astype(str).eq("True").sum()))
                for s, x in c.groupby("site")}
    # groove pockets by peptide presence (descriptive, from site_by_chain.csv)
    gr = sb[sb.site == "AF-2 groove"]
    gp = gr.groupby(["group", "has_site_partner"]).has_pocket.agg(["sum", "size"])
    prob = pl.groupby("label").probability.median()
    rank1 = pl[pl["rank"] == 1]
    lp_extra = pl[(pl.label == "ligand pocket") & (pl["rank"] > 1)]

    L = {s: row(s, "share with pocket", "chain") for s in SITES}
    E = {s: row(s, "share with pocket", "entry") for s in SITES}
    R = {s: row(s, "best rank", "chain") for s in SITES}
    P = {s: row(s, "best probability", "chain") for s in SITES}
    Re = {s: row(s, "best rank", "entry") for s in SITES}
    Pe = {s: row(s, "best probability", "entry") for s in SITES}
    sl, se = st.loc["ligand"], st.loc["entry"]
    part_b_ok = bool(sc.passed.astype(str).eq("True").all())

    def share_txt(s):
        return (f"chains cl7 {f2(L[s].value_cl7)} ({int(L[s].k_cl7)}/{int(L[s].n_cl7)}) vs cl4 {f2(L[s].value_cl4)} "
                f"({int(L[s].k_cl4)}/{int(L[s].n_cl4)}); entries {f2(E[s].value_cl7)} vs {f2(E[s].value_cl4)}")

    gpk = {k: (int(v["sum"]), int(v["size"])) for k, v in gp.iterrows()}
    v = [
        f"# SV11 verdict ([DEFAULT] site sets and labels; 'differs' = q < 0.05 at chain and entry level): the rank-1 "
        f"P2Rank pocket is the ligand pocket in {int((rank1.label == 'ligand pocket').sum())}/{len(rank1)} chains, "
        f"so cl7 and cl4 do not differ there (best probability median {f2(P['ligand pocket'].value_cl7)} vs "
        f"{f2(P['ligand pocket'].value_cl4)}, q {g(P['ligand pocket'].q)} chain / {g(Pe['ligand pocket'].q)} entry).",
        f"- AF-2 groove: {'differs' if differs('AF-2 groove', 'share with pocket') else 'no difference'} - a groove "
        f"pocket in {share_txt('AF-2 groove')} (q {g(L['AF-2 groove'].q)} / {g(E['AF-2 groove'].q)}); P2Rank saw no "
        f"peptide; in cl7, groove pockets occur in {gpk.get(('cl7', True), (0, 0))[0]}/"
        f"{gpk.get(('cl7', True), (0, 0))[1]} peptide-bound chains and {gpk.get(('cl7', False), (0, 0))[0]}/"
        f"{gpk.get(('cl7', False), (0, 0))[1]} without a peptide.",
        f"- Dimer interface: {'differs' if differs('dimer interface', 'share with pocket') else 'no difference'} - "
        f"{share_txt('dimer interface')} (q {g(L['dimer interface'].q)} / {g(E['dimer interface'].q)}); cl4 pockets "
        f"rank better (median {g(R['dimer interface'].value_cl4)} vs {g(R['dimer interface'].value_cl7)}, q "
        f"{g(R['dimer interface'].q)} / {g(Re['dimer interface'].q)}).",
        f"- Groove and dimer pockets are low-probability P2Rank pockets (median probability "
        f"{prob.get('AF-2 groove'):.3f} and {prob.get('dimer interface'):.3f} vs {prob.get('ligand pocket'):.3f} for "
        f"the ligand pocket); the groove rank/probability contrasts rest on {int(R['AF-2 groove'].n_cl4)} cl4 chains in "
        f"{int(Re['AF-2 groove'].n_cl4)} entries.",
        (f"- Part B: steroid ring system in {int(st.n_ligands_flagged_all.iloc[0])} of {int(st.n_ligands_all.iloc[0])} "
         f"ligands; ligands cl7 {int(sl.steroid_cl7)}/{int(sl.n_cl7)} ({f2(sl.share_cl7)}) vs cl4 "
         f"{int(sl.steroid_cl4)}/{int(sl.n_cl4)} ({f2(sl.share_cl4)}), q {g(sl.q)}; entries {int(se.steroid_cl7)}/"
         f"{int(se.n_cl7)} ({f2(se.share_cl7)}) vs {int(se.steroid_cl4)}/{int(se.n_cl4)} ({f2(se.share_cl4)}), q "
         f"{g(se.q)} (ring system only, no action class).") if part_b_ok else
        "- Part B: BLOCKED (SMARTS control check failed, steroid_smarts_check.csv).",
    ]
    (OUT / "verdict.md").write_text("\n".join(v) + "\n")

    env = [f"# SV11 environment, {datetime.now().astimezone().strftime('%Y-%m-%dT%H:%M%z')}",
           f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}"]
    for name in ("gemmi", "numpy", "pandas", "scipy", "statsmodels", "rdkit", "matplotlib"):
        env.append(f"{name} {metadata.version(name)}")
    env += ["P2Rank: not run by SV11 (X08 predictions read; results/p2rank/run_info.json: "
            f"{chk['p2rank_run_info']['p2rank_version']})",
            "nothing installed, upgraded or downgraded by SV11", "network: none"]
    (OUT / "env.txt").write_text("\n".join(env) + "\n")

    k = chk["counts"]
    sha_ok = sum(bool(x["equal"]) for x in chk["s17_sha_check"])
    q_tested = len(tested)
    lig_set = sr[(sr.site == "ligand pocket") & sr.in_site.astype(str).eq("True")]
    line1 = (
        "Not done / substituted: nothing installed or substituted as a tool (gemmi, RDKit, statsmodels present). "
        "(1) Ligand-pocket site set: Q18 results/chain_ligand_of_interest.csv pocket_residues used as they are (step "
        "input line) - scripts/ligand_of_interest.py CONTACT_CUTOFF_A = 4.5 A, heavy atoms, over its 19 canonical "
        "pocket positions only - not 4 A over all residues; the AF-2 groove and dimer sets use 4.0 A on all residues of "
        "the SV00 structures. (2) Declared before running where the run file is open: chain tests = GEE (Gaussian, "
        "independence, entry clusters; share as a linear-probability model, best rank / best probability as average "
        "ranks among the chains with a pocket at the site); entry level = SV10 group rule (>= 80 % of kept chains), "
        "entry has the site if >= 50 % of its chains do (Fisher), entry medians (Mann-Whitney); BH within Part A "
        f"({q_tested} tested rows, both levels) and within Part B (2 tests); 'differs' = q < 0.05 at both levels; ties "
        f"joined by ' + ' ({n['n_tie_pockets']} occurred); SAS points = those of the best-ranked pocket with the label; "
        "bound peptide copy and partner copies grouped as in SV00 (chain, symmetry operation, ASU/symmetry). (3) Kept "
        f"as computed, flagged: the chain-level AF-2 groove best rank / probability tests rest on "
        f"{int(R['AF-2 groove'].n_cl4)} cl4 chains ({int(Re['AF-2 groove'].n_cl4)} entries) that all have rank "
        f"{g(R['AF-2 groove'].value_cl4)}, so their robust p (q {g(R['AF-2 groove'].q)}, {g(P['AF-2 groove'].q)}) is "
        f"anti-conservative; the entry-level Mann-Whitney gives q {g(Re['AF-2 groove'].q)}, {g(Pe['AF-2 groove'].q)}. "
        "(4) data/clean was regenerated by SV00 with the S17 script; P2Rank (X08) read the earlier copy - "
        f"{sha_ok}/{len(chk['s17_sha_check'])} S17 evidence hashes equal and 0 P2Rank residues lie outside the files' "
        "polymer. (5) Supporting files not named in the run file: chain_site_contacts.csv, p2rank_input_check.json, "
        "numbers_partA.json, steroid_smarts_check.csv. The SV11 card and RUN_STATUS row were set up by the "
        "orchestrator, not by this step.")

    md = [line1, "", "# SV11 - Pocket sites and steroid flag", "",
          "**Verdict** " + v[0].lstrip("# ").replace("SV11 verdict ", ""), *v[1:], "",
          "## Inputs (rule 1.10)", "",
          "- X08 card: results/p2rank/<entry>_<label_asym>.pdb_predictions.csv and _residues.csv (788 + 788; "
          "evidence scripts/p2rank_batch.py, run_info.json: "
          f"{chk['p2rank_run_info']['p2rank_version']}, {chk['p2rank_run_info']['n_input_pdb']} inputs, exit code "
          f"{chk['p2rank_run_info']['exit_code']}); X07 card: docs/p2rank_check.md.",
          "- S17 card: data/clean/*.pdb (scripts/clean_chain_pdbs.py); Q18 card: results/chain_ligand_of_interest.csv "
          "and scripts/ligand_of_interest.py (pocket_residues, 4.5 A, 19 canonical positions).",
          "- SV00 card: results/sv/SV00/chain_table.csv (s08_cluster, auth_asym_id, uniprot_author_offset = "
          "chain_offsets.csv), dimer_partner.csv, chain_peptide.csv, the restored mmCIFs in data/raw/ and "
          "data/sv/lbd_extract (chain roles, residue mapping).",
          "- SV10 card: results/sv/SV10/ligand_druglikeness.csv (ligand groups), entry_druglikeness.csv (entry groups; "
          f"recomputed entry groups equal for {n['entry_groups_equal_sv10']}/{n['n_entries']} entries); Q08 card: "
          "results/ligands.csv (SMILES).", "",
          "## Step 1 - what P2Rank saw", "",
          "- scripts/p2rank_batch.py: one P2Rank `predict` run on a dataset listing every data/clean/*.pdb (default "
          "model, -visualizations 0); no other structure file is passed.",
          "- scripts/clean_chain_pdbs.py: each file holds the polymer atoms of the chain's own label_asym_id and its "
          "Q18 ligand instance(s) as HETATM after TER; waters, ions, additives and every other chain (partner "
          "receptor copies and coactivator peptides are separate chains) are not written.",
          f"- Files (p2rank_input_check.json): {k['one_chain_id']}/{k['files']} files have one chain id (= the author "
          f"chain: {k['chain_id_is_author']}); HETATM after TER equals the chain's Q18 ligand(s) in "
          f"{k['het_after_ter_only_q18']}/{k['files']} ({k['apo_files_no_het']} apo files have none); "
          f"{k['n_p2rank_residues']} residues in P2Rank residues.csv and {k['n_pocket_tokens']} pocket residue ids, "
          f"{k['n_p2rank_residues_not_in_polymer'] + k['n_pocket_tokens_not_in_polymer']} outside the file's polymer; "
          f"residue chain-id mismatches {n['residue_chain_mismatch']}. Conclusion: P2Rank saw only the receptor chain "
          "and its ligand of interest - no peptide, no partner chain.",
          f"- Numbering: UniProt = author + SV00 offset; equal to the SV00 per-residue mapping for "
          f"{n['offset_vs_residue_mapping']['equal']}/{n['offset_vs_residue_mapping']['residues']} observed residues.",
          "", "## Step 2 - site residue sets ([DEFAULT] support >= 50 %, fixed before labelling)", ""]
    for s in SITES:
        x = sr[(sr.site == s) & sr.in_site.astype(str).eq("True")]
        md.append(f"- {s}: {len(x)} residues, {int(x.n_chains_relevant.iloc[0])} relevant chains; "
                  + ", ".join(f"{r.wt_aa}{r.uniprot_position} ({r.support:.2f})" for r in x.itertuples()))
    md += [f"- Contacts on the SV00 structures (chain_site_contacts.csv): bound peptide (<= 30 residues; all 419 "
           f"non-receptor polymers in these entries are 8-16-residue peptides) in {contacts['AF-2 groove'][0]} chains, "
           f"equal to SV00 peptide_bound for {contacts['AF-2 groove'][1]}/788; dimer partner in "
           f"{contacts['dimer interface'][0]} chains, H11 atom count equal to SV00 for "
           f"{contacts['dimer interface'][1]}/788.",
           f"- Overlap between sets: {'; '.join(f'{k2}: {len(v2)}' for k2, v2 in n['overlap'].items())}. "
           f"Ligand-pocket set = {len(lig_set)} of the 19 Q18 positions.", "",
           "## Steps 3-4 - pocket labels", "",
           f"- {n['n_pockets']} pockets in {n['n_chains']} chains (per chain {n['n_pockets_per_chain']['min']}-"
           f"{n['n_pockets_per_chain']['max']}, median {n['n_pockets_per_chain']['median']:g}). Labels: "
           + ", ".join(f"{k2} {v2}" for k2, v2 in n["label_counts"].items()) + f"; ties {n['n_tie_pockets']}.",
           f"- Rank-1 pocket = ligand pocket in {int((rank1.label == 'ligand pocket').sum())}/{len(rank1)} chains "
           f"(share of its residues in the set: median {rank1.share.median():.2f}, min {rank1.share.min():.2f}); "
           f"{len(lp_extra)} more ligand-pocket pockets at ranks {sorted(lp_extra['rank'].unique().tolist())} "
           f"(SAS points {int(lp_extra.sas_points.min())}-{int(lp_extra.sas_points.max())}).",
           "- Median P2Rank probability by label: " + ", ".join(f"{k2} {v2:.3f}" for k2, v2 in prob.items()) + ".",
           "", "## Step 5 - cl7 vs cl4 per site (site_tests.csv)", "",
           "| site | outcome | chain cl7 | chain cl4 | q chain | entry cl7 | entry cl4 | q entry | differs |",
           "|---|---|---|---|---|---|---|---|---|"]
    for s in SITES:
        for o in OUTCOMES:
            a, b = row(s, o, "chain"), row(s, o, "entry")
            if o == "share with pocket":
                ca = f"{f2(a.value_cl7)} ({int(a.k_cl7)}/{int(a.n_cl7)})"
                cb = f"{f2(a.value_cl4)} ({int(a.k_cl4)}/{int(a.n_cl4)})"
                ea = f"{f2(b.value_cl7)} ({int(b.k_cl7)}/{int(b.n_cl7)})"
                eb = f"{f2(b.value_cl4)} ({int(b.k_cl4)}/{int(b.n_cl4)})"
            else:
                ca, cb = f"{a.iqr_cl7} n {int(a.n_cl7)}", f"{a.iqr_cl4} n {int(a.n_cl4)}"
                ea, eb = f"{b.iqr_cl7} n {int(b.n_cl7)}", f"{b.iqr_cl4} n {int(b.n_cl4)}"
            qa = g(a.q) if a.tested else a.note
            qb = g(b.q) if b.tested else b.note
            md.append(f"| {s} | {o} | {ca} | {cb} | {qa} | {ea} | {eb} | {qb} | {'yes' if differs(s, o) else 'no'} |")
    md += ["", "Chain level: GEE with entry clusters (share: linear-probability model; rank / probability: average ranks "
           "among chains with such a pocket). Entry level: Fisher (entry has the site if >= 50 % of its chains do) and "
           "Mann-Whitney on entry medians. BH over all tested Part A rows. 'other surface' rows in site_tests.csv are "
           "descriptive (tested = False).",
           f"- Descriptive (site_by_chain.csv, has_site_partner): AF-2 groove pockets in cl7 chains with a bound "
           f"peptide {gpk.get(('cl7', True), (0, 0))[0]}/{gpk.get(('cl7', True), (0, 0))[1]}, without "
           f"{gpk.get(('cl7', False), (0, 0))[0]}/{gpk.get(('cl7', False), (0, 0))[1]}; cl4 with "
           f"{gpk.get(('cl4', True), (0, 0))[0]}/{gpk.get(('cl4', True), (0, 0))[1]}, without "
           f"{gpk.get(('cl4', False), (0, 0))[0]}/{gpk.get(('cl4', False), (0, 0))[1]}.",
           "", "## Part B - steroid ring system", "",
           "- SMARTS check (steroid_smarts_check.csv): " + ", ".join(
               f"{r.ligand_id} {'match' if str(r.observed_match) == 'True' else 'no match'} "
               f"({'pass' if str(r.passed) == 'True' else 'FAIL'})" for r in sc.itertuples()) + ".",
           f"- Flagged {int(st.n_ligands_flagged_all.iloc[0])} of {int(st.n_ligands_all.iloc[0])} SV10 ligands "
           f"(SMILES parsed for {int(ls.flag_reason.isna().sum())}): "
           + ", ".join(sorted(ls.ligand_id[ls.steroid_ring_system.astype(str) == "True"])) + ".",
           f"- Ligand level: cl7 {int(sl.steroid_cl7)}/{int(sl.n_cl7)} ({f2(sl.share_cl7)}) vs cl4 "
           f"{int(sl.steroid_cl4)}/{int(sl.n_cl4)} ({f2(sl.share_cl4)}), Fisher OR {sl.odds_ratio_cl7_vs_cl4:.2f}, "
           f"p {g(sl.p)}, q {g(sl.q)}. Entry level: cl7 {int(se.steroid_cl7)}/{int(se.n_cl7)} ({f2(se.share_cl7)}) "
           f"vs cl4 {int(se.steroid_cl4)}/{int(se.n_cl4)} ({f2(se.share_cl4)}), OR {se.odds_ratio_cl7_vs_cl4:.2f}, "
           f"p {g(se.p)}, q {g(se.q)}. cl4 steroid ligands: {sl.steroid_ligands_cl4}.",
           "- The flag is its own column (steroid_ring_system); ligand_chemotype from SV00 is copied unchanged "
           "(ligand_chemotype_sv00). It names the ring system only.", "",
           "## Figure (backup)", "",
           "- figures/sv/SV11/SV11_site_share.png / .svg (panel figures/panels/sv/SV11/SV11_site_share__share.png), "
           "from site_tests.csv. Caption: share of chains with a pocket at the site - " + "; ".join(
               f"{s}: cl7 {f2(L[s].value_cl7)}, cl4 {f2(L[s].value_cl4)}" for s in SITES) + ".", "",
           "## Files", ""]
    files = [OUT / x for x in ("site_residues.csv", "pocket_labels.csv", "site_by_chain.csv", "site_tests.csv",
                               "verdict.md", "ligand_steroid.csv", "steroid_tests.csv", "steroid_smarts_check.csv",
                               "chain_site_contacts.csv", "p2rank_input_check.json", "numbers_partA.json", "env.txt")]
    files += [ROOT / "figures/sv/SV11/SV11_site_share.png", ROOT / "figures/sv/SV11/SV11_site_share.svg",
              ROOT / "figures/panels/sv/SV11/SV11_site_share__share.png"]
    files += [ROOT / "scripts/sv" / x for x in ("SV11_contacts.py", "SV11_sites.py", "SV11_steroid.py",
                                                "SV11_figures.py", "SV11_summary.py")]
    md += [f"- {p.relative_to(ROOT)} ({size(p)})" for p in files] + [
           "- results/sv/SV11/SUMMARY.md (this file)",
           "- docs/EXPLORATION.md: appended section 'Suggested, not done (SV11, 2026-11-08)' (rule 1.1; existing text "
           "unchanged)", "",
           "Run order: .venv/bin/python -B -I scripts/sv/SV11_contacts.py; SV11_sites.py; SV11_steroid.py; "
           "SV11_figures.py; SV11_summary.py"]
    (OUT / "SUMMARY.md").write_text("\n".join(md) + "\n")
    print((OUT / "SUMMARY.md").read_text())


if __name__ == "__main__":
    main()
