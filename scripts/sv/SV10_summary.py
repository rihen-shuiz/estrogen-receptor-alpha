"""SV10 (RUN_SV1028): write results/sv/SV10/SUMMARY.md (rule 1.7) from the saved SV10 outputs only.
Run after scripts/sv/SV10_druglikeness.py: .venv/bin/python -B -I scripts/sv/SV10_summary.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV10"


def size(path):
    n = path.stat().st_size
    return f"{n / 1e6:.2f} MB" if n >= 1e6 else f"{n / 1e3:.1f} kB"


def main():
    n = json.loads((OUT / "numbers.json").read_text())
    w = pd.read_csv(OUT / "druglikeness_by_group.csv", keep_default_na=False, na_values=[""])
    verdict = (OUT / "verdict.md").read_text().strip().splitlines()
    ha_l, ha_e = n["heavy_atoms_ligand"], n["heavy_atoms_entry"]

    def pq(x):
        return "NA" if pd.isna(x) else f"{x:.2g}"

    line1 = (
        f"Not done / substituted: nothing substituted (RDKit and statsmodels present; {n['n_parsed']} of "
        f"{n['n_ligands']} SMILES parsed, {n['n_unparsed']} empty rows). (1) results/potency_by_cluster.csv is not an "
        "`out` of a card: it is evidence of the X16 DONE line in tracker/esr1.log and of tracker/prov/X15.json and "
        "X16.json (cards X15/X16 list figures/potency_by_cluster.png and results/potency_tests.csv). (2) Declared before "
        "running where the run file is open: PAINS = one count over PAINS_A + PAINS_B + PAINS_C (per-family counts in "
        "ligand_druglikeness.csv); QED ALERTS is treated as an alert count (Fisher); the size-check OLS ranks the count "
        "values themselves; entry group = the same >= 80 % rule on the entry's kept chains; entry ligand = the "
        f"ligand_id_primary of most of its chains ({n['n_entries_ligand_tie']} ties resolved by the first chain); the md "
        "table shows rho at ligand level (entry rho in the csv); a no-variation rule was declared and not triggered. "
        f"(3) {len(n['qed_property_identical_to'])} QED properties equal a Lipinski/Veber row for every ligand "
        f"({', '.join(f'{q} = {b}' for q, b in n['qed_property_identical_to'].items())}); kept as written, so each BH "
        f"family has {n['n_parameters']} rows. (4) Supporting files not named in the run file: entry_druglikeness.csv, "
        "numbers.json.")

    beyond = w[w.label == "adds beyond size"]
    distinct = beyond[beyond.identical_to.fillna("") == ""]
    md = [line1, "", "# SV10 - Ligand drug-likeness, cl7 vs cl4 with a size check", "",
          f"**Verdict ([DEFAULT] size check): {len(beyond)} of {len(w)} parameters add beyond size "
          f"({len(distinct)} distinct quantities: {', '.join(distinct.parameter)}); {len(n['restates_size'])} restate "
          "size.** These parameters describe oral drug-likeness of the crystallised ligands, not binding, conformation "
          "or activity.", "",
          "## Inputs (rule 1.10)", "",
          "- Q08 card out results/ligands.csv: canonical_isomeric_smiles for all "
          f"{n['smiles_column_counts'].get('canonical_isomeric_smiles', 0)} ligands (smiles fallback not needed).",
          "- Y09 card out results/ext/ligand_descriptors.csv: ligand set, structural_majority (>= 80 % rule), "
          f"heavy_atoms (RDKit heavy atoms of the SMILES equal it for {n['heavy_atoms_rdkit_equal_y09']} of "
          f"{n['n_ligands']}). Groups recomputed from the SV00 chain table agree for "
          f"{n['ligand_groups_recomputed_equal']} of {n['ligand_groups_total']} ligands.",
          f"- results/potency_by_cluster.csv (X15/X16 evidence): pchembl for {n['n_potency_rows']} ligands, one value each.",
          "- SV00 card out results/sv/SV00/chain_table.csv: entry_id, s08_cluster and ligand_id_primary (first row of "
          "Q18 results/chain_ligand_of_interest.csv per chain) for the entry level.", "",
          "## Key numbers", "",
          f"- Ligands: {n['n_ligands']} described; groups " + ", ".join(f"{k} {v}" for k, v in
                                                                        n['ligand_group_counts'].items())
          + f". Compared: cl7 {n['n_ligands_compared']['cl7']}, cl4 {n['n_ligands_compared']['cl4']}; with pChEMBL "
          f"cl7 {n['n_ligands_with_pchembl_compared']['cl7']}, cl4 {n['n_ligands_with_pchembl_compared']['cl4']} "
          f"(of {n['n_ligands_with_pchembl']}).",
          f"- Entries: {n['n_entries']}; groups " + ", ".join(f"{k} {v}" for k, v in n['entry_group_counts'].items())
          + f"; {n['n_entries_without_ligand']} without a ligand of interest; compared cl7 "
          f"{n['n_entries_compared']['cl7']}, cl4 {n['n_entries_compared']['cl4']} "
          f"({n['n_distinct_ligands_entry_level']} distinct ligands; {n['n_entries_multiple_ligands']} entries hold "
          f">1 ligand of interest; in {n['n_entries_group_differs_from_ligand_group']} the entry group differs from its "
          "ligand's group).",
          f"- Size: heavy atoms median [IQR] cl7 {ha_l['cl7']['median']:.0f} [{ha_l['cl7']['q1']:.0f}-"
          f"{ha_l['cl7']['q3']:.0f}] vs cl4 {ha_l['cl4']['median']:.0f} [{ha_l['cl4']['q1']:.0f}-"
          f"{ha_l['cl4']['q3']:.0f}] (ligands, Mann-Whitney p {n['heavy_atoms_mw_p']['ligand']:.2g}); entries "
          f"{ha_e['cl7']['median']:.0f} vs {ha_e['cl4']['median']:.0f} (p {n['heavy_atoms_mw_p']['entry']:.2g}).",
          f"- Unadjusted: q < 0.05 for {n['n_q_lt_005']['ligand']} of {n['n_parameters']} parameters at ligand level, "
          f"{n['n_q_lt_005']['entry']} at entry level. Size-adjusted q < 0.05: {n['n_adj_q_lt_005']['ligand']} "
          f"ligand, {n['n_adj_q_lt_005']['entry']} entry.",
          f"- Alerts (all {n['n_ligands']} ligands): Brenk >= 1 in {n['n_ligands_brenk_ge1']}, PAINS >= 1 in "
          f"{n['n_ligands_pains_ge1']} (" + ", ".join(f"{k} x{v}" for k, v in n['pains_names_all'].items()) + ").",
          "- LE and LipE use pChEMBL, which mixes IC50, Ki and Kd from different assays.", "",
          "| parameter | cl7 (ligands) | cl4 (ligands) | q ligand | q entry | rho HA | size-adj q lig / entry | direction at equal HA | label |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in w.itertuples():
        rho = "NA" if pd.isna(r.lig_rho_heavy_atoms) else f"{r.lig_rho_heavy_atoms:+.2f}"
        direction = r.adj_direction.replace(" at equal heavy atoms", "") if r.label == "adds beyond size" else ""
        md.append(f"| {r.parameter}{' (= ' + r.identical_to + ')' if isinstance(r.identical_to, str) and r.identical_to else ''} "
                  f"| {r.lig_cl7} | {r.lig_cl4} | {pq(r.lig_q)} | {pq(r.ent_q)} | {rho} | {pq(r.lig_adj_q)} / "
                  f"{pq(r.ent_adj_q)} | {direction} | {r.label} |")
    md += ["", "Full table with definitions, entry-level n and medians, p values and test statistics: "
           "druglikeness_by_group.csv / druglikeness_by_group.md.", "", "## verdict.md", ""]
    md += [f"> {x}" for x in verdict]
    files = [OUT / f for f in ("ligand_druglikeness.csv", "entry_druglikeness.csv", "druglikeness_by_group.csv",
                               "druglikeness_by_group.md", "verdict.md", "numbers.json", "env.txt")]
    files += [ROOT / "scripts/sv/SV10_druglikeness.py", ROOT / "scripts/sv/SV10_summary.py"]
    md += ["", "## Files", ""] + [f"- {f.relative_to(ROOT)} ({size(f)})" for f in files]
    md += ["- results/sv/SV10/SUMMARY.md (this file)", "- No figure (run file: 'No figure').", "",
           "Run order: .venv/bin/python -B -I scripts/sv/SV10_druglikeness.py; scripts/sv/SV10_summary.py"]
    (OUT / "SUMMARY.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
