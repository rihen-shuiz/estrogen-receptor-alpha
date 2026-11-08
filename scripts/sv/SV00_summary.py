"""SV00: write results/sv/SV00/SUMMARY.md (rule 1.7) from the saved SV00 files only.
Line 1 (not done / substituted) is passed in NOT_DONE below and must list everything substituted in this step.
Run: .venv/bin/python -I scripts/sv/SV00_summary.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV00"
NOT_DONE = (
    "Not done / substituted: (1) scripts/fetch_data.py was not run as a whole - it rewrites data/raw/download_manifest.csv "
    "and writes data/ccd and data/lbd, outside rule 1.2 - so its download function (esr_analysis.download_structures."
    "download_one: same URL, plain .cif format, validation) filled data/raw, the manifest went to data/sv/manifest.csv, "
    "data/ccd and data/lbd were not rebuilt, and the Q04 extraction was regenerated into data/sv/lbd_extract; "
    "(2) the RCSB mmCIFs carry no _pdbx_sifts_unp_segments, so offsets come from the section-3 SIFTS fallback (PDBe API), "
    "and for 1 chain (6IAR|A|1, SIFTS segment lengths disagree) from the deposited _struct_ref_seq; (3) the original "
    "download date is not recorded in the data and was taken as 2026-10-14 (Q03 DONE log line that recorded "
    "data/raw/sha256sums.txt); (4) the optional fresh RCSB search for newer entries was not run; (5) pH stated only in "
    "the crystal-growth text (65 entries that Y02 parsed) stays NaN, as section 3 reads only _exptl_crystal_grow.pH "
    "and the Data API.")


def rd(p, **kw):
    return pd.read_csv(p, keep_default_na=False, na_values=[""], low_memory=False, **kw)


def size(p):
    p = ROOT / p
    if p.is_dir():
        files = [f for f in p.rglob("*") if f.is_file()]
        return f"{len(files)} files, {sum(f.stat().st_size for f in files) / 1e6:.1f} MB"
    return f"{p.stat().st_size / 1e3:.1f} kB" if p.stat().st_size < 1e6 else f"{p.stat().st_size / 1e6:.1f} MB"


def main():
    gate = json.loads((OUT / "gate.json").read_text())
    counts = rd(OUT / "counts.csv")
    man = rd(ROOT / "data/sv/manifest.csv")
    t = rd(OUT / "chain_table.csv")
    miss = rd(OUT / "missing_by_column.csv")
    clean = rd(OUT / "metadata_cleaning.csv")
    cross = rd(OUT / "api_crosscheck.csv")
    offs = rd(OUT / "chain_offsets.csv")
    exc = rd(OUT / "excluded_chains.csv")
    cdiff = rd(OUT / "coordinate_diff_by_chain.csv")
    mask = json.loads((OUT / "mask_check.json").read_text())
    ext = rd(OUT / "extract_compare.csv")
    vcx = rd(OUT / "value_crosscheck.csv")
    dim = rd(OUT / "dimer_partner.csv")
    clean_dir = ROOT / "data/clean"
    n_clean = len(list(clean_dir.glob("*.pdb"))) if clean_dir.exists() else 0

    L = [NOT_DONE, "", "# SV00 - Data restore and shared chain table", "",
         f"**Verdict:** reproduction gate {gate['verdict']} (ARI at 1.5 A = {gate['ari_by_cut']['1.5']}, max |dRMSD| = "
         f"{gate['max_abs_delta_rmsd_A']:.2g} A, {gate['n_differing_chains']} differing chains); chain table {len(t)} rows; "
         f"missing values reported per column. Done-when criteria met.", ""]
    L += ["## Restore (data/sv/manifest.csv)", "",
          f"- Frozen entry list: data/pdb_ids_esr1.txt + data/pdb_ids_esr2.txt ({len(man)} IDs; tracker Q01/Q02 out, Q03 "
          f"log ev data/raw/download_manifest.csv + data/raw/sha256sums.txt). No RCSB search was run.",
          f"- Downloaded OK: {int(man.download_ok.astype(str).eq('True').sum())}/{len(man)}; required entries (hold the 788 "
          f"kept chains) {int(man[man.required.astype(str).eq('True')].download_ok.astype(str).eq('True').sum())}/"
          f"{int(man.required.astype(str).eq('True').sum())}.",
          f"- sha256 equal to data/raw/sha256sums.txt: {int(man.sha_matches_frozen.astype(str).eq('True').sum())}/{len(man)}; "
          f"latest revision after the original download date ({man.original_download_date.iloc[0]}, Q03 log line): "
          f"{int(man.revised_after_original.astype(str).eq('True').sum())} entries "
          f"({', '.join(man.entry[man.revised_after_original.astype(str).eq('True')].tolist()) or 'none'}).",
          f"- Changed files (sha256 differs): "
          f"{', '.join(man.entry[~man.sha_matches_frozen.astype(str).eq('True')].tolist()) or 'none'}.", ""]
    L += ["## Counts vs the original run (results/sv/SV00/counts.csv)", "", "| quantity | run file | frozen tables | regenerated | match |",
          "|---|---|---|---|---|"]
    L += [f"| {r.quantity} | {r.run_file} | {r.frozen_tables} | {r.regenerated} | {r.match} |" for r in counts.itertuples()]
    reasons = exc[exc.present_in != "frozen_only"].eligibility_reason.value_counts().to_dict()
    L += ["", f"Excluded ESR1 chains with reasons (results/sv/SV00/excluded_chains.csv): {len(exc)} rows; regenerated "
              f"reasons {reasons}.", ""]
    L += ["## Reproduction gate (results/sv/SV00/gate.json)", "",
          f"- Code path: Q04 extraction (esr_analysis.extract_lbd.inventory_entry) -> S03 fit_all (core "
          f"{gate['core']}, reference {gate['reference_chain']}) -> S04 compute_distances -> S08 structural_groups.",
          f"- Chains: recomputed {gate['n_chains_recomputed']}, frozen {gate['n_chains_frozen']}, compared "
          f"{gate['n_chains_compared']}; distance exclusions {gate['n_distance_exclusions']}.",
          f"- Kept chains with C-alpha coordinates different from the frozen residue mapping: "
          f"{gate['n_chains_coordinates_differ']} (results/sv/SV00/coordinate_diff_by_chain.csv; max |dxyz| "
          f"{np.nanmax(cdiff.max_abs_dxyz_A.astype(float)):.3g} A).",
          f"- Aligned coordinates vs results/alignment/aligned_coordinates.npz: max |d| {gate['aligned_coords_max_abs_diff_A']:.3g} A, "
          f"same gaps {gate['aligned_coords_same_gaps']}.",
          f"- **ARI at 1.5 A = {gate['ari_by_cut']['1.5']}**; ARI at 1.0/2.0/2.5 A = {gate['ari_by_cut']['1.0']}/"
          f"{gate['ari_by_cut']['2.0']}/{gate['ari_by_cut']['2.5']}; labels identical per cut {gate['labels_identical_by_cut']}.",
          f"- **Max |dRMSD| = {gate['max_abs_delta_rmsd_A']:.3g} A** ({gate['n_pairs_delta_gt_1e-6']} pairs differ by > 1e-6 A).",
          f"- Clusters per cut recomputed {gate['clusters_per_cut_recomputed']}, frozen {gate['clusters_per_cut_frozen']}.",
          f"- Differing chains: {gate['n_differing_chains']} ({', '.join(gate['differing_chains']) or 'none'}); "
          f"remediated: {len(gate['remediated_chains'])} ({', '.join(gate['remediated_chains']) or 'none'}).",
          f"- **Verdict: {gate['verdict']}.**", ""]
    L += ["## Regenerated data", "",
          "- Q04 extraction re-run on the restored files vs the frozen tables (extract_compare.csv): " + "; ".join(
              f"{Path(r.frozen).name} {'byte-identical' if str(r.byte_identical) == 'True' else 'max numeric diff %.2g (%s), %d text cells differ' % (r.max_abs_numeric_diff, r.numeric_columns_differing, r.text_cells_differing)}"
              for r in ext.itertuples()) + ".",
          f"- data/clean: {n_clean} PDB files written by the existing scripts/clean_chain_pdbs.py (S17); the folder did not exist.",
          f"- Backbone N/CA/C/O of the 788 kept chains: data/sv/backbone/backbone.npz (deposited frame, ESR1 305-550) and "
          f"backbone_atoms.csv.gz.",
          f"- Observation mask results/sv/SV00/observed_mask.npy {mask['shape']}: {mask['observed_ca']} observed C-alpha; "
          f"rows differing from the frozen mapping: {mask['rows_differ_from_frozen_mapping']}.", ""]
    vc = lambda c: t[c].value_counts(dropna=False).to_dict()  # noqa: E731
    L += ["## Chain table (results/sv/SV00/chain_table.csv)", "",
          f"- {len(t)} rows x {len(t.columns)} columns; sources in column_sources.md; missing values per column in "
          f"missing_by_column.csv.",
          f"- ligand_chemotype: {vc('ligand_chemotype')}; carboxylic acid on {int(t.ligand_carboxylic_acid.astype(str).eq('True').sum())} "
          f"chains, basic amine on {int(t.ligand_basic_amine.astype(str).eq('True').sum())}.",
          f"- Dimer partner source: {vc('dimer_partner_source')}; partner is the chain's own symmetry image: "
          f"{int(t.dimer_partner_self_image.astype(str).eq('True').sum())}; H11 heavy atoms in contact with the partner "
          f"min {int(dim.dimer_partner_h11_atoms.min())}, median {dim.dimer_partner_h11_atoms.median():g}; second-best "
          f"copy above half the best: {int((dim.dimer_second_best_h11_atoms > 0.5 * dim.dimer_partner_h11_atoms).sum())}; "
          f"partner not a kept chain (no cluster): {', '.join(t.dimer_partner_chain_key[t.dimer_partner_cluster.isna()]) or 'none'}.",
          f"- Chains with two distinct ligand-of-interest codes (heavy atoms and flags from the primary, first code): "
          f"{int((t.n_ligand_codes > 1).sum())} ({', '.join(t.chain_key[t.n_ligand_codes > 1])}).",
          f"- Chain-level bound peptide: {int(t.peptide_bound.astype(str).eq('True').sum())} chains "
          f"({vc('peptide_bound_source')}); entry-level peptide_present (Y02): "
          f"{int(t.peptide_present_entry.astype(str).eq('True').sum())} chains.",
          f"- UniProt-author offset: {offs.uniprot_author_offset.value_counts(dropna=False).to_dict()}; source "
          f"{offs.offset_source.value_counts().to_dict()}; constant over the LBD {offs.offset_sifts_constant.astype(str).value_counts().to_dict()}.",
          f"  Deposited _struct_ref_seq differs from SIFTS for "
          f"{int((offs.offset_sifts.notna() & (offs.offset_sifts != offs.offset_struct_ref_seq)).sum())} chains in "
          f"{offs[offs.offset_sifts.notna() & (offs.offset_sifts != offs.offset_struct_ref_seq)].entry_id.nunique()} entries "
          f"(struct_ref_seq offsets {offs[offs.offset_sifts.notna() & (offs.offset_sifts != offs.offset_struct_ref_seq)].offset_struct_ref_seq.value_counts().to_dict()}); "
          f"SIFTS equals the existing residue mapping (esr1_position - auth_seq_id) for "
          f"{int((offs.offset_sifts == offs.offset_existing_mapping).sum())}/{int(offs.offset_sifts.notna().sum())} chains.",
          f"- multi_valued entries: {int(t.drop_duplicates('entry_id').multi_valued.astype(str).eq('True').sum())}; "
          f"remediated chains: {int(t.remediated.astype(str).eq('True').sum())}.", "",
          "Metadata cleaning (metadata_cleaning.csv):", "", "| field | range | missing in mmCIF | out of range -> NaN | filled from API | missing final |",
          "|---|---|---|---|---|---|"]
    L += [f"| {r.field} | {r.plausible_range if isinstance(r.plausible_range, str) else ''} | {r.n_missing_in_mmcif} | {r.n_out_of_range_set_nan} "
          f"{('(' + r.out_of_range_entries + ')') if isinstance(r.out_of_range_entries, str) and r.field != 'multi_valued' else ''} | "
          f"{r.n_filled_from_api} | {r.n_missing_final} |" for r in clean.itertuples() if r.field != "multi_valued"]
    L += ["", "Value cross-checks, entry level (value_crosscheck.csv): " + "; ".join(
        f"{r.column} vs {r.compared_with}: {r.n_agree}/{r.n_both} agree, only other {r.n_only_other}" for r in vcx.itertuples())]
    L += ["", "Cross-checks (api_crosscheck.csv): " + "; ".join(f"{r.field} {r.n_agree}/{r.n_both} agree" for r in cross.itertuples()), ""]
    L += ["Missing values per column (columns with any missing):", "", "| column | missing | fraction |", "|---|---|---|"]
    L += [f"| {r.column} | {r.n_missing} | {r.fraction_missing} |" for r in miss.itertuples() if r.n_missing > 0]
    files = ["data/sv/manifest.csv", "data/sv/fetch.log", "data/sv/lbd_extract", "data/sv/lbd_chains.csv",
             "data/sv/backbone", "data/sv/gate", "data/sv/api_cache", "data/sv/sifts_cache", "data/clean",
             "results/sv/SV00/chain_table.csv", "results/sv/SV00/column_sources.md",
             "results/sv/SV00/missing_by_column.csv", "results/sv/SV00/observed_mask.npy",
             "results/sv/SV00/observed_mask_rows.csv", "results/sv/SV00/observed_mask_cols.csv",
             "results/sv/SV00/mask_check.json", "results/sv/SV00/gate.json", "results/sv/SV00/gate_chains.csv",
             "results/sv/SV00/counts.csv", "results/sv/SV00/excluded_chains.csv",
             "results/sv/SV00/coordinate_diff_by_chain.csv", "results/sv/SV00/rmsd_recomputed.npy",
             "results/sv/SV00/rmsd_recomputed_index.csv", "results/sv/SV00/clusters_recomputed.csv",
             "results/sv/SV00/entry_metadata_mmcif.csv", "results/sv/SV00/metadata_cleaning.csv",
             "results/sv/SV00/api_crosscheck.csv", "results/sv/SV00/chain_offsets.csv",
             "results/sv/SV00/dimer_partner.csv", "results/sv/SV00/chain_peptide.csv", "results/sv/SV00/extract_compare.csv",
             "results/sv/SV00/value_crosscheck.csv", "results/sv/SV00/gate_stdout.txt", "results/sv/SV00/env.txt",
             "data/sv/fetch_stdout.txt", "data/sv/extract_stdout.txt", "data/sv/clean_chain_pdbs_stdout.txt"]
    L += ["", "## Files", ""] + [f"- {f} ({size(f)})" for f in files if (ROOT / f).exists()]
    L += ["- data/raw/*.cif (" + str(len(list((ROOT / 'data/raw').glob('*.cif')))) + " files added; download_manifest.csv "
          "and sha256sums.txt untouched)",
          "- scripts: scripts/sv/SV00_fetch.py, SV00_extract.py, SV00_gate.py, SV00_api_test.py, SV00_metadata.py, "
          "SV00_extract_compare.py, SV00_contacts.py, SV00_chain_table.py, SV00_env.py, SV00_summary.py"]
    (OUT / "SUMMARY.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
