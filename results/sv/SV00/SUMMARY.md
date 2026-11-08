Not done / substituted: (1) scripts/fetch_data.py was not run as a whole - it rewrites data/raw/download_manifest.csv and writes data/ccd and data/lbd, outside rule 1.2 - so its download function (esr_analysis.download_structures.download_one: same URL, plain .cif format, validation) filled data/raw, the manifest went to data/sv/manifest.csv, data/ccd and data/lbd were not rebuilt, and the Q04 extraction was regenerated into data/sv/lbd_extract; (2) the RCSB mmCIFs carry no _pdbx_sifts_unp_segments, so offsets come from the section-3 SIFTS fallback (PDBe API), and for 1 chain (6IAR|A|1, SIFTS segment lengths disagree) from the deposited _struct_ref_seq; (3) the original download date is not recorded in the data and was taken as 2026-10-14 (Q03 DONE log line that recorded data/raw/sha256sums.txt); (4) the optional fresh RCSB search for newer entries was not run; (5) pH stated only in the crystal-growth text (65 entries that Y02 parsed) stays NaN, as section 3 reads only _exptl_crystal_grow.pH and the Data API.

# SV00 - Data restore and shared chain table

**Verdict:** reproduction gate pass (ARI at 1.5 A = 1.0, max |dRMSD| = 9.2e-15 A, 0 differing chains); chain table 788 rows; missing values reported per column. Done-when criteria met.

## Restore (data/sv/manifest.csv)

- Frozen entry list: data/pdb_ids_esr1.txt + data/pdb_ids_esr2.txt (526 IDs; tracker Q01/Q02 out, Q03 log ev data/raw/download_manifest.csv + data/raw/sha256sums.txt). No RCSB search was run.
- Downloaded OK: 526/526; required entries (hold the 788 kept chains) 351/351.
- sha256 equal to data/raw/sha256sums.txt: 526/526; latest revision after the original download date (2026-10-14, Q03 log line): 0 entries (none).
- Changed files (sha256 differs): none.

## Counts vs the original run (results/sv/SV00/counts.csv)

| quantity | run file | frozen tables | regenerated | match |
|---|---|---|---|---|
| esr1_entries | 487 | 487 | 487 | True |
| esr1_chains | 901 | 901 | 901 | True |
| esr1_chain_entries | 444 | 444 | 444 | True |
| kept_chains | 788 | 788 | 788 | True |
| kept_entries | 351 | 351 | 351 | True |
| excluded_chains | 113 | 113 | 113 | True |
| esr2_entries | 39 | 39 | 39 | True |
| esr2_chains | 86 | 86 | 86 | True |

Excluded ESR1 chains with reasons (results/sv/SV00/excluded_chains.csv): 113 rows; regenerated reasons {'non_LBD_fragment': 109, 'short_construct_LBD': 2, 'incomplete_observed_LBD': 2}.

## Reproduction gate (results/sv/SV00/gate.json)

- Code path: Q04 extraction (esr_analysis.extract_lbd.inventory_entry) -> S03 fit_all (core alignment_core_v3, reference 7NFB|B|1) -> S04 compute_distances -> S08 structural_groups.
- Chains: recomputed 788, frozen 788, compared 788; distance exclusions 0.
- Kept chains with C-alpha coordinates different from the frozen residue mapping: 0 (results/sv/SV00/coordinate_diff_by_chain.csv; max |dxyz| 0 A).
- Aligned coordinates vs results/alignment/aligned_coordinates.npz: max |d| 3.55e-14 A, same gaps True.
- **ARI at 1.5 A = 1.0**; ARI at 1.0/2.0/2.5 A = 1.0/1.0/1.0; labels identical per cut {'1.0': True, '1.5': True, '2.0': True, '2.5': True}.
- **Max |dRMSD| = 9.21e-15 A** (0 pairs differ by > 1e-6 A).
- Clusters per cut recomputed {'1.0': 16, '1.5': 7, '2.0': 4, '2.5': 3}, frozen {'1.0': 16, '1.5': 7, '2.0': 4, '2.5': 3}.
- Differing chains: 0 (none); remediated: 0 (none).
- **Verdict: pass.**

## Regenerated data

- Q04 extraction re-run on the restored files vs the frozen tables (extract_compare.csv): residue_mapping.csv byte-identical; receptor_chains.csv byte-identical; entries.csv max numeric diff 7.1e-15 (mean_B_factor_A2), 0 text cells differ; lbd_chains.csv max numeric diff 1.4e-14 (mean_CA_B_A2), 0 text cells differ.
- data/clean: 788 PDB files written by the existing scripts/clean_chain_pdbs.py (S17); the folder did not exist.
- Backbone N/CA/C/O of the 788 kept chains: data/sv/backbone/backbone.npz (deposited frame, ESR1 305-550) and backbone_atoms.csv.gz.
- Observation mask results/sv/SV00/observed_mask.npy [788, 246]: 181829 observed C-alpha; rows differing from the frozen mapping: 0.

## Chain table (results/sv/SV00/chain_table.csv)

- 788 rows x 53 columns; sources in column_sources.md; missing values per column in missing_by_column.csv.
- ligand_chemotype: {'other': 514, 'amine': 208, 'acid': 37, 'both': 22, 'none': 7}; carboxylic acid on 59 chains, basic amine on 230.
- Dimer partner source: {'asu': 742, 'symmetry': 46}; partner is the chain's own symmetry image: 36; H11 heavy atoms in contact with the partner min 42, median 59; second-best copy above half the best: 0; partner not a kept chain (no cluster): 5T92|B|1, 5T97|A|1.
- Chains with two distinct ligand-of-interest codes (heavy atoms and flags from the primary, first code): 7 (5TLX|A|1, 5TLX|B|1, 6DFN|B|1, 6DFN|D|1, 8EV1|A|1, 8EV2|A|1, 8EV2|B|1).
- Chain-level bound peptide: 419 chains ({'asu': 417, 'none': 369, 'symmetry': 2}); entry-level peptide_present (Y02): 419 chains.
- UniProt-author offset: {0.0: 783, -1000.0: 3, -2000.0: 1, 304.0: 1}; source {'sifts_pdbe_api': 787, 'struct_ref_seq': 1}; constant over the LBD {'True': 777, 'False': 10}.
  Deposited _struct_ref_seq differs from SIFTS for 113 chains in 54 entries (struct_ref_seq offsets {-173.0: 113}); SIFTS equals the existing residue mapping (esr1_position - auth_seq_id) for 787/787 chains.
- multi_valued entries: 2; remediated chains: 0.

Metadata cleaning (metadata_cleaning.csv):

| field | range | missing in mmCIF | out of range -> NaN | filled from API | missing final |
|---|---|---|---|---|---|
| ph | 3.0-11.0 | 109 | 0  | 0 | 109 |
| growth_temp_K | 270.0-315.0 | 35 | 0  | 0 | 35 |
| ambient_temp_K | 50.0-320.0 | 7 | 0  | 0 | 7 |
| r_work |  | 8 | 0  | 0 | 8 |

Value cross-checks, entry level (value_crosscheck.csv): ph vs Y02 ph (also parsed from crystal-growth text): 242/242 agree, only other 65; growth_temp_K vs Y02 temp_K: 316/316 agree, only other 0; ambient_temp_K vs Y02 diffraction_temp_K: 344/344 agree, only other 0; r_work vs Y02 r_work: 343/343 agree, only other 0; resolution_A vs mmCIF _refine.ls_d_res_high: 351/351 agree, only other 0; r_free vs mmCIF _refine.ls_R_factor_R_free: 351/351 agree, only other 0; esr1_chains_asu vs regenerated Q04 table: ESR1 chains in model 1: 351/351 agree, only other 0; release_year vs mmCIF initial release (revision ordinal 1): 351/351 agree, only other 0; space_group vs mmCIF _symmetry.space_group_name_H-M: 351/351 agree, only other 0

Cross-checks (api_crosscheck.csv): ph 242/242 agree; growth_temp_K 316/316 agree; ambient_temp_K 344/344 agree; resolution 351/351 agree; r_free 351/351 agree; r_work 343/343 agree; space_group 351/351 agree; offset_sifts_vs_struct_ref_seq 674/787 agree; offset_sifts_vs_existing_mapping 787/787 agree

Missing values per column (columns with any missing):

| column | missing | fraction |
|---|---|---|
| peptide_bound_auth_asym_id | 369 | 0.4683 |
| peptide_bound_length | 369 | 0.4683 |
| ligand_id_primary | 7 | 0.0089 |
| ligand_heavy_atoms | 7 | 0.0089 |
| ligand_basic_amine | 7 | 0.0089 |
| ligand_carboxylic_acid | 7 | 0.0089 |
| r_work | 9 | 0.0114 |
| ph | 285 | 0.3617 |
| growth_temp_K | 58 | 0.0736 |
| ambient_temp_K | 20 | 0.0254 |
| offset_constant | 1 | 0.0013 |
| multi_valued_fields | 782 | 0.9924 |
| dimer_partner_cluster | 2 | 0.0025 |

## Files

- data/sv/manifest.csv (142.4 kB)
- data/sv/fetch.log (11.1 kB)
- data/sv/lbd_extract (3 files, 26.1 MB)
- data/sv/lbd_chains.csv (5.2 MB)
- data/sv/backbone (2 files, 22.2 MB)
- data/sv/gate (4 files, 10.8 MB)
- data/sv/api_cache (9 files, 2.6 MB)
- data/sv/sifts_cache (351 files, 0.4 MB)
- data/clean (789 files, 122.2 MB)
- results/sv/SV00/chain_table.csv (231.3 kB)
- results/sv/SV00/column_sources.md (4.1 kB)
- results/sv/SV00/missing_by_column.csv (1.5 kB)
- results/sv/SV00/observed_mask.npy (194.0 kB)
- results/sv/SV00/observed_mask_rows.csv (14.1 kB)
- results/sv/SV00/observed_mask_cols.csv (1.9 kB)
- results/sv/SV00/mask_check.json (0.1 kB)
- results/sv/SV00/gate.json (1.0 kB)
- results/sv/SV00/gate_chains.csv (0.2 kB)
- results/sv/SV00/counts.csv (0.3 kB)
- results/sv/SV00/excluded_chains.csv (15.8 kB)
- results/sv/SV00/coordinate_diff_by_chain.csv (27.7 kB)
- results/sv/SV00/rmsd_recomputed.npy (5.0 MB)
- results/sv/SV00/rmsd_recomputed_index.csv (14.1 kB)
- results/sv/SV00/clusters_recomputed.csv (18.1 kB)
- results/sv/SV00/entry_metadata_mmcif.csv (52.7 kB)
- results/sv/SV00/metadata_cleaning.csv (0.3 kB)
- results/sv/SV00/api_crosscheck.csv (0.2 kB)
- results/sv/SV00/chain_offsets.csv (55.8 kB)
- results/sv/SV00/dimer_partner.csv (50.1 kB)
- results/sv/SV00/chain_peptide.csv (54.2 kB)
- results/sv/SV00/extract_compare.csv (0.6 kB)
- results/sv/SV00/value_crosscheck.csv (0.6 kB)
- results/sv/SV00/gate_stdout.txt (1.4 kB)
- results/sv/SV00/env.txt (0.4 kB)
- data/sv/fetch_stdout.txt (0.3 kB)
- data/sv/extract_stdout.txt (0.2 kB)
- data/sv/clean_chain_pdbs_stdout.txt (0.1 kB)
- data/raw/*.cif (526 files added; download_manifest.csv and sha256sums.txt untouched)
- scripts: scripts/sv/SV00_fetch.py, SV00_extract.py, SV00_gate.py, SV00_api_test.py, SV00_metadata.py, SV00_extract_compare.py, SV00_contacts.py, SV00_chain_table.py, SV00_env.py, SV00_summary.py
