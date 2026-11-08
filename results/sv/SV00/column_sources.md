# SV00 chain_table.csv: column sources

788 rows (kept chains, matrix order), 53 columns. Missing values per column: results/sv/SV00/missing_by_column.csv.

| column(s) | source |
|---|---|
| chain_key, entry_id | results/rmsd_entry_index.csv (S05), matrix row order |
| auth_asym_id | results/lbd_chains.csv (Q04) |
| s08_cluster | results/entry_clusters.csv cut_1.5 (S08, 1.5 A average linkage) |
| qc_label, tier_resolution, tier_rfree, tier_b | results/qc_tiers.csv (K05) |
| h12_label, n_h12_ca | results/h12_metric.csv (S16) |
| peptide_present_entry, peptide_class_entry | results/ext/chain_metadata.csv (Y02; entry level: any non-receptor polypeptide of <= peptide_max_length = 50 residues in the entry) |
| peptide_bound, peptide_bound_source, peptide_bound_auth_asym_id, peptide_bound_length | SV00 results/sv/SV00/chain_peptide.csv (SV00_contacts.py): peptide copy (same definition, <= 50 residues) with the most receptor heavy atoms within 4.0 A, ASU or symmetry copy |
| construct_mutations, mutation_status | results/ext/chain_metadata.csv mutations, mutation_status (Y02, from the residue-level mapping of results/lbd_chains.csv) |
| ligand_codes, n_ligand_codes, ligand_status | results/chain_ligand_of_interest.csv (Q18); codes = distinct ligand_id per chain; status of the first row |
| ligand_id_primary | first ligand-of-interest row per chain (= results/ext/structural_features.csv ligand_id, the convention of scripts/ext_common.chain_table) |
| ligand_heavy_atoms, ligand_basic_amine | results/ext/structural_features.csv (Y11) from results/ext/ligand_descriptors.csv (Y09): RDKit heavy atoms; basic amine = SMARTS [NX3;H2,H1,H0;!$(NC=O);!$(NS=O);!$(N-a);!$(N-C=N);!$(N#*);!$(N=*)]([C;!$(C=*)]) (existing definition) |
| ligand_carboxylic_acid | SV00: RDKit SMARTS [CX3](=O)[OX1H0-,OX2H1] on results/ligands_of_interest.csv canonical_isomeric_smiles of ligand_id_primary (acid or carboxylate only) |
| ligand_chemotype | SV00: amine / acid / both / other / none from ligand_basic_amine and ligand_carboxylic_acid; none = no ligand of interest |
| resolution_A, r_free, mean_b_entry_A2 | results/qc_tiers.csv (K05; values of results/entries.csv, Q05) |
| r_work | SV00 results/sv/SV00/entry_metadata_mmcif.csv: mmCIF _refine.ls_R_factor_R_work (first row); Data API only where the mmCIF has none (r_work_source) |
| mean_ca_b_chain_A2 | results/lbd_chains.csv mean_CA_B_A2 (Q04, chain mean C-alpha B) |
| ph, growth_temp_K | SV00 entry_metadata_mmcif.csv: _exptl_crystal_grow.pH / .temp, crystal_id 1, '?'/'.' -> NaN, [DEFAULT] plausible pH 3-11, 270-315 K, else NaN; Data API fill only where missing |
| ambient_temp_K | SV00 entry_metadata_mmcif.csv: _diffrn.ambient_temp, diffrn id 1, [DEFAULT] 50-320 K |
| space_group | results/entries.csv space_group (Q05) |
| esr1_chains_asu, release_year | results/ext/entry_metadata.csv chains_asu, release_year (Y02) |
| tls_refined | SV00 entry_metadata_mmcif.csv: category _pdbx_refine_tls exists in the mmCIF |
| uniprot_author_offset, offset_constant, offset_source | SV00 results/sv/SV00/chain_offsets.csv: SIFTS (PDBe API mappings/uniprot, P03372 segments, UniProt 305-550) minus author number (_pdbx_poly_seq_scheme); offset = most common per-residue value, offset_constant = one value over the LBD |
| multi_valued, multi_valued_fields | SV00 entry_metadata_mmcif.csv: several rows in _exptl_crystal_grow, _diffrn or _refine (id 1 / first row used) |
| latest_revision_date | data/sv/manifest.csv (_pdbx_audit_revision_history of the restored file) |
| remediated | SV00 results/sv/SV00/gate.json remediated_chains (only filled when the gate is 'tolerated') |
| dimer_partner_* | SV00 results/sv/SV00/dimer_partner.csv (SV00_contacts.py): receptor copy with the most H11 (497-530) heavy atoms within 4.0 A; source asu / symmetry / none; self_image = the chain's own symmetry image |
| dimer_partner_peptide_bound | peptide_bound of the partner chain (chain_peptide.csv) |
| dimer_partner_cluster | results/entry_clusters.csv cut_1.5 of the partner chain (blank if the partner is not a kept chain) |
