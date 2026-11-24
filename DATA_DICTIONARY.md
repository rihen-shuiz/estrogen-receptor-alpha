# ESR1 LBD Dataset — Data Dictionary (`esr1_dataset_final.csv`)

Total records: **377** (`is_real_ligand == True`: **368**).
Distribution of `action_class`: `antagonist`: 119, `agonist`: 104, `partial_agonist`: 79, `unknown`: 64, `non_ligand`: 8, `mixed`: 2, `not_applicable`: 1.

## CSV Columns
- `pdb_id`: PDB structure code.
- `ligand_id`: 3-letter (or extended) PDB chemical component (CCD) ligand code.
- `doi`: DOI of the primary publication (empty for entries without a DOI).
- `action_class`: functional class on ERα (`antagonist`, `agonist`, `partial_agonist`, `mixed`, `unknown`, `non_ligand`, `not_applicable`).
- `drug_classes`: pharmacological classes (`SERM`, `SERD`, `SERCA`, `CERAN`, `DMERI`), pipe-delimited (`|`).
- `compound_id`: compound identifier in the source publication/library.
- `mutations`: raw ERα LBD construct mutation string (empty for WT).
- `mutations_norm`: normalized mutations (`WT`, `YES (unspecified)`, or residue-number-sorted `A123B` list).
- `mutation_group`: coarse mutation category for plotting (`WT`, `Y537S`, `L372S, L536S`, `C381S, C417S, C530S`, `C381S, C417S, C530S, L536S`, `YES (unspecified)`, `Other`).
- `is_mutant`: `True` if the structure carries mutations, `False` otherwise.
- `is_real_ligand`: `False` for `non_ligand`/`not_applicable` (ions, amino acids, apo structures), `True` otherwise.
- `has_peptide`: boolean flag (`True`/`False`) indicating presence of a coactivator/corepressor peptide.
- `inchikey`: ligand InChIKey.
- `evidence_sources`: all supporting evidence sources, pipe-delimited (`|`).
- `primary_source`: highest-priority evidence source (`full_text` > `abstract` > `derived_from_data` > `chembl` > `ligand_identity` > `literature` > `pdb_metadata` > `none`; recalculated entries from `10.1038/nchembio.1214` use `derived_from_data`).
- `luc_efficacy_pct`: ERα transcriptional efficacy (% of E2, from `MCF-7 Luc` or `L-Luc ERα-WT`).
- `proliferation_pct`: cell proliferation (% of E2, from `MCF-7 cell proliferation` or `Proliferation`).
- `pdb_title`: PDB entry title.
- `note`: curation notes (including mutant-structure vs WT-assay mismatches and recalculation notes).

## Plotting and Filtering Notes
1. Filter `is_real_ligand == True` to exclude non-ligands (ions, amino acids, apo structures).
2. Nine entries carry `SERD` in `drug_classes` while `action_class == "unknown"` (degradation demonstrated without a separate antagonism/anti-proliferation assay).
3. `derived_from_data` thresholds: `>=80%` = `agonist`, `20-80%` = `partial_agonist`, `<20%` with proliferation inhibition (`<0%`) = `antagonist`.
