# Curated ligand action classes: `esr1_dataset_final.csv`

One row per PDB entry × ligand (CCD code): **377** rows, **368** with `is_real_ligand == True`.
`action_class` counts: antagonist 119, agonist 104, partial_agonist 79, unknown 64, non_ligand 8, mixed 2,
not_applicable 1. The labels were assigned by hand from the primary literature of each entry.

## Columns
- `pdb_id`: PDB entry code.
- `ligand_id`: PDB chemical component (CCD) code of the ligand.
- `doi`: DOI of the primary article (empty when the entry has none).
- `action_class`: functional class at ERα (`agonist`, `partial_agonist`, `antagonist`, `mixed`, `unknown`,
  `non_ligand`, `not_applicable`).
- `drug_classes`: pharmacological classes (`SERM`, `SERD`, `SERCA`, `CERAN`, `DMERI`), separated by `|`.
- `compound_id`: compound name in the article or library.
- `mutations`: ERα LBD mutations as written in the source (empty for wild type).
- `mutations_norm`: normalised mutations (`WT`, `YES (unspecified)` or a list such as `A123B` sorted by residue).
- `mutation_group`: coarse category for plots (`WT`, `Y537S`, `L372S, L536S`, `C381S, C417S, C530S`,
  `C381S, C417S, C530S, L536S`, `YES (unspecified)`, `Other`).
- `is_mutant`: `True` if the structure carries mutations.
- `is_real_ligand`: `False` for `non_ligand` / `not_applicable` (ions, amino acids, apo entries).
- `has_peptide`: whether a coactivator or corepressor peptide is present.
- `inchikey`: InChIKey of the ligand.
- `evidence_sources`: every source of evidence, separated by `|`.
- `primary_source`: the main source by priority (`full_text` > `abstract` > `derived_from_data` > `chembl` >
  `ligand_identity` > `literature` > `pdb_metadata` > `none`; entries recomputed from 10.1038/nchembio.1214 are
  `derived_from_data`).
- `luc_efficacy_pct`: ERα transcriptional efficacy (% of E2; MCF-7 luciferase or L-Luc ERα-WT).
- `proliferation_pct`: cell proliferation (% of E2; MCF-7 proliferation).
- `pdb_title`: PDB entry title.
- `note`: remarks, including mismatches between a mutant structure and wild-type assay data, and recomputed rows.

## Use in plots
1. Keep `is_real_ligand == True` to drop ions, amino acids and apo entries.
2. Nine rows carry `SERD` in `drug_classes` but `action_class == "unknown"` (degradation shown without a separate
   antagonism or antiproliferation test).
3. `derived_from_data` labels follow fixed thresholds: efficacy ≥ 80% = agonist, 20–80% = partial_agonist,
   < 20% with proliferation below 0% = antagonist.
