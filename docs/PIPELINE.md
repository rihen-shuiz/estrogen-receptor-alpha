# Pipeline (run order)

One command: `bash run_all.sh` (stages `upstream`, `core`, `ext`, `docs`; see the header of `run_all.sh`).
Interpreter: Python 3.11 with `requirements-lock.txt`. Optional external tools: P2Rank 2.5.1 (pocket step X08) and
open-source PyMOL 3.1 (structure images, always through `scripts/pymol_safe.py`). Seed 42 everywhere.
Detailed description of every step: `docs/METHODOLOGY.md`. Task IDs refer to `tracker/esr1.tasks`.

## Stage 1 - upstream (`bash run_all.sh upstream`)

| step | script | produces |
|---|---|---|
| download mmCIF, CCD, data/lbd | `scripts/fetch_data.py` | `data/raw/*.cif`, `data/ccd/`, `data/lbd/` |
| rebuild the upstream tables and compare with the frozen ones | `scripts/rebuild_all.py` (`rebuild_*.py`) | `results/rebuild/*`, `results/rebuild/rebuild_report.csv` |

The frozen upstream tables (`results/lbd_chains.csv`, `results/entries.csv`, `results/ligands.csv`, fingerprints,
K-means, numbering map, ESR2 entries, ChEMBL tables) were produced by `esr_analysis/` modules; the rebuild reproduces
them (identical or equivalent except the stochastic UMAP / t-SNE coordinates).

## Stage 2 - core pipeline (`bash run_all.sh core`)

1. Ligand of interest per chain (Q18) -> cleaned chain PDBs (S17) -> chemistry on ligands of interest (C11).
2. Keep-all labels: attrition funnel (K01), helix table (K02), normalised B per helix (K03), core choice (K04 ->
   `alignment_core_v3.json`), QC tiers (K05).
3. Structure: superposition on core v3 (S03) -> RMSD matrix 788 x 788 (S04/S05) -> H12 metric (S16) -> average-linkage
   clusters at 1.0/1.5/2.0/2.5 A (S08) -> stats, dendrogram, medoids, chain status, QC delta, sensitivity (S09-S14, K06).
4. Redundancy: entropy (R01), discovery curve (R02/R03), saturation (R04), QC comparison (R06), entry level (K07),
   core sensitivity (K09), representative set (S15).
5. ESR1/ESR2 and extensions: displacements (E05-E10), ESR2 matrix / clusters / PCA / basins (X02-X05), potency
   (X15/X16), pockets (X08/X09/X12), PyMOL scenes (F01/F02), provisional September-matrix checks (P01-P04).

## Stage 3 - extended analyses (`bash run_all.sh ext`)

In this order (each script: results in `results/ext/`, figures in `figures/ext/`, a figure card and a numbers file):

`ext_metadata` (Y02) -> `ext_embeddings` (Y03) -> `ext_zoo` (Y04) -> `ext_metric_check` (Y42) -> `ext_compactness`
(Y05) -> `ext_agreement` (Y06) -> `ext_dendrogram` (Y07) -> `ext_chem_zoo` (Y08) -> `ext_ligand_descriptors` (Y09) ->
`ext_displacement_profile` (Y10) -> `ext_structural_features` (Y11) -> `ext_phipsi` (Y12) -> `ext_supervised` (Y13) ->
`ext_split_hypotheses` (Y14) -> `ext_outliers` (Y15) -> `ext_outlier_quality` (Y16) -> `ext_metadata_tests` (Y17) ->
`ext_colour_by` (Y18) -> `ext_interactive` (Y19) -> `ext_coord_pca` (Y20) -> `ext_bfactors` (Y21) -> `ext_esr2` (Y22)
-> `ext_joint_esr1_esr2` (Y23) -> `ext_selectivity` (Y24) -> `ext_chem_struct_link` (Y25) -> `ext_redundancy` (Y26) ->
`ext_time` (Y27) -> `ext_robustness` (Y28) -> `ext_esr2_core_fraction` (Y29) -> `ext_pymol_renders` (Y30) ->
`ext_ncs_copies` (Y43) -> `ext_h12_contacts` (Y44) -> `ext_cl4_substructure` (Y45).

Shared helpers: `scripts/ext_common.py` (paths, loaders, palette, figure saving, card registry, BH, Cramer's V),
`scripts/pymol_safe.py` (PyMOL runner).

## Stage 4 - documents (`bash run_all.sh docs`)

`scripts/ext_build_index.py` -> `docs/figures/INDEX.md`, `docs/gallery.html`;
`scripts/collect_numbers.py` (+ `collect_numbers_ext.py`) -> `results/numbers.json`, `results/ext/numbers_ext.json`,
and `docs/METHODOLOGY.md`, `docs/RESULTS_SUMMARY.md` rendered from `docs/templates/`.

## Checks

- `python tracker/tl.py validate` and `python tracker/tl.py status` - task state; `tl.py check ID` re-runs a task's checks.
- `python scripts/check_evidence.py` (tasks up to 2026-10-20) and `python scripts/check_evidence_ext.py`
  (final-run tasks) - independent recounts of the logged claims.
