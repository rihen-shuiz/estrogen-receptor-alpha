# Repository map

What matters for reading and re-running the project, and what is archive. Start with `docs/RESULTS_SUMMARY.md`,
`docs/gallery.html` and `docs/METHODOLOGY.md`.

## What matters

| path | what it is |
|---|---|
| `run_all.sh`, `requirements-lock.txt` | one command to reproduce everything; exact package versions |
| `tracker/` | `esr1.tasks` (spec of every task, exact-path checks), `esr1.log` (what was done, evidence), `tl.py` (the only tool that writes the log) |
| `scripts/` | every analysis step: core tasks (`align_selected_chains.py`, `build_selected_rmsd.py`, `cluster_selected.py`, ...), upstream rebuild (`rebuild_*.py`), final-run analyses (`ext_*.py`), helpers (`ext_common.py`, `pymol_safe.py`), numbers (`collect_numbers*.py`), evidence checks (`check_evidence*.py`) |
| `esr_analysis/` | the package with the scientific functions (structure reading, LBD extraction, mapping, alignment, RMSD, clustering, PyMOL scenes) and `resources/config/` (all thresholds: `analysis_settings.json`, `alignment_core_v3.json`, `qc_tiers.json`) |
| `data/` | inputs: `pdb_ids_esr1.txt`, `pdb_ids_esr2.txt`, `raw/` (mmCIF, git-ignored, sha256 tracked), `clean/` (one PDB per chain), `chembl_esr1_activities.csv`, `selectivity_residues.csv` (draft, E01), `cache/` (RCSB API JSON, ESR2 mmCIFs) |
| `results/` (top level) | frozen outputs of the core tasks: `rmsd_entry.npy` + `rmsd_entry_index.csv` (the 788 x 788 matrix), `entry_clusters.csv`, `chain_status.csv`, `lbd_chains.csv`, `qc_tiers.csv`, `h12_metric.csv`, `helix_bfactor_by_chain.csv`, `ligands_of_interest.csv`, `rmsd_esr2.npy`, ... |
| `results/alignment/` | superposed C-alpha coordinates (`aligned_coordinates.npz`) and transforms |
| `results/ext/` | every result table of the final run (`zoo_runs.csv`, `structural_features.csv`, `split_hypotheses.csv`, `outlier_flags.csv`, `entry_metadata.csv`, `rmsd_joint.npy`, ...), `figure_cards.csv`, `numbers_ext.json` |
| `results/rebuild/` | upstream tables rebuilt from raw inputs and `rebuild_report.csv` (frozen vs rebuilt) |
| `figures/` | figures of the core tasks; `figures/ext/` = final-run figures (PNG + SVG), `figures/ext/pymol_panels/` raw PyMOL renders |
| `docs/` | `RESULTS_SUMMARY.md`, `METHODOLOGY.md`, `PIPELINE.md`, `PLAN_EXTENDED.md`, `EXPLORATION.md`, `HUMAN_STEPS.md`, `gallery.html`, `interactive_embedding.html`, `figures/INDEX.md`, `numbers_ext.md`, `EVIDENCE.md` / `EVIDENCE_EXT.md`, decision drafts (`chemistry_k_choice.md`, `decision_weeks9_12.md`), `templates/` (sources of the generated documents), `manuscript/` |
| `notebooks/worked_example.ipynb` | a small worked example on saved outputs |

## Archive (kept, not part of the current pipeline)

| path | what it is |
|---|---|
| `results/archive_v2_strict/` | the earlier strict one-chain-per-entry run (351 chains), kept for comparison (K00) |
| `results/provisional/`, `figures/provisional_*` | checks on the September matrix before S04 was redone (P01-P04) |
| `results/task_outputs/`, `task_outputs/`, `results/figures/`, `results/reference_render/`, `reference_render/`, `recipes/` | the September report package (blocks 01-13 and its eleven figures) |
| `verification/`, `VERIFICATION.txt`, `BLOCK_IO.txt`, `PACKAGE_MANIFEST.json`, `INPUT_MANIFEST.csv`, `FIGURE_MANIFEST.csv`, `README.txt`, `ESR1_Function_Blocks_and_Monthly_Checklist.txt`, `ESR1_September_Report.pdf` | September package documentation and verification records |
| `evidence_pack_20261020.zip`, `tracker/bundle_*.zip` | earlier review packs (see W08 in `docs/HUMAN_STEPS.md`) |
| `tracker/esr1.tasks.bak*`, `update_tasks_keepdata.py` | backups and the one-off script of the keep-all task revision |
| `tests/test_reported_contracts.py` | tests of the September package contracts |
