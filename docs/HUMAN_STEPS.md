# Remaining human steps (final run, 2026-10-25)

Everything a machine can do is done (see `python tracker/tl.py status`). The steps below need a person. For each:
what is ready, and the exact action left. Log with `python tracker/tl.py log <ID> DONE --who human --ev <files>`.

| ID | owner | what is ready | action left |
|---|---|---|---|
| Q10, Q11 | teammate | nothing: these labels are made by hand by a teammate; by project rule no labels were searched, drafted or pre-filled | write `results/ligand_action_labels.csv` (cols ligand_id, action, source), log Q10 then Q11 |
| C08 | chemistry lead | `docs/chemistry_k_choice.md` (section "Auto-selection"), `results/ext/chem_zoo_ranking.csv`, `results/ext/chem_labels_best.csv` | accept the auto-selected K = 6 (MACCS, k-means on MDS) or pick another run; log C08 DONE --who human |
| E01 | D | `data/selectivity_residues.csv` has 2 draft rows (L384M, M421I; Pike et al. 1999). `results/ext/receptor_residue_diff.csv` lists every ESR1/ESR2 sequence difference with its backbone shift | confirm the citation; add >= 2 more positions with real references (the card needs >= 4 rows); log E01 |
| W01 | team | `docs/decision_weeks9_12.md` (draft, recommendation: implement) | choose an option and sign the Decision line |
| W05 | team | all work is in the working tree (nothing committed by Claude) | review and commit to the chosen branch; push |
| W08 | team | `evidence_pack_20261020.zip`, `tracker/bundle_*.zip`, `final_pack_20261025.zip` are the packs on disk | keep the final pack, move older packs/duplicates to `archive/`, confirm |
| W17 | team | storyline of 10 figures in `docs/RESULTS_SUMMARY.md`, gallery `docs/gallery.html` | present to the supervisor, write `docs/supervisor_feedback.md` (>= 5 lines), log W17 |

## What unblocks what

- **Q10 (labels)** unblocks S11 (PCA by ligand class), R05 (diversity per class), C09, L01, X10, X11, then T01-T03,
  X06, R07 and W11. The extended analyses already test the agonist/antagonist question indirectly (peptide, construct,
  ligand basic amine, H12 geometry); with labels, `scripts/ext_split_hypotheses.py` can add the direct test.
- **C08** unblocks L02-L06 and W12. The ligand-level chemistry x structure test is already done on the auto-selected
  clusters (Y25), so L02/L03 can reuse `results/ext/chem_struct_link.csv` once K is signed.
- **E01** unblocks E02 -> E09 -> E11 -> W13, and F03 -> F04 -> F05.

## Notes for the reviewer

- `W04` evidence path `data/lbd/entries.csv` is a git-ignored rebuild product that is not on disk now;
  `python scripts/fetch_data.py` recreates it. The W04 check itself passes.
- `data/raw` now holds all 526 mmCIFs (the 175 missing ones were fetched by `scripts/rebuild_fetch_raw.py`; every
  sha256 matches `data/raw/sha256sums.txt`).
- Corrections to earlier DONE notes are logged as NOTE lines in `tracker/esr1.log` (Y10, Y15, Y17, Y27).
