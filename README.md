# ESR1 September analysis

Implementation and saved results for blocks 01–13 in
[ESR1_September_Report.pdf](ESR1_September_Report.pdf). These blocks cover
structure collection, receptor mapping, quality flags, exact ligand identities,
supported action labels, reference selection, alignment, RMSD, structural
groups, PCA, and regional B-factor and coordinate variation summaries.
Later checklist tasks are outside this repository's implemented scope.

## Contents

- `esr_analysis/`: analysis modules, plotting and rendering helpers, and bundled inputs.
- `results/task_outputs/block_01` through `block_13`: executed offline task results.
- `results/figures/`: eleven PNG/PDF figures with text captions.
- `results/reference_render/`: independently rendered 7NFB reference image and scene.
- `requirements.txt`: core Python dependencies.

## Run

Use Python 3.12. From this repository's root:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m esr_analysis.report_figures --output RUN_OUTPUT/figures
```

Each task accepts an explicit output directory. For example:

```bash
.venv/bin/python -m esr_analysis.extract_lbd --output RUN_OUTPUT/block_04
.venv/bin/python -m esr_analysis.calculate_rmsd --output RUN_OUTPUT/block_10
.venv/bin/python -m esr_analysis.regional_variation --output RUN_OUTPUT/block_13
```

Task modules, in report order:

```text
query_structures
download_structures
extract_metadata
extract_lbd
flag_quality
build_ligand_table
assign_action_classes
choose_reference
align_structures
calculate_rmsd
cluster_structures
structural_pca
regional_variation
```

Use `python -m esr_analysis.<module> --help` for options. Commands default to
bundled frozen prerequisites; running these commands in sequence does not
automatically connect their outputs. `--source-folder` selects explicit upstream
inputs, which must contain all prerequisites for that task. The figure command
also reads frozen products by default.

Only four raw mmCIFs are bundled: 3Q95, 7NFB, 7BAA, and 6CHZ. Full raw
reconstruction requires the remaining structures; `python scripts/fetch_data.py`
downloads them into `data/raw`, fetches missing CCD files into `data/ccd`, and
rebuilds `data/lbd`. Block 09 defaults to two
chains; the full saved ensemble contains 788 aligned chains. Block 06 also
resolves action evidence. Block 12 uses saved PCA signs and structural labels.
Block 13 recomputes profiles and displacements, but its regional producer reads
profiles and displacements from the selected source folder. These commands
reproduce saved analyses, not a complete automated raw-to-report pipeline.

Live queries and downloads require an explicit `--live` option and create
separate outputs. They were not used for the saved offline results.

PyMOL is a separate application. To reproduce the reference image:

```bash
.venv/bin/python -m esr_analysis.render_reference \
  --pymol /absolute/path/to/pymol --output RUN_OUTPUT/reference_render
```

MDAnalysis is optional for a common-atom RMSF cross-check. Neither PyMOL nor
MDAnalysis is installed by `requirements.txt`. Their absence does not prevent
the core offline analysis or figure generation.

## Reproduced results

The 2026-10-04 run used Python 3.12.3. Core numerical results matched the report:
487 master entries, 444 verified ESR1 entries, 901 ESR1 chain-models, 788 eligible
copies, 351 selected entries, and 247 selected entries passing combined QC.
The 400 exact ligand identities retain 3 supported agonist labels, 4 supported
degrader labels, and 393 unknown labels. The reference is `7NFB|B|1`.

The 351-entry matrix contains 61,425 distinct pairs. Average linkage at 1.5 Å
produces nine groups with counts 2, 2, 2, 124, 1, 6, 212, 1, 1. PCA fractions
are 0.9061001199700581 and 0.07493449329557664. Reported H6, H12*, H5, and H8
regional values reproduced.

All thirteen task commands, their CLI help, eleven figure producers, and the
eleven existing focused tests passed. Thirty-three baseline numerical comparisons
passed. A separate connected recomputation from saved unaligned mappings refitted
all 788 chains and reproduced downstream distances, groups, PCA, ensemble
profiles, displacements, and regional values. This does not establish full raw
collection or selection reproduction. MDAnalysis 2.10.0 passed the optional
cross-check; actual PyMOL rendering also passed.

Ten regenerated PNGs matched the previously supplied PNGs exactly. The alignment
figure differs in overall title and layout; its underlying coordinate stages
reproduce. The supplied report PDF is retained, not rebuilt or certified for
exact layout reproduction. Tests and verification harnesses are intentionally
excluded from the GitHub upload.

Coordinate diversity, B-factors, and contact annotations are descriptive
measurements. They do not establish biological causality. H12* 536–545 is a
provisional regional annotation, distinct from configured H12 538–547.

## Known limitations

- Coordinates are C-alpha only. No all-atom RMSD validation is in this repository yet.
- The reference `7NFB|B|1` is an engineered PROSS variant with 28 declared substitutions, not wild-type ESR1. See `docs/reference.md`.
- The ligand table includes crystallization ions and additives (for example MG, CL, GOL, EDO). The top-ranked "ligands" are mostly these. See `results/ligands_ranked.csv`.
- Chemistry clusters have low silhouette (0.10 to 0.13 for K=3 to 7). The choice of K is provisional until the human sign-off in `docs/chemistry_k_choice.md`.
- The ESR2 validation set is 39 entries. It is kept separate from the 351-entry ESR1 selection.
- Cross-receptor ligand matches (`results/matched_pairs.csv`) use identical ligand codes only, with at least 10 heavy atoms.
- Not yet in the repository: P2Rank pocket analysis and the manuscript text. The ChEMBL activity table is a separate pending task.

## Sources

Scientific code and resources originate in the existing ESR1 project. Source
and evidence fields remain in bundled resources. RCSB PDB mmCIFs, canonical
UniProt sequences, CCD descriptors, and receptor-specific action annotations
support the reported analysis. Dependency licenses remain with their respective
distributions. No new license is applied to the project's code or data.
