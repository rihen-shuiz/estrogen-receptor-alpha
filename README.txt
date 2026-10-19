ESR1 FIRST-MONTH STANDALONE

This package isolates thirteen scientific task files, their required helpers,
frozen inputs and eleven local-data report figures. It excludes fingerprints,
chemical clustering/UMAP, medoids, matched receptor pairs, selectivity, entropy,
QC reclustering, sensitivity, ChEMBL activities, pockets and later-report code.
Existing exact-compound action annotations are evidence for block 07; they are
not new ChEMBL acquisition or later expanded manual curation.

INSTALL

Use Python 3.12 for the tested constraints (NumPy and SciPy require >=3.12).
Source syntax supports >=3.10 with compatible older dependency versions, but
those combinations were not tested. The actual clean test used Linux x86_64
and Python 3.12.3; see VERIFICATION.txt.
From the extracted ESR1_First_Month_Standalone/ directory:
  python -m venv .venv
  .venv/bin/python -m pip install -c constraints.txt .
For the included focused checks, install pytest separately:
  .venv/bin/python -m pip install -c constraints.txt pytest
  .venv/bin/python -m pytest -q tests
The tested versions were actually downloaded and installed in a fresh temporary
environment. Core requirements are minimal; constraints.txt records that tested
environment rather than claiming all recorded tools are required dependencies.
No umap-learn is retained. The wheel includes its resources, so default input
paths remain valid after installation and when commands run from another folder.

DEFAULT OFFLINE REPORT ROUTE

  .venv/bin/python -m esr_analysis.report_figures --output RUN_OUTPUT/figures

The command reads only installed frozen resources. It makes eleven PNG/PDF
figures and truthful TXT captions without live queries, downloads or a full
structural fit. FIGURE_MANIFEST.csv maps each panel to blocks and exact files.
The reference panel uses the full genuine saved 7NFB image. A separately verified
new PyMOL render is included; it is not a cached-image-only render claim.

The named seven-page PDF, LaTeX source ZIP and report figure recipe were not
present in the workspace or nested ZIPs. The user instructed proceeding with
local contents. Numerical reproduction is verified against frozen inputs and
the explicit supplied acceptance values. Exact missing-report colours, limits,
crop/layout and appearance cannot be certified; visual reference matching is
explicitly skipped. Available local scene appearance and captions are preserved.

TASK EXECUTION — ONE REAL COMMAND PER BLOCK

Each module supports --help, --input-root, --source-folder, --settings and an
explicit --output. Defaults resolve from installed package resources, never an
author's home directory or the shell's working directory. Relative explicit
output arguments resolve normally from the invoking shell. These commands were
executed, not merely checked for importability:
  .venv/bin/python -m esr_analysis.query_structures --output RUN_OUTPUT/block_01
  .venv/bin/python -m esr_analysis.download_structures --output RUN_OUTPUT/block_02
  .venv/bin/python -m esr_analysis.extract_metadata --output RUN_OUTPUT/block_03
  .venv/bin/python -m esr_analysis.extract_lbd --output RUN_OUTPUT/block_04
  .venv/bin/python -m esr_analysis.flag_quality --output RUN_OUTPUT/block_05
  .venv/bin/python -m esr_analysis.build_ligand_table --output RUN_OUTPUT/block_06
  .venv/bin/python -m esr_analysis.assign_action_classes --output RUN_OUTPUT/block_07
  .venv/bin/python -m esr_analysis.choose_reference --output RUN_OUTPUT/block_08
  .venv/bin/python -m esr_analysis.align_structures --output RUN_OUTPUT/block_09
  .venv/bin/python -m esr_analysis.calculate_rmsd --output RUN_OUTPUT/block_10
  .venv/bin/python -m esr_analysis.cluster_structures --output RUN_OUTPUT/block_11
  .venv/bin/python -m esr_analysis.structural_pca --output RUN_OUTPUT/block_12
  .venv/bin/python -m esr_analysis.regional_variation --output RUN_OUTPUT/block_13

The commands default to bundled upstream prerequisites, so one block does not
silently execute another. To consume a new upstream output, pass its folder with
--source-folder and supply the other documented prerequisites there. BLOCK_IO.txt
lists signatures, method/tool, inputs and outputs. A missing file is reported
with its prerequisite path; no implicit download or baseline reconstruction runs.

RAW AND LIVE ROUTES

The promised raw pilot is exactly 3Q95, 7NFB, 7BAA and 6CHZ, with all deposited
chains/models retained. Four actual hashed mmCIFs are bundled. Blocks 02/03/04
run directly on them. Example selecting one real raw entry:
  .venv/bin/python -m esr_analysis.extract_lbd --ids 7NFB --output RUN_OUTPUT/7NFB_domains
All 400 required local ESR1 CCDs and action annotation resources are bundled,
so blocks 06/07 reproduce the whole report identity population offline.
Full 487-entry raw re-inventory is not promised: the other raw mmCIFs are not
repacked. To run more entries, explicitly supply their local --raw-folder or
download exact IDs into a separate directory:
  .venv/bin/python -m esr_analysis.download_structures --live --ids 7NFB --output RUN_OUTPUT/live_download
  .venv/bin/python -m esr_analysis.query_structures --live --output RUN_OUTPUT/live_query
Network collection was not needed or executed for report reproduction. Live
query files are labelled separately and never replace the frozen report cohort.
The frozen ID snapshot is derived from the supplied ID list with its source hash;
it is not falsely presented as an original cached HTTP response or retrieval date.

PYMOL — SEPARATE RENDERING APPLICATION

requirements.txt does not install PyMOL. Install an appropriate PyMOL executable
separately. The tested temporary environment installed pymol-open-source 3.2.0a0:
  .venv/bin/python -m pip install -c constraints.txt pymol-open-source
  .venv/bin/python -m esr_analysis.render_reference --pymol /absolute/path/to/pymol --output RUN_OUTPUT/reference_render
To use that newly rendered full image in the report figures:
  .venv/bin/python -m esr_analysis.report_figures --reference-image RUN_OUTPUT/reference_render/reference_context.png --output RUN_OUTPUT/rendered_figures
recipes/reference_context.pml is relocatable. Direct PyMOL execution requires
explicit ESR1_REFERENCE_CIF and ESR1_REFERENCE_OUTPUT environment variables;
the Python render command supplies them. The recipe applies the original saved
whole-object transform and checks every transformed coordinate. It retains the
actual ASU: 7NFB B blue, other polymers grey, peptides orange, GEN green across
copies. No LaTeX crop is applied to the full molecular image. The original PML
in resource provenance retains its producing-workspace paths and is not executed.

BLOCK 13 AND OPTIONAL MDANALYSIS

Population RMSF uses available experimental observations after the saved core
fit; missing atoms remain NaN. This is not an MD simulation. Within-chain CA B z
uses ddof=0, B>=0 and >=80% valid regional coverage. H12* 536–545 remains a
provisional annotation, distinct from configured H12 538–547.
MDAnalysis is optional for a separate common-atom numerical check. It was absent
in the clean core environment; that check is reported UNAVAILABLE, not passed.
The NumPy calculation and its reproduction were verified. To request the optional
cross-check, install MDAnalysis separately and rerun block 13 into a new output.
The original frozen common-atom trajectory is historical evidence, not proof of
a new MDAnalysis execution. LaTeX is optional and is not needed for analysis or
ordinary PNG/vector-PDF generation.

INPUTS, RESULTS AND EVIDENCE

The full master remains 487 entries; 351 selected chains and 400 ESR1 identities
are explicit subsets. Some fields retain original absolute source paths solely
as provenance; runtime paths use declared resources or arguments. Raw source
paths in copied entry records are portable, with original paths retained separately.
Provisional and unknown action/annotation statuses remain explicit. No later
manual action expansion or new biological claim is introduced.
INPUT_MANIFEST.csv records every package input's purpose and SHA-256.
verification/SOURCE_PROVENANCE.json records original source hashes and explicit
row/column selections. Figures and verification records describe the actual
executed clean-package commands. The original workspace and scientific results
were not modified. See VERIFICATION.txt for numerical tolerances and limitations.

DELIVERED EXECUTION RECORDS

figures/ contains all eleven regenerated PNG/PDF figures and TXT captions.
reference_render/ contains the actual independently reproduced 7NFB render.
task_outputs/block_01 through block_13 contain the outputs of the verified real
commands, including explicitly copied prerequisites for block 13. These are
execution records; default reproducible inputs remain the installed resources.
The monthly checklist groups the thirteen tasks, one checked main file each.
verification/run_tasks.py and compare_outputs.py retain the exact temporary
verification harnesses; their /tmp paths record this run, not package runtime
requirements. Use the portable README task commands for a different extraction.
