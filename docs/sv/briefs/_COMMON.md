# Common instructions for every RUN_SV1028 step subagent

You run ONE step of docs/sv/RUN_SV1028.md for the ESR1 project at /Users/ayan/personal/estrogen-receptor-alpha
(git branch `sv1028-run`). Your step text, §1 rules and §2 order (and §3 data access when needed) are in
`docs/sv/briefs/<ID>.md`, extracted verbatim from the run file. Follow them exactly. The full run file is
docs/sv/RUN_SV1028.md if you need context, but do only your step.

Operating mode (from §0 of the run file):
- No questions to the user. Apply [DEFAULT] rules as written. If a rule cannot be applied or gives a degenerate
  result, finish what you can, mark the affected part BLOCKED in your SUMMARY with the reason and what was tried.
- Scripts go to `scripts/sv/<ID>_*.py`, shared code to `scripts/sv/common/`. Every number in a verdict or summary
  must come from a saved file. Write `results/sv/<ID>/SUMMARY.md` (rule 1.7 format) and `results/sv/<ID>/env.txt`.
- Read-only: existing files in results/, figures/, data/clean/, esr_analysis/. Write only where rule 1.2 allows.
- Find inputs from earlier work through the `out` paths of the cards in tracker/esr1.tasks (rule 1.10); do not guess.

Environment on this machine:
- Python: ALWAYS `/Users/ayan/personal/estrogen-receptor-alpha/.venv/bin/python` (3.11, packages from
  requirements-lock.txt + pydssp). See docs/sv/env_pre.txt. System python3 is 3.9 without packages.
- Do NOT install packages (pip/uv installs are blocked on this machine). If something is missing, use the fallback
  your step names and record it in line 1 of the SUMMARY; if there is none, mark that part BLOCKED.
- Not available: mkdssp (use pydssp fallback), conda, Homebrew, rsvg-convert; cairosvg is installed but libcairo is
  missing, so it cannot render.
- Network: RCSB (files/data) works. Be polite (<= 5 req/s, retries with backoff) and cache responses as §3 says.
- Long commands: run them in the background or in chunks; a single foreground shell call times out after 10 min.

Project rules you must keep:
- Do not run git commit/push, do not delete existing files, do not edit tracker/esr1.log by hand.
  Do not log your card in the tracker; the orchestrator does that after checking your outputs.
- CLAUDE.md: Q10 and Q11 (ligand action labels) are done by hand by a teammate. Do not search for labels for them,
  draft them, pre-fill results/ligand_action_labels.csv, or log MANUAL for them.
- Do not recompute a saved matrix unless your step says so (SV00 gate and SV05 do).

When done, reply with the full text of results/sv/<ID>/SUMMARY.md, the list of outputs with their sizes, and the
status you propose for RUN_STATUS (DONE / BLOCKED / partially BLOCKED with what).

Notes from earlier steps:
- Run Python with `-B` (e.g. `.venv/bin/python -B script.py`) so no __pycache__ is written into esr_analysis/ or scripts/.
- Shared figure helper: scripts/sv/common/sv_fig.py (project palette; saves PNG + SVG and one cropped panel per axes
  into figures/panels/sv/<ID>/, rule 1.8). Reuse it.
- SV00 (DONE): results/sv/SV00/chain_table.csv (788 kept chains, column sources in column_sources.md), observed_mask.npy
  (+ _rows/_cols csv), backbone of every kept chain in data/sv/backbone/backbone.npz (deposited frame, UniProt 305-550),
  core-v3 transforms of the gate fit in data/sv/gate/alignment_transforms.json, restored mmCIFs in data/raw/, cleaned
  chain PDBs in data/clean/. No chain is flagged remediated.
- Rule 1.1 vs 1.2: ideas outside your step may be APPENDED to docs/EXPLORATION.md under a heading
  "Suggested, not done" (rule 1.1 names that file explicitly); never edit its existing text.
- results/medoids.csv (S12): its cl4 medoid 2QE4|B|1 has no modelled H12. Where a step needs complete-H12 medoids,
  the run file names 5TM1|A|1 (cl7) and 7UJ7|C|1 (cl4) (the Y10 pair; SV06 used them).
- SV02 (DONE, pydssp fallback, 'H' = any helix): consensus helices results/sv/SV02/helix_boundaries.csv (H1 312-321,
  H3 342-361, H4/5 372-394, H7 422-437, H8 442-455, H9 473-492, H10/11 497-526, H12 538-544), per-chain table
  sse_per_chain.csv, phi/psi phipsi_per_chain.csv. Data-derived H12 window 538-544 is a PROPOSAL; primary analyses keep 538-547.
- SV01a (DONE): results/sv/SV01a/ (model code scripts/sv/SV01a_completeness.py, model_data.csv, block definitions).
- SV01b (DONE): H12 mobility z per chain results/sv/SV01b/mobility_per_chain.csv; lattice contacts
  results/sv/SV01b/lattice_contacts_per_chain.csv / _per_residue.csv; shared module scripts/sv/common/lattice_contacts.py;
  model code scripts/sv/SV01b_models.py, model_data.csv (Q, X with lattice terms, B blocks).
- Core v3 = esr_analysis/resources/config/alignment_core_v3.json core_ranges (338-362, 371-396, 420-439, 441-456, 496-531).
- SV03 (DONE): core-v3 frame frozen in data/sv/SV03/core_v3_transforms.csv with scripts/sv/common/core_frame.py
  (load_transforms, to_core, load_ca_core); axis fit scripts/sv/common/helix_geometry.py (fit_axis); H12 map
  scripts/sv/common/h12_map.py (plot_h12_map, add_legends; marker column + hollow markers are parameters);
  per-chain geometry results/sv/SV03/geometry_per_chain.csv (H12 window 538-547; 'H11' = SV02 H10/11 497-526).
- PRE addition (DONE): ChEMBL API reachable; PyMOL NOT available (SV13 Part C uses its matplotlib fallback);
  data/sv/action_class_teammate.csv does not exist. Web tools (WebSearch/WebFetch) work for literature steps.
- SV11 (DONE): ligand-pocket site residues results/sv/SV11/site_residues.csv; steroid flag ligand_steroid.csv.
- SV10 (DONE): ligand groups (>= 80 % rule) results/sv/SV10/ligand_druglikeness.csv / entry_druglikeness.csv.
- SV12 (DONE): provisional classes results/sv/SV12/action_class.csv (entry, ccd, class, source, conflict, doi) and
  chain_class.csv (788 chains, class, cluster, primary_set); coverage results/sv/SV12/class_coverage.csv: only agonist
  (27 entries) and SERM (20) pass the >= 5-entry gate; SERD 4, antagonist 0, partial_agonist 0, other 0; 311 rows
  unlabelled. connectivity_only rows (primary_set False) stay out of primary analyses. ChEMBL activities on CHEMBL206:
  results/sv/SV12/chembl_activities.csv; E2 parent CHEMBL135. Shared HTTP client: scripts/sv/common/cached_http.py.
- SV14 (DONE): per-chain pocket/ligand volume and filling results/sv/SV14/pocket_per_chain.csv (primary set excludes
  ligand_incomplete); class styles scripts/sv/common/class_style.py (class colours never equal cluster colours).
- SV15 (DONE): per-ligand dpX to E2 results/sv/SV15/affinity_per_ligand.csv (relative for 14 ligands; agonists with
  dpX: DES +0.63, J3Z -0.78; 4OH absolute only; EST is the reference).
