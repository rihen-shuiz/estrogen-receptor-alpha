# Brief for step VIS1 (RUN_SV1028)

Extracted verbatim from docs/sv/RUN_SV1028.md: §1, §2 and the step text.

## 1. Rules for every step

1. **Scope.** Do only what the step says: no extra analyses, figures or refactors. Ideas go to `docs/EXPLORATION.md` under "Suggested, not done".
2. **Write locations.**
   - Read-only: existing files in `results/`, `figures/`, `data/clean/` and `esr_analysis/`.
   - New outputs go only to `results/sv/<ID>/`, `figures/sv/<ID>/` (PNG + SVG), `figures/panels/`, `visuals/`, `data/sv/`, `docs/sv/` and `scripts/sv/`.
   - Exception: SV00 may add *missing* files to `data/raw/` and `data/clean/`, never overwrite existing ones.
   - Steps that edit a doc say so explicitly.
3. **Decision rules.** Apply [DEFAULT] rules exactly as written. Never choose a new rule after seeing a result.
4. **Statistics.**
   - Run every test per chain and per entry.
   - Per-chain models use entry-clustered errors (GEE, or a bootstrap over entries).
   - When one step runs several tests, apply BH correction.
5. **Environment.**
   - Record package versions in `results/sv/<ID>/env.txt`.
   - Install only missing packages. Never upgrade or downgrade a package the existing pipeline uses, because that can break reproduction.
   - If a tool cannot be installed, use the fallback the step names and record it. Never emulate a tool silently.
6. **Tracker.**
   - Register each step as a card in `tracker/esr1.tasks`, in the existing card format. IDs: SV00, SV01a, SV01b, SV01c, SV02–SV15, SV08b, VIS1, VIS2. Check that each ID is unused first.
   - Write log entries only through `tracker/tl.py`.
7. **SUMMARY.md format.**
   - Line 1: anything not done, stubbed or substituted (tool, version, threshold, missing input), or "Nothing not done or substituted".
   - Then key numbers, the verdict and the file list.
8. **Figures.**
   - Colours: cl7 blue, cl4 orange, minor clusters green, as in the existing figures.
   - Each figure is drawn by a script from a saved CSV.
   - Titles and captions state what the data show. Never write an expected outcome into a figure.
   - Ligand markers use the `ligand_chemotype` column. In SV13–SV15 they use `action_class` (SV12), styled by `scripts/sv/common/class_style.py`, with unlabelled ligands in grey. Class colours never reuse the cluster colours.
   - Every SV figure script also saves each graph as its own PNG in `figures/panels/sv/<ID>/`. A graph is one axes together with its labels, legend and colorbar; crop it with the axes' tight bounding box.
9. **Ligand action classes.**
   - SV00–SV11 (including SV08) stay free of action classes, as written: describe ligands there only by chemistry (basic amine, carboxylic acid, heavy atoms).
   - SV08b and SV12–SV15 use action classes, taken only from `results/sv/SV12/action_class.csv`. SV09, VIS2 and FINAL may report those steps' class results.
   - Never assign a class from structure, cluster, chemistry or a figure. The class must stay independent of the outcomes it is tested against.
   - Wording: crystal structures show H12 order and position, not receptor degradation or cellular activity. Verdicts and captions say "H12 disorder" or "H12 mobility", never "kills the receptor".
10. **Inputs from earlier work.**
    - Find inputs through the output paths in the cards of `tracker/esr1.tasks`: Y07 heatmap, Y10 displacement profile, Y14, Y16 outliers, the shared MDS coordinates, S08 clusters, core v3, the helix table, the ligand table. Do not guess paths.
    - If an input is missing, the step is BLOCKED.

## 2. Order and dependencies

| # | Step | Needs |
|---|---|---|
| 1 | PRE — preflight | — |
| 2 | SV00 — data restore, shared chain table | PRE (network) |
| 3 | SV01a — completeness, first pass | uses the SV00 table if it exists, otherwise existing metadata |
| 4 | SV06 — outliers on the map | existing results only |
| 5 | SV10 — ligand drug-likeness | existing results only |
| 6 | SV11 — pocket sites, steroid flag | SV00 coordinates for Part A (see SV11) |
| 7 | SV02 — consensus secondary structure | SV00 coordinates |
| 8 | SV01b — completeness, mobility, lattice contacts | SV00 coordinates |
| 9 | SV03 — helix geometry | SV02 |
| 10 | SV04 — helix midpoints (Ye) | SV02 (reuses SV03 code if DONE) |
| 11 | SV05 — core v4 | SV02 |
| 12 | SV07 — spike audit | SV00, SV02 (contact columns from SV01b if DONE) |
| 13 | SV01c — completeness with the data-derived H12 window | SV01b, SV02 |
| 14 | SV08 — literature | web tools |
| 15 | PRE addition — preflight for the 2026-11-04 steps | — |
| 16 | SV12 — ligand IDs, action classes, ChEMBL activities | network (ChEMBL), unless the teammate file exists |
| 17 | SV14 — ligand size vs pocket size, filling | SV00 coordinates (SV11 site residues and SV12 classes if DONE) |
| 18 | SV15 — affinities across assay systems | SV12 (SV14 if DONE) |
| 19 | SV13 — H12 by action class | SV12, SV01b, SV02, SV03 (SV15 for C2b if DONE) |
| 20 | SV08b — literature for the class questions | web tools |
| 21 | SV09 — integration | every DONE step |
| 22 | VIS1 — cut multi-graph images | after all figure steps |
| 23 | VIS2 — collect valuable visuals | VIS1, SV09 |
| 24 | FINAL — checks and report | all |

- When a step is BLOCKED, mark its dependants SKIPPED ("needs <ID>") and run everything else. A dependency in brackets is optional: the step runs without it and says so in line 1 of its SUMMARY.
- Pending human decisions do not block anything:
  - The **H12 window** stays 538–547 for all primary analyses. SV02 proposes a data-derived window and SV01c tests it as a sensitivity analysis. The human signs it off later.
  - The **superposition core** stays v3 for SV03, SV04 and SV07. SV05 only compares v4 with v3.

### VIS1 — Cut multi-graph images into single graphs

Many PNG files contain several graphs. Split them so that each graph is its own file. Do not change the originals.

1. **Inventory.**
   - Every PNG in the project, except those under `data/`, `figures/panels/`, `visuals/`, virtual environments and hidden folders.
   - Images identical by SHA-1 (e.g. `results/task_outputs/` vs `task_outputs/`) are cut once; the manifest lists every path.
   - Skip images smaller than 200 px on either side.
2. **Definition of a graph:** one plot with everything needed to read it: its axes, tick labels, axis titles, own legend, own colorbar, and attached marginal axes (dendrogram, side bars, twin axes). Never separate these from their plot.
3. **Classify each image.**
   - `single` (one graph, incl. PyMOL renders) — not cut;
   - `multi` — cut.
4. **Method, in order of preference:**
   - **(a) SV figures:** already saved per graph by their scripts (rule 1.8). Only list them.
   - **(b) An SVG with the same stem exists** (matplotlib output):
     - Each axes is a top-level `<g id="axes_N">` group. Decide which groups form one graph by looking at the PNG.
     - For each graph, *hide* (do not delete) the other groups, so that shared glyph and clip-path definitions keep working.
     - Render at the PNG's pixel size with cairosvg or `rsvg-convert`, then crop the white margin, leaving 10 px padding.
   - **(c) PNG only:** split along full-width and full-height white gutters, recursively (XY-cut). Strips that hold only text (suptitle, panel letters) are not graphs.
5. **Shared figure-level legends or colour keys.**
   - If it sits next to a graph, include it in that graph's crop.
   - Otherwise save it once as `<stem>__legend.png` and reference it in the manifest for each graph that needs it.
6. **Only cropping.** No pixel edits and no redrawing of data. Rendering the original SVG is allowed.
7. **Naming:** `figures/panels/<parent>__<stem>__<nn>_<slug>.png`, where slug = the graph's title or axis label (≤ 40 characters, lowercase, `_`).
8. **Manifest:** `figures/panels/manifest.csv` with columns source_path(s), sha1, panel_path, method (sv_script / svg / png_gutter / single / not_cut), panel_title, needs_legend, notes.
9. **Check.**
   - For each multi image, build a contact sheet in `figures/panels/_check/`: the original on the left, the numbered panels on the right. Look at every sheet.
   - Each panel must contain the whole plot (axes, ticks, labels, legend, colorbar), with nothing cut off and no fragments of neighbours.
   - Fix failures. If a figure still fails, mark it `not_cut` with the reason; its original is then used as is.
   - Report the counts: single / cut into N panels / not_cut.
