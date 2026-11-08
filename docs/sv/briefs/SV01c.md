# Brief for step SV01c (RUN_SV1028)

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

### SV01 — H12 completeness: artifact or biology

Ask: *"check distribution of completeness, and if completeness correlates with any experimental conditions/parameters, or is it biology"*.

- The cluster is **not** a predictor: cluster labels come from the same coordinates whose completeness is the outcome. Analyse within cl7 and within cl4 instead.
- Outputs go to `results/sv/SV01a/`, `results/sv/SV01b/` and `results/sv/SV01c/`.

**SV01a — first pass**

Input: the SV00 chain table if it exists. Otherwise use the existing per-chain metadata from the ext run; the missing columns then move to SV01b.

1. **Distribution.**
   - per-residue observation over 305–550, by cluster and by construct mutation;
   - complete / partial / not modelled, per cluster, per chain and per entry;
   - for partial H12, which end is missing (N-side or C-side).
2. **Within-entry concordance:** in multi-chain entries, how often chains differ in H12 label and in cluster.
3. **Model.**
   - Outcome: H12 complete. [DEFAULT] binary, on the existing window 538–547. Run a sensitivity analysis with the fraction observed. Not multinomial: only 23 chains have no H12 modelled, too few per cluster.
   - Predictor blocks:
     - **Q** (quality): resolution, R-free, mean B
     - **X** (experiment): pH, growth temperature, ambient temperature, space group (P21 vs other), ESR1 chains in the ASU, dimer source
     - **B** (biology): peptide, partner has peptide, Y537S, L536S, D538G, other mutations, basic amine, carboxylic acid, ligand heavy atoms
   - Fit within cl7 and within cl4: at chain level with entry-clustered errors, and at entry level.
   - Report the McFadden gain of each block entered first and entered last, and odds ratios with 95 % CI.
   - If a predictor is constant within a cluster, or perfectly separates the outcome there, drop it in that cluster and list it (it is not a block failure).
4. **[DEFAULT] verdict.**
   - "artifact" if the Q+X gain exceeds the B gain in both entry orders and in both clusters;
   - "biology" if the reverse holds;
   - otherwise "mixed".

**SV01b — with coordinates**

5. **Local mobility as a second outcome.**
   - z-score the B of residues 528–537 and of the observed H12 residues against the same chain's core v3 B: `(B_res − mean_core) / sd_core`.
   - Model it with the same blocks (linear model, entry-clustered errors).
   - It is a response, not a predictor: as a predictor it would absorb the signal being explained.
   - Flag TLS-refined entries and run a sensitivity analysis without them.
6. **Lattice contacts** (gemmi): heavy-atom contacts (≤ 4 Å) of residues 528–550 with other chains, split into ASU partner vs symmetry mate. Save this as `scripts/sv/common/lattice_contacts.py`; SV07 reuses it.
7. Add the lattice terms to block X, refit both outcomes and re-issue the verdict.
8. **Only if the verdict is "mixed":** for each partial or not-modelled H12:
   - superpose the H12 of the cl7 medoid 5TM1|A|1 and of the cl4 medoid 7UJ7|C|1 onto the chain's core;
   - count clashes with lattice neighbours, to see whether the crystal blocks one of the two positions.

**SV01c — sensitivity to the H12 window**

9. If the SV02 data-derived H12 window differs from 538–547, re-run SV01b (both outcomes) on it. Report the verdict side by side with the primary verdict. If the windows are equal, mark SV01c SKIPPED ("same window").

**Outputs:** `completeness_profile.csv`, `concordance.csv`, `model_blocks.csv`, `odds_ratios.csv`, `verdict.md` (≤ 10 lines).

**Figures:**
- **Forest plot**, one per outcome (Y14 style). SV01a draws the completeness plot; SV01b redraws it and adds the mobility plot.
  - odds ratios (or coefficients) with 95 % CI, predictors grouped and coloured by block Q / X / B;
  - panels for cl7 and cl4; chain level as coloured dots, entry level as grey dots;
  - beside it, a bar chart of each block's gain, entered first and entered last.
- **Completeness profile:** observed fraction per residue 520–550, cl7 vs cl4.

**This subagent runs only part SV01c of SV01** (outputs in results/sv/SV01c/).
