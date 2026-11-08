# RUN_SV1028 — ESR1 LBD: the supervisor's requests (meeting 2026-10-28), one unattended run

You are Claude Code in the ESR1 project root. Execute every step of this file in order. Do not ask for confirmation. This file replaces `docs/PLAN_SV1028.md`.

Updated 2026-11-04: steps PRE addition, SV12–SV15 and SV08b added (section "Addition 2026-11-04", before SV09); rules 1.6, 1.8, 1.9, §2 and §3 changed; the new supervisor notes are in Appendix C.

---

## 0. Operating mode

1. **No questions, no approvals.**
   - Where this file gives a [DEFAULT] rule, apply it.
   - Where a rule cannot be applied or gives a degenerate result, mark the step BLOCKED with the reason and what was tried. Skip only the steps that depend on it (§2) and continue with the rest.
   - Stop the whole run only before an action outside this project folder or one that needs a password.
2. **Status file `docs/sv/RUN_STATUS.md`.**
   - One row per step: step | status (TODO / RUNNING / DONE / BLOCKED / SKIPPED) | one-line result | main outputs.
   - Create it in PRE with every step TODO. Update it after every step.
3. **Resume.** If `docs/sv/RUN_STATUS.md` already exists, continue from the first step that is not DONE, BLOCKED or SKIPPED. Never redo a DONE step.
4. **One subagent per step** from SV00 onward, run one at a time.
   - Give the subagent §1, §3 (only if the step needs data access), the step text verbatim and the input paths you found.
   - The subagent writes `results/sv/<ID>/SUMMARY.md` in the format of rule 1.7 and returns it.
   - The main session only orchestrates: it checks that the step's listed outputs exist and are not empty, updates RUN_STATUS and moves on.
5. **Code location.** Scripts for a step go to `scripts/sv/<ID>_*.py`; shared code goes to `scripts/sv/common/`. Every number in a verdict or summary comes from a saved file.

---

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

---

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

---

## 3. Data access (RCSB)

- **Downloads.** Use the project's `scripts/fetch_data.py` first: it worked in the original run. Keep its file format, so that the Q03 check finds the files.
  - Direct download fallback: `https://files.rcsb.org/download/{ID}.cif.gz` (or `.cif`).
- **Metadata.** Prefer reading the mmCIF file with gemmi:
  - crystal growth: `_exptl_crystal_grow.pH`, `.temp`
  - data collection: `_diffrn.ambient_temp`, `_symmetry.space_group_name_H-M`
  - refinement: `_refine.ls_d_res_high`, `.ls_R_factor_R_free`, `.ls_R_factor_R_work`
  - TLS refined = the `_pdbx_refine_tls` category exists
  - revisions: `_pdbx_audit_revision_history`
  - UniProt mapping: `_struct_ref_seq` and, if present, `_pdbx_sifts_unp_segments`
  - B and occupancy from `_atom_site`; handle altlocs as the existing code does.
- **Cleaning metadata values.**
  - `?` and `.` become NaN.
  - When there are several rows, use crystal / diffraction id 1 and flag the entry `multi_valued`.
  - [DEFAULT] plausibility: growth pH in 3–11, growth temperature in 270–315 K, ambient temperature in 50–320 K. Values outside these ranges become NaN and are counted in the summary.
- **Data API** (GraphQL, POST `https://data.rcsb.org/graphql`): only for fields missing from the mmCIF and as a cross-check.
  - Use the query in Appendix B, in batches of ≤ 50 entries.
  - Cache the raw JSON in `data/sv/api_cache/`.
  - Test the query once on 5TM1 first. If a field errors, introspect the schema; do not guess field names.
  - If `docs/RCSB_API/` exists, its query templates list more fields. The example responses in it are for 1DB1 (VDR), not ESR1 data.
- **Ligand SMILES** (if not already in the ligand table): `https://data.rcsb.org/rest/v1/core/chemcomp/{CCD}` → `rcsb_chem_comp_descriptor`.
- **SIFTS fallback:** `https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{id}`.
- **Primary-citation DOI:** mmCIF `_citation.pdbx_database_id_DOI` where `_citation.id` = `primary`; Data API `rcsb_primary_citation` only if the mmCIF has none.
- **ChEMBL API** (`https://www.ebi.ac.uk/chembl/api/data/`, JSON): resources `molecule` (filter `molecule_structures__standard_inchi_key`), `molecule_form`, `mechanism`, `target`, `document` (filter `doi`), `activity`, `assay`.
  - Test once on EST first, with the InChIKey from its RCSB chemcomp record. If a filter or field errors, read the returned JSON or the resource schema; do not guess field names.
  - Before using CHEMBL206, confirm that it is the single-protein target for UniProt P03372.
  - Page through all results (`page_meta.next`) and cache every response in `data/sv/chembl_cache/`.
- **Be polite to the servers:** ≤ 5 requests per second, retry with backoff, and cache every response.

---

## 4. Steps

### PRE — Preflight

1. Read `docs/RESULTS_SUMMARY.md`, `docs/METHODOLOGY.md`, `docs/limitations.md` and the card format in `tracker/esr1.tasks`.
2. Environment.
   - Find the Python environment the existing pipeline uses.
   - Check that these are available: numpy, pandas, scipy, scikit-learn, statsmodels, matplotlib, gemmi, rdkit, Pillow, and cairosvg or `rsvg-convert`.
   - Check for `mkdssp` ≥ 4 (conda-forge package `dssp` if conda exists).
   - Install what is missing (rule 1.5) and record everything in `docs/sv/env_pre.txt`.
3. Network: send one test request each to files.rcsb.org and data.rcsb.org.
   - If RCSB cannot be reached, mark SV00 BLOCKED and its dependants SKIPPED.
   - If web search is unavailable, SV08 will be SKIPPED.
4. Snapshot: write the SHA-1 of every file under `results/`, `figures/`, `data/clean/` and `esr_analysis/` to `docs/sv/readonly_snapshot.csv`. FINAL checks against it.
5. Create `docs/sv/RUN_STATUS.md`.
6. Supervisor feedback:
   - Write Appendix A verbatim to `docs/supervisor_feedback.md`. If the file exists, append it under a heading dated 2026-10-28.
   - Log W17 as done via `tl.py`.

---

### SV00 — Data restore and shared chain table

Prerequisite: `data/raw` is empty, which is why Q03 and S17 fail.

1. **Entry list.** Use the **frozen entry list of the original run**. Look for it in this order:
   - the manifest or list written by the original fetch (`data/`, `results/`, `esr_analysis/resources/`);
   - the tracker log entry of the fetch card;
   - the entries in `results/entry_clusters.csv`.

   Never define the set with a fresh RCSB search. A fresh search may be used only to report how many newer entries exist; they are not used.
2. **Restore and check the files.**
   - **Required:** every entry that holds one of the 788 kept chains (351 entries).
   - **Optional:** the rest of the frozen list, so that Q03 passes.
   - Write `data/sv/manifest.csv` with entry, checksum, download date and latest revision date.
3. **Counts.** Check against the original run and report any mismatch; a mismatch in the optional part does not block.
   - ESR1: 487 entries; 901 chains in 444 entries; 788 kept chains in 351 entries; 113 excluded, each with its reason.
   - ESR2: 39 entries, 86 chains.
4. **Regenerate missing data.**
   - Regenerate missing `data/clean/` files with the existing code.
   - Extract the full backbone (N, CA, C, O) of each kept chain to `data/sv/backbone/`.
5. **Reproduction gate.** Recompute the RMSD matrix with the existing code path, cluster at 1.5 Å, and compare with `results/rmsd_entry.npy` and `results/entry_clusters.csv`.
   - Report the ARI and the maximum |ΔRMSD|.
   - **Pass:** ARI = 1.0.
   - **[DEFAULT] tolerated:** at most 8 chains (1 %) differ, **and** every differing chain belongs to an entry revised after the original download date.
     - Flag those chains `remediated` and exclude them from all SV analyses.
     - List them in the SUMMARY.
   - **Otherwise:** SV00 is BLOCKED for coordinates. Still write the chain table (step 6), because SV01a can use its metadata.
6. **Chain table** `results/sv/SV00/chain_table.csv`: one row per kept chain. Reuse existing columns and record the source of every column in `column_sources.md`.
   - **Required columns:** chain id, entry, S08 cluster, QC tier, H12 label, peptide, construct mutations, ligand code(s), ligand heavy atoms, basic-amine flag (existing definition), resolution, R-free, R-work, mean B, pH, growth temperature, ambient temperature, space group, ESR1 chains in the ASU, release year, TLS-refined flag, UniProt–author numbering offset (SIFTS), `multi_valued`, `remediated`.
   - **New columns:**
     - carboxylic-acid flag, from SMARTS `[CX3](=O)[OX1H0-,OX2H1]` (acid or carboxylate only, no bioisosteres);
     - `ligand_chemotype` = `amine` / `acid` / `both` / `other` / `none` (no ligand), built from the two flags;
     - dimer partner = the chain contacting H11 (residues 497–530, heavy atoms ≤ 4 Å), with source `asu` / `symmetry` / `none`;
     - the partner's peptide and cluster.
7. **Observation mask.** Save `results/sv/SV00/observed_mask.npy` (chains × residues 305–550, Cα observed) with its row and column index files.

**Done when** the gate passes or is tolerated, the chain table has 788 rows, and missing values are reported per column.

---

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

---

### SV02 — Consensus secondary structure

Ask: *"consensus alignment of secondary structures"*.

1. Run DSSP on the full backbone of every kept chain with `mkdssp` ≥ 4.
   - Fallback: `pydssp`, noted in line 1. pydssp gives only H / E / loop, so the 3₁₀ and π fractions are then reported as not available.
2. Map to UniProt numbering with the SV00 offsets.
3. For each position 305–550, among the chains that observe it: fraction α (H), 3₁₀ (G), π (I) and other. Compute this globally and separately for cl7 and cl4.
4. **[DEFAULT] consensus helix:**
   - at least 4 consecutive positions with H in ≥ 80 % of observing chains;
   - positions observed in < 50 % of chains are marked "low support" and never define a boundary;
   - helices are named by overlap with the existing helix table (H1, H3–H12; ERα has no H2).
5. **Per-chain table** (chain × helix): start, end, observed helical residues, and whether the start or end is censored (adjacent to missing residues).
6. **Data-derived H12 window** = the global consensus H12 boundaries. Compare it with 538–547, 537–546 and 536–545. It is proposed only; a human signs it off.

**Outputs:** `sse_per_position.csv`, `consensus_{global,cl7,cl4}.csv`, `helix_boundaries.csv`, `sse_per_chain.csv`, `h12_window.md`.

**Figures:**
- **Main — helix track:**
  - x = residue 305–550 (UniProt), y = fraction α among observing chains;
  - lines for cl7, cl4 and all chains; the 80 % rule as a dashed line; the three historical H12 windows marked;
  - lower panel: fraction of chains observing each position, with low-support positions shaded.
- **Supplementary:** position × chain heatmap, sorted by cluster.

**Addition (2026-11-02) — Ramachandran plots per residue**

Ask: show the backbone shape of single key residues across all chains. One plot per residue; each point is that residue in one chain.

7. **Angles.**
   - Compute φ and ψ for every kept chain and every position 305–550 directly from the SV00 backbone (N, CA, C), with gemmi or a NumPy dihedral. Use UniProt numbering via the SV00 offsets.
   - A residue gets φ only if the previous residue is bonded to it (C–N ≤ 2.0 Å), and ψ only if the next residue is bonded to it.
   - DSSP is not needed for the angles. If mkdssp ran, compare its φ/ψ with these on 5 chains and record the largest difference.
   - Save `results/sv/SV02/phipsi_per_chain.csv`: chain, entry, cluster, position, phi, psi.
8. **Positions, fixed before plotting [DEFAULT]:**
   - main: 537 (H11–H12 hinge);
   - backup: 533, 534, 535, 536, 538, 539, 540 (hinge); 329, 359, 457 (`top_positions` in `results/ext/numbers_Y12.json`); 333 and 467 (the spike positions of SV07).
   - Do not add or drop positions after seeing the plots.
9. **One plot per position.**
   - x = φ, y = ψ, both −180° to 180°, equal aspect. One point per chain that has both angles at that position.
   - Colours as in rule 1.8 (cl7 blue, cl4 orange, minor clusters green); small, semi-transparent points; legend with n per cluster.
   - No allowed-region contours: the project has no reference distribution. Faint text labels "α" near (−60°, −45°) and "β" near (−120°, +130°) are allowed.
   - Title: the residue number and its helix or loop from the helix table.
   - Save each plot as PNG + SVG in `figures/sv/SV02/rama/`, and its panel copy in `figures/panels/sv/SV02/` (rule 1.8).
10. **Summary table** `results/sv/SV02/rama_summary.csv`, one row per position:
    - chains and entries per cluster that have both angles;
    - circular mean and circular SD of φ and ψ per cluster;
    - separation = √(Δφ² + Δψ²) between the cl7 and cl4 circular means, each Δ wrapped to ±180°;
    - overlap = share of cl7 and cl4 chains closer to the other cluster's mean than to their own, at chain level and at entry level (an entry's angles = circular mean over its chains).
11. **Check.** At 537 the circular-mean ψ per cluster must agree within 10° with the values recorded in `docs/EXPLORATION.md` (+157° in cl7, −47° in cl4). If they do not, state it in line 1 of the SUMMARY and report the new values.
12. **Resume.** If `docs/sv/RUN_STATUS.md` shows SV02 DONE but these outputs are missing, run 7–11 alone and add the status row "SV02 Ramachandran addition".

- SV09 treats the residue-537 plot as a main figure; the other 12 plots are backup.
- **Outputs (addition):** `phipsi_per_chain.csv`, `rama_summary.csv`, 13 plots.

---

### SV03 — Helix geometry

Ask: *"geometrical approach, more sensitive for things like the length of the helix (could affect ligand)"*.

Inputs: the SV02 tables; chains superposed on core v3.

1. **Helix length** per chain × consensus helix.
   - A terminus is censored when it is next to missing residues. A real end means the next 2 residues are observed and not helical.
   - Compare only uncensored lengths, and report the censoring fraction per helix.
2. **Helix axis** [DEFAULT]:
   - needs ≥ 6 observed helical Cα;
   - points = centroids of consecutive 4-residue windows; axis = principal line through them;
   - also fit H11 as two halves to get its kink angle.
3. **Descriptors.** Angles and distances do not depend on the frame; centroids are reported in the core frame.
   - angles H12–H3, H12–H11, H12–H4/5, H11–H3;
   - distances from the H12 centroid to the H3 axis and to the H11 axis;
   - H11 kink angle and H11–H12 hinge angle.
4. **Tests.**
   - (a) cl7 vs cl4 for each descriptor. H12 separation is expected; the informative part is the other helices.
   - (b) Within each cluster: association with ligand heavy atoms, basic amine, carboxylic acid, peptide and mutations. Use rank tests at chain and entry level, with BH correction.
5. **[DEFAULT] "more sensitive than RMSD"** means a descriptor separates a subgroup inside cl7 or cl4 with |rank-biserial r| ≥ 0.3 and q < 0.05 at both levels. Otherwise report "no added sensitivity".
6. **Figure — H12 map:**
   - x = distance from the H12 centroid to the H11 axis, y = angle H12–H3;
   - colour = cluster, marker = `ligand_chemotype`;
   - one point per chain that has an H12 axis.
   - Put the plotting function in `scripts/sv/common/` so that SV04 can reuse it.

**Outputs:** `geometry_per_chain.csv`, `tests.csv`, the H12 map, descriptor distributions by cluster.

---

### SV04 — Helix midpoints for incomplete helices (Ye)

Ask: *"compare incomplete residues by taking the middle atoms based on start and end of the helix"*.

1. Work in the core frame, on consensus helices with ≥ 6 observed Cα [DEFAULT]. Fit the axis as in SV03.
2. **Estimate midpoints.**
   - Project the observed Cα onto the axis and regress the projection on residue number (slope = rise per residue).
   - Use the fit to estimate the axis points of the consensus start and end residues.
   - Midpoint = the axis point at the middle residue.
   - Each chain is represented by the midpoint and unit axis vector of each helix.
3. **Validate first.**
   - On chains with complete helices, delete 1–6 residues from either end and measure the midpoint error and the axis-angle error against the full helix.
   - [DEFAULT] use a helix for a chain only if its observed length gives a median midpoint error ≤ 1.0 Å.
4. **Cluster on midpoints.**
   - Chain distance = RMS of midpoint distances over the helices both chains have; report the axis angle separately.
   - Average linkage.
   - ARI against S08 at K = 2 (cl7 and cl4 chains only) and at the cut that gives 7 clusters.
5. For partial-H12 chains: where they land (cl7-like or cl4-like), and whether this agrees with their S08 cluster.
6. **Figure:** the SV03 H12 map redrawn with the partial-H12 chains added as hollow markers at their estimated positions.

**Outputs:** `validation_error.csv` with the error curve, `midpoints.csv`, `ari.csv`, `partial_h12_assignment.csv`, `h12_map_partial.png`.

---

### SV05 — Stable core v4 from complete helices

Ask: *"take helices without the problem of completeness, find a stable core"*.

1. **Candidates:** positions in the SV02 global consensus helices that are observed in ≥ 98 % of kept chains [DEFAULT], excluding 528–550.
2. Superpose all chains on the candidates iteratively, starting from v3, until the mean shift is < 0.01 Å. Compute the spread per position.
3. **[DEFAULT] core v4:**
   - core v4 = positions with spread ≤ 0.5 Å;
   - if fewer than 80 positions remain, BLOCKED;
   - also report the 123 lowest-spread positions (the size of v3) for a like-for-like comparison.
4. **Compare with v3:** overlap (Jaccard) and which helices differ.
5. **Recompute with v4:**
   - the 788 × 788 RMSD with v4 superposition, using the existing code path that produced `core_sensitivity.csv`;
   - clusters at 1.0 / 1.5 / 2.0 / 2.5 Å, with ARI against S08 at each cut;
   - cluster sizes;
   - a Y10-style displacement profile on v4.

**Outputs:** `core_v4.json`, `spread.csv`, `ari.csv`, `rmsd_v4.npy`, one figure.

---

### SV06 — Outliers on the clustering map

Ask: *"color outliers in clustering (but no pattern is also possible)"*.

Inputs: the existing outlier flags (Y16) and the shared MDS coordinates. No new outlier methods.

1. **Figure, two panels:**
   - **A:** the shared MDS map with the 22 consensus outliers labelled (PDB|chain) and coloured by reason (quality / conformation, from Y16).
   - **B:** for each outlier, RMSD to its own medoid vs RMSD to the other medoid, with the diagonal drawn. Points near the diagonal sit between the clusters.
   - The reason counts (6 / 16) are already in Y16; do not redraw them.
   - Also add the same outlier marks as a side bar on the Y07 heatmap.
2. Test the number of outliers per cluster against expectation (hypergeometric).
3. "No pattern" is a valid result.

**Outputs:** the figure, `outliers_on_map.csv`, `verdict.md` (≤ 3 lines).

---

### SV07 — Spike audit on Y10 (slide 5)

Ask: *"check spikes on figure in slide 5"*.

1. **[DEFAULT] spike** = a position outside 528–550 where the cluster-mean or medoid-pair displacement is > 1.0 Å.
   - Expected: the H1–H3 loop near 333, the H8–H9 loop / start of H9 near 467, and medoid-line peaks around 400–460.
   - List every position that meets the rule.
2. **For each spike:**
   - chains observing it, per cluster (SV00 mask);
   - the displacement recomputed on chains that observe it in both clusters;
   - the top entries driving the mean, and their numbering offsets;
   - SSE context (SV02) and B z-score;
   - lattice and partner contacts (SV01b module; "not available" if SV01b is not DONE).
3. **[DEFAULT] verdict per spike**, one of four:
   - sparse-sample effect (gone when restricted);
   - driven by ≤ 3 entries;
   - numbering error;
   - genuine loop difference (persists across many entries). For these, test against peptide, ligand chemotype and dimer source.

**Outputs:** `spikes.csv`, `verdict.csv`.

**Figure — spike profile:**
- the Y10 profile with two curves on the same axis: all chains (as in Y10) and only chains observing the position in both clusters;
- the observation fraction in a panel below;
- each spike labelled with its verdict.

This is a diagnostic only; the main analysis keeps all chains.

---

### SV08 — Literature

Ask: *"check literature if only H12 is significant"*.

**Tools:** use WebSearch and WebFetch. If they are unavailable, write the questions to `docs/sv/literature_TODO.md` and mark the step SKIPPED.

**Questions:**
1. Besides H12, which ERα LBD regions are reported to differ between crystal structures (H3, H11, H11–H12 loop, H1–H3 loop, H8–H9 loop)? What are those differences associated with: coactivator peptide, construct mutations, ligand size?
2. How is missing H12 density reported and interpreted in ERα LBD crystal structures?
3. How do Y537S, L536S and D538G act on H12?
4. Is flexibility known near residues 333 and 467?

**Rules:**
- Open every source you cite and take each claim from the opened text.
- Skip any claim that is framed by ligand action class (rule 1.9).
- Stop a question once another source would not change the answer.

**Output:** `docs/sv/literature.md`.
- A table: question | claim | source (DOI or URL, opened) | which of our results it bears on.
- A list of sources that could not be opened.

Nothing enters SV09 unless its source was opened.

---

### SV10 — Ligand drug-likeness table (addition 2026-11-02)

Ask (lecture): describe the crystallised ligands with standard drug-likeness parameters (Lipinski, HBA, ALERTS and similar) and compare cl7 with cl4.

**Running this addition**
- It needs only existing results. Run it right after SV06.
- Add the row SV10 to `docs/sv/RUN_STATUS.md` (in PRE, or when this addition is first read if the file already exists), and register the card SV10 in `tracker/esr1.tasks` (rule 1.6; check the ID is unused).
- Tool: RDKit (already in the PRE list). Record its version in `env.txt`.

**Inputs** (confirm the paths through the tracker cards, rule 1.10)
- SMILES: `results/ligands.csv` (`canonical_isomeric_smiles`, else `smiles`).
- Group and size: `results/ext/ligand_descriptors.csv` (`structural_majority`, `heavy_atoms`). Groups are those of the existing ligand analysis: a ligand belongs to cl7 or cl4 when ≥ 80 % of its chains are in that group. Ligands that are mixed, minor or without a matched chain get descriptors but are not compared.
- Potency: `results/potency_by_cluster.csv` (`pchembl`, one value per ligand).

1. **Descriptors per ligand** (RDKit, definitions fixed):
   - Lipinski: MW (`Descriptors.MolWt`), cLogP (`Crippen.MolLogP`), HBD (`Lipinski.NumHDonors`), HBA (`Lipinski.NumHAcceptors`); violations = count of MW > 500, cLogP > 5, HBD > 5, HBA > 10.
   - Veber: rotatable bonds (`rdMolDescriptors.CalcNumRotatableBonds`), TPSA (`rdMolDescriptors.CalcTPSA`); violations = count of rotatable bonds > 10, TPSA > 140.
   - QED: `QED.qed` and its eight properties from `QED.properties` (MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, ALERTS).
   - Structural alerts: QED ALERTS; Brenk and PAINS (A, B, C) match counts with `FilterCatalog`, and the names of the matched alerts.
   - Efficiency, only for ligands with pChEMBL: LE = 1.37 × pChEMBL / heavy atoms; LipE = pChEMBL − cLogP.
   - A SMILES that does not parse gets an empty row with the reason; count these in the SUMMARY.
2. **Comparison cl7 vs cl4**, at two levels:
   - ligand level (primary): one row per ligand;
   - entry level: one row per entry, with the descriptors of its ligand of interest and the entry's majority group.
   - Chain level is not used: a ligand's descriptors are the same in all its chains.
   - Continuous parameters: median [IQR] per group, Mann–Whitney. Counts of violations or alerts: share of ligands with ≥ 1, Fisher's exact test. BH over all parameters within each level.
3. **[DEFAULT] size check.** Most of these parameters grow with ligand size, which already differs between the groups.
   - For each parameter, report its Spearman ρ with heavy atoms.
   - Fit rank(parameter) ~ group + rank(heavy atoms) by OLS and report the group p, BH-corrected.
   - Label a parameter "adds beyond size" if this size-adjusted q < 0.05 at both levels; otherwise "restates size".
4. **Wording.** Describe ligands by chemistry only (rule 1.9). These parameters describe oral drug-likeness, not binding, conformation or activity; the set includes the natural hormone, fragments and tool compounds. Potency mixes assay types; say so next to LE and LipE.

**Outputs** (`results/sv/SV10/`)
- `ligand_druglikeness.csv`: one row per ligand, all descriptors, alert names and group.
- `druglikeness_by_group.csv` and the same table as `druglikeness_by_group.md`: one row per parameter with its definition; cl7 n and median [IQR] or share; cl4 the same; p and q at ligand level; p and q at entry level; Spearman ρ with heavy atoms; size-adjusted q; label.
- `verdict.md` (≤ 5 lines): which parameters add beyond size.
- No figure.

SV09 mentions SV10 in one sentence and links the table; it is not a main item.

---

### SV11 — Pocket sites and steroid flag (addition 2026-11-02)

Ask (lecture notes): label each predicted pocket by the site it occupies, and add a ligand chemical-class label.

**Running this addition**
- Run it right after SV10.
- Part A needs the SV00 structures for the groove and dimer-interface residue sets; everything else needs only existing results. If SV00 is BLOCKED for coordinates, run the rest, label pockets only as "ligand pocket" or "other surface", and say so in line 1 of the SUMMARY.
- Add the row SV11 to `docs/sv/RUN_STATUS.md` and register the card SV11 in `tracker/esr1.tasks` (rule 1.6; check the ID is unused).

**Part A — which site does each P2Rank pocket occupy?**

Inputs (confirm the paths through the tracker cards, rule 1.10):
- every predicted pocket per chain, all ranks: `results/p2rank/<entry>_<chain>.pdb_predictions.csv` (`residue_ids` in author numbering);
- ligand contacts per chain: `results/chain_ligand_of_interest.csv` (`pocket_residues`); use them as they are and record the distance cutoff from the script that wrote them;
- the SV00 chain table (peptide, dimer partner), the SV00 structures, and the SV00 numbering offsets (author → UniProt).

1. **Check the P2Rank input.** Confirm from `scripts/p2rank_batch.py` and `scripts/clean_chain_pdbs.py` that P2Rank saw only the receptor chain and its ligand of interest, with no peptide and no partner chain. Write what you found in the SUMMARY.
2. **Site residue sets, fixed before any pocket is labelled [DEFAULT].** A residue belongs to a site if it is in contact (any heavy atom ≤ 4 Å) in ≥ 50 % of the relevant chains:
   - ligand pocket: contact with the ligand of interest, among chains that have one;
   - AF-2 coactivator groove: contact with a coactivator peptide chain, among chains that have one (peptide = a non-receptor polymer chain of ≤ 30 residues), computed on the SV00 structures;
   - dimer interface: contact with the SV00 dimer partner, among chains that have one.
   - A residue may belong to more than one site. Save each set with the support (share of chains) of every residue.
3. **Label every pocket** (all ranks, UniProt numbering): the site that holds the largest share of the pocket's residues, if that share is ≥ 0.5 [DEFAULT]; otherwise "other surface". On a tie, record both sites.
4. **Per chain and site:** whether any pocket has that label, its best rank, best probability and SAS points.
5. **cl7 vs cl4, per site:** the share of chains with a pocket there, and the best rank and probability. Test per chain (entry-clustered errors) and per entry, with BH correction (rule 1.4).

Outputs (`results/sv/SV11/`): `site_residues.csv`, `pocket_labels.csv` (chain, rank, label, share, probability, SAS points), `site_by_chain.csv`, `site_tests.csv`, `verdict.md` (≤ 5 lines).

Figure: share of chains with a pocket at each site, cl7 vs cl4 (bars, colours as in rule 1.8). The caption states the shares only.

**Part B — steroid flag**

6. Flag ligands that contain a steroid ring system, using the SMILES in `results/ligands.csv` as in SV10 and this SMARTS (carbon atoms, any bond order):
   `[#6]1~[#6]~[#6]~[#6]2~[#6](~[#6]~1)~[#6]~[#6]~[#6]1~[#6]~2~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]~1~2`
   - Check it before use: EST must match; DES, GEN, OHT and RAL must not. If the check fails, mark Part B BLOCKED.
7. Compare the share of steroid ligands in cl7 vs cl4 at ligand and entry level (Fisher's exact test), with the groups used in SV10.
8. Keep the flag as its own column; do not change `ligand_chemotype` from SV00. The label names the ring system only (rule 1.9).

Outputs (`results/sv/SV11/`): `ligand_steroid.csv`, `steroid_tests.csv`.

SV09 gives Part A one paragraph and Part B one sentence; the SV11 figure is backup.

---

### Addition 2026-11-04 — action classes, pocket filling, affinities (overview)

Asks (meeting notes, Appendix C):
- ligand size vs pocket size, with apo as the pocket 100 % and E2 as the ligand 100 %, then the filling percentage;
- comparing affinities that come from different assay systems;
- the 6 action classes, taken from the papers (DOIs), shown on the plots;
- whether antagonists destabilise H12, with pictures of H12 and the ligand per class, and whether structure links to function;
- whether H12 differs between agonist subtypes (the VDR superagonist point);
- SMILES / InChIKey and the ChEMBL API.

Not part of this run: affinity prediction from folding. No step waits for it.

**Running this addition**
- Steps: PRE addition, SV12, SV14, SV15, SV13, SV08b, in the §2 order.
- Add their rows to `docs/sv/RUN_STATUS.md` when this addition is first read, and register the cards (rule 1.6; check each ID is unused).
- If SV09, VIS1, VIS2 or FINAL are already DONE at that moment, add the rows "SV09 addition", "VIS1 addition", "VIS2 addition" and "FINAL addition" after SV08b. Run them only for the new steps, following "Integration of the addition" below and the original step texts. Never redo the DONE parts.

### PRE addition — preflight for the new steps

Run it before SV12, also when PRE is already DONE.

1. Check for PyMOL (the open-source build is fine); install it if missing (rule 1.5). If it cannot be installed, SV13 Part C uses the fallback named there.
2. Network: one test request to the ChEMBL API (§3). If ChEMBL cannot be reached and `data/sv/action_class_teammate.csv` does not exist, mark SV12 BLOCKED.
3. Append Appendix C verbatim to `docs/supervisor_feedback.md` under the heading "Meeting notes received 2026-11-04".
4. Record versions and the network results in `docs/sv/env_pre_1104.txt`.

---

### SV12 — Ligand identifiers, action classes and ChEMBL activities

Ask: *"get the SMILES, InChIKey from the API and use the ChEMBL API"*; *"6 different action class types, based on their DOIs"*.

**Unit:** one row per entry × ligand of interest (`results/chain_ligand_of_interest.csv`). Each chain takes the label of its own ligand of interest.

1. **Identifiers.**
   - Per CCD code: stereo SMILES and InChIKey from the RCSB chemcomp record (§3).
   - Per entry: the primary-citation DOI (§3).
2. **ChEMBL match.**
   - Full InChIKey → `molecule` → its parent (`molecule_hierarchy`).
   - No full match: match on the first InChIKey block (connectivity) and flag `connectivity_only`. Stereoisomers can differ in activity, so these rows are left out of the primary analyses of SV13 and SV15 and used only in a sensitivity run.
   - No match: `no_chembl`.
3. **Label set [DEFAULT]**, provisional, 6 classes: `agonist`, `partial_agonist`, `SERM`, `SERD`, `antagonist`, `other`. Rows without a label are `unlabelled`, which is not a class.
4. **Sources, in order of precedence.**
   - **(a) Teammate file** `data/sv/action_class_teammate.csv` (entry, ccd, action_class, doi), if it exists at run time. Use its labels and class names as they are. The ChEMBL labels then serve only as a cross-check (`class_agreement.csv`).
   - **(b) ChEMBL mechanism** (curated; mostly drugs and clinical compounds): records for the parent whose target contains UniProt P03372 (CHEMBL206, or a family or complex target containing it). Map [DEFAULT]:
     - AGONIST → agonist; PARTIAL AGONIST → partial_agonist;
     - mechanism text contains "degrad" or "downregulat", or the action type is DEGRADER → SERD;
     - MODULATOR, mechanism text containing "modulator", or both agonist-type and antagonist-type records → SERM;
     - ANTAGONIST, INVERSE AGONIST, NEGATIVE ALLOSTERIC MODULATOR → antagonist;
     - any other action type → other;
     - several records: SERD > SERM > antagonist > partial_agonist > agonist > other.
   - **(c) The entry's own paper (DOI).** ChEMBL `document` with the entry's DOI → activities of the ligand (parent or its forms) on ESR1 targets in that document with pChEMBL ≥ 6 [DEFAULT] → their assay descriptions, lower-cased, matched by these patterns:
     - degradation: `degrad|downregulat|down-regulat`;
     - antagonism: `antagonis|inhibition of (e2|estradiol|17.?beta-estradiol)[- ]induced`;
     - agonism: `(?<!ant)agonis|transactivation|transcriptional activation`, only in descriptions with no antagonism match ("agonist" is part of "antagonist", hence the lookbehind).
     - Class: degradation → SERD; agonism and antagonism → SERM; antagonism only → antagonist; agonism only → agonist; nothing → unlabelled. Route (c) cannot give partial_agonist or other.
     - **Validation gate [DEFAULT]:** compare (c) with (b) on rows that both label, counting partial_agonist as agonist. Use route (c) only if it agrees on ≥ 80 % of ≥ 10 such rows; otherwise its labels are not used, and line 1 says so.
   - **Conflicts:** the higher source wins. Flag `class_conflict` and list the conflicts in the SUMMARY.
5. **Activities for SV15.** For every matched parent, fetch all ChEMBL activities on CHEMBL206. Keep: activity, assay and document IDs, assay type, assay description, standard type / relation / value / units, pChEMBL, data-validity comment, potential-duplicate flag.
6. Map E2 (`EST`) to its ChEMBL parent; SV15 uses it as the reference.

**Coverage gate [DEFAULT]** (used by SV13–SV15): a class comparison is tested only if each side has ≥ 5 entries. Smaller classes are drawn but not tested.

**Outputs** (`results/sv/SV12/`): `ligand_ids.csv` (ccd, smiles, inchikey, ChEMBL parent, match type), `entry_doi.csv`, `action_class.csv` (entry, ccd, class, source, conflict flag, doi), `class_coverage.csv` (chains and entries per class, source and cluster), `route_c_validation.csv`, `chembl_activities.csv`, and `class_agreement.csv` if the teammate file exists.

---

### SV14 — Ligand size vs pocket size, and pocket filling

Ask: *"correlation between ligand size and pocket size; apo is the baseline, normalised to 100 %; the natural ligand binding ER (E2) is 100 % for ligand volume; then the percentage of filling"*.

Inputs: SV00 structures and chain table; the ligand of interest per chain (`results/chain_ligand_of_interest.csv`); the SV11 ligand-pocket residue set (`site_residues.csv`; if SV11 is not DONE, compute it with SV11 rule 2, first bullet); SV12 classes if DONE. Work in the core v3 frame (the SV03 superposition, or recomputed with its code). "Wild-type" below means no construct mutations in the SV00 chain table.

1. **Atoms and grid [DEFAULT].**
   - Receptor chain only: heavy atoms, highest-occupancy altloc. No ligand, waters, ions, peptide or other chains (as P2Rank saw it, SV11).
   - Bondi radii: C 1.70, N 1.55, O 1.52, S 1.80, F 1.47, Cl 1.75, Br 1.85, I 1.98, other 1.80 Å. No hydrogens on either side, so ligand and pocket volumes follow the same convention.
   - Grid spacing 0.4 Å, over the box of the ligand-pocket residues + 6 Å, wide enough for ligand side chains that leave the pocket.
2. **Pocket volume.**
   - Probe centres = grid points ≥ r_vdW + 1.4 Å from every atom. Cavity = grid points within 1.4 Å of a probe centre (solvent-excluded space).
   - Buried = ≥ 10 of 14 rays (±x, ±y, ±z and the 8 cube diagonals) hit a protein vdW sphere within 10 Å.
   - Pocket = the 6-connected component of buried cavity points that overlaps the ligand of interest. For chains without a ligand, the component holding the grid point nearest the E2 reference centroid (median EST centroid over wild-type EST chains).
   - Volume = number of points × 0.064 Å³.
   - Flag `pocket_open` when the component borders non-buried cavity points, and record how many.
   - Sensitivity: buriedness ≥ 8 and ≥ 12; report the Spearman ρ of each with the primary volume.
3. **Ligand volume.** vdW volume of the modelled heavy atoms of the ligand of interest, on the same grid. Flag `ligand_incomplete` when fewer heavy atoms are modelled than the CCD has; these chains stay out of the primary analysis.
4. **Filling.** filling % = 100 × (ligand volume inside the pocket component) / pocket volume. Also report the share of ligand volume outside the pocket (side chains leaving it).
5. **References [DEFAULT].**
   - Ligand 100 % = the median EST ligand volume over wild-type EST chains.
   - Apo entries = `ligand_chemotype` = none and no non-water, non-ion residue of ≥ 6 heavy atoms within 4 Å of the ligand-pocket residues. Report how many pass (about 7 entries expected).
   - Pocket 100 % = the apo baseline: median of entry medians over apo entries, with its IQR.
   - If fewer than 5 apo entries pass, the pocket 100 % is the median pocket of wild-type EST chains instead; line 1 says so.
   - The E2-bound pocket is reported as a second reference in every case.
6. **Analysis.**
   - Spearman ρ of ligand volume vs pocket volume at ligand level (median per ligand) and entry level: all chains, within cl7, within cl4, and per class (coverage gate).
   - Filling % by cluster and by class: Mann–Whitney or Kruskal–Wallis at entry level, bootstrap over entries at chain level, BH.

**Outputs** (`results/sv/SV14/`): `pocket_per_chain.csv` (chain, pocket volume, sensitivity volumes, `pocket_open`, ligand volume, inside share, filling %, % of each reference), `references.json`, `correlations.csv`, `tests.csv`, `verdict.md` (≤ 5 lines).

**Figures:**
- **Main:** x = ligand volume (% of E2), y = pocket volume (% of apo). Colour = cluster, marker = class (`ligand_chemotype` if SV12 is not DONE). Draw the apo-baseline IQR as a band around 100 % and dashed lines of equal filling at 25 / 50 / 75 %.
- Filling % by class and by cluster (box with points, entry level).

---

### SV15 — Affinities across assay systems

Ask: *"different affinities come from different assays (in vitro, …) and give different IC50 …; think how to compare across systems"*.

Inputs: SV12 `chembl_activities.csv` and ligand map; SV14 if DONE.

1. **Clean [DEFAULT].** Target CHEMBL206 only; standard relation "="; pChEMBL present; no data-validity comment; not a potential duplicate; no `connectivity_only` rows in the primary set.
2. **Systems.**
   - **Binding:** assay type B with Ki, Kd or IC50.
   - **Functional:** assay type F. Listed in the table only, never compared across classes: an agonist's EC50 and an antagonist's IC50 measure opposite effects.
   - Assay types A, T, P and U are dropped.
3. **Relative affinity to E2 in the same assay (primary).** For every binding assay that measured both the ligand and E2: ΔpX = pX(ligand) − pX(E2), and RBA = 100 × 10^ΔpX. In one competition assay the IC50 ratio follows the Ki ratio, so Ki, Kd and IC50 assays all qualify. Per ligand: the median ΔpX and the number of assays.
4. **Fallback.** Ligands with no same-assay E2 pair get the median binding pChEMBL (Ki/Kd before IC50), flagged `absolute`. They are drawn as hollow markers and never enter the tests.
5. **Check.** For ligands with ≥ 2 paired assays, compare the SD of ΔpX across assays with the SD of absolute pX across the same assays. This shows how much the E2 reference reduces the spread between assay systems.
6. **Comparisons** (E2 itself excluded, being the reference): ΔpX by cluster (the ligand groups of SV10), by class (SV12, coverage gate), and Spearman ρ with filling % and ligand volume (if SV14 is DONE). Ligand level is primary, entry level secondary; BH.
SV15 supersedes `results/potency_by_cluster.csv` only inside the SV section; the earlier file and its numbers stay untouched.

**Outputs** (`results/sv/SV15/`): `affinity_per_ligand.csv`, `assay_spread.csv`, `tests.csv`, `verdict.md` (≤ 5 lines).

**Figures:** ΔpX per ligand by cluster and by class (strip + box; hollow = `absolute`); ΔpX vs filling % if SV14 is DONE.

---

### SV13 — H12 by action class

Ask: *"can an antagonist kill the receptor by heavily destabilising H12? Show H12 and the ligand and compare several class types. Is the structural difference linked to function?"*; *"show the 6 classes in the plot"*; the VDR superagonist point (no H12 difference found between agonist types).

Inputs: SV12 `action_class.csv`; SV00 chain table; SV01b outcomes (H12 complete, fraction observed, H12 mobility z); SV02 φ/ψ at 537; SV03 geometry and uncensored H12 lengths; the shared MDS coordinates; the core v3 superposition (from SV03; recompute with its code if the transforms were not saved); SV15 ΔpX for C2b if DONE. H12 window: 538–547 (§2).

**Part A — classes on the existing maps**
1. `class_by_cluster.csv`: chains and entries per class × S08 cluster.
2. Redraw the shared MDS map and the SV03 H12 map with colour = cluster and marker = class (rule 1.8). Reuse the SV03 plotting function.

**Part B — tests, fixed in advance**
3. **C1 — destabilisation:** SERD vs SERM, and antagonist vs SERM.
   - Outcomes: H12 complete, H12 fraction observed, H12 mobility z.
   - Strata: all chains, and within cl4.
4. **C2 — agonist subtypes (the VDR superagonist point)**, within cl7:
   - **C2a:** agonist vs partial_agonist;
   - **C2b:** among agonists, Spearman ρ of each outcome with ΔpX (SV15; agonists that bind more strongly than E2 are the closest analogue of superagonists);
   - Outcomes: H12 mobility z, angle H12–H3, distance from the H12 centroid to the H11 axis, uncensored H12 length, ψ at 537 (as the wrapped difference from the cl7 circular mean).
5. **Models** (C1, C2a): outcome ~ class + Q + X (the SV01b blocks, with lattice terms) + construct mutations + peptide. Fit at chain level with entry-clustered errors, and at entry level. Peptide and mutations are included because they are spread unevenly across classes and act on H12 themselves. A covariate that is constant or perfectly separating in a stratum is dropped there and listed, as in SV01. C2b uses rank correlation at entry level and at chain level (bootstrap over entries).
6. Report each class effect (OR or coefficient with 95 % CI, or ρ) and its p; BH within each comparison.
7. **[DEFAULT] verdict per comparison:** "linked" if the class term (or ρ) has q < 0.05 at both chain and entry level; "not linked" otherwise; "not testable" if the coverage gate fails. Wording as in rule 1.9.

**Part C — pictures of H12 and the ligand**
8. **One view for all panels.** Orient on core v3 plus the H12 of the cl7 medoid 5TM1|A|1, save the view in `results/sv/SV13/view.txt`, and use it unchanged for every panel.
9. **One panel per class:**
   - every chain of the class superposed on core v3, shown as thin, semi-transparent Cα traces of residues 528–550;
   - the class representative [DEFAULT]: the chain with the lowest median RMSD (existing matrix) to the other chains of its class, ties broken by resolution. Draw it as a cartoon with H12 coloured by B z-score and its ligand of interest as sticks;
   - unmodelled residues are simply absent: no rebuilding and no interpolation;
   - the caption states the number of chains, entries and chains with complete H12.
10. One overlay panel with all class representatives in the same view.
11. Each panel as its own PNG (rule 1.8), plus one grid image.
12. Fallback without PyMOL: Cα traces projected with matplotlib onto the same saved view, stated in line 1.

**Outputs** (`results/sv/SV13/`): `class_by_cluster.csv`, `tests.csv`, `representatives.csv`, `view.txt`, `verdict.md` (≤ 6 lines).

**Figures:** a forest plot of the class effects (C1, C2a; chain level coloured, entry level grey) with the C2b ρ values beside it; the per-class H12 panels; the two redrawn maps.

---

### SV08b — Literature for the class questions

Same tools and rules as SV08, except that claims framed by action class are allowed (rule 1.9).

Questions:
1. In VDR crystal structures with superagonists vs agonists, does the H12 position differ? (The professor's point; verify it before it goes on a slide.)
2. Is H12 disorder or destabilisation reported for SERDs compared with SERMs in ERα (crystal, HDX, NMR or MD), and is it linked to receptor degradation?
3. Which volumes are reported for the ERα ligand pocket and for E2, and with which method? (Context for the SV14 references.)

**Output:** `docs/sv/literature_class.md`, in the SV08 format.

---

### Integration of the addition
- **SV09:** one paragraph per new step. The main figures for the next meeting are, in addition to those listed in SV09: the SV13 per-class H12 panels, the SV13 forest plot and the SV14 main figure. The SV13 maps and SV15 figures are backup. Literature from SV08b enters only if its source was opened.
- **VIS2:** always include those three main figures.
- **FINAL, "Decisions for the team":** replace the provisional ChEMBL classes with the teammate's labels when they are ready, then re-run SV13–SV15.

---

### SV09 — Integration

1. Write one paragraph and one figure per step, from the SUMMARY and verdict files.
   - Main figures for the next meeting: the SV01 forest plot, the SV02 helix track, the SV03 H12 map with the SV04 overlay, the SV07 spike profile and the SV06 outlier panels. All other figures are backup.
   - Add the text as a new section "Supervisor requests 2026-10-28" at the end of `docs/RESULTS_SUMMARY.md`. Leave the earlier text and numbers untouched.
   - Add new limitations to `docs/limitations.md`.
2. Check by script that every number in the new section matches the CSV it comes from. Write the check to `results/sv/SV09/number_check.csv`.

---

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

---

### VIS2 — Collect the valuable visuals into one folder

1. **Pool:** the panels from VIS1, the `single` images and the `not_cut` originals.
2. **[DEFAULT] valuable** = the graph supports a number or claim stated in `docs/RESULTS_SUMMARY.md` (including the new SV section) or in an SV verdict file.
   - Keep one graph per claim: the clearest one.
   - Prefer the keep-all (788-chain) version over the September (351-chain) version, and the declared 1.5 Å threshold over the alternative thresholds.
   - Always include the five main SV figures (SV09) that exist.
   - Exclude: QC or diagnostic graphs tied to no claim, duplicates, superseded versions, and legend-only files. A needed legend goes next to its graph with the same number and the suffix `b`.
3. **Copy, never move,** into `visuals/` (one flat folder, outside `results/`).
   - Name files `NN_<slug>.png`, numbered in the order of the story in `RESULTS_SUMMARY.md`.
4. **`visuals/INDEX.md`** — a table: NN | file | what it shows (one line, factual) | claim supported (the number or sentence, with its RESULTS_SUMMARY section or verdict file) | source path.
5. **Check by script:**
   - every INDEX claim text exists in the cited document;
   - every file in `visuals/` is listed in INDEX and opens;
   - no two files have the same SHA-1.

---

### FINAL — Checks and report

1. Run `tl.py validate`; it must pass.
2. **Read-only check:** recompute the SHA-1 of `docs/sv/readonly_snapshot.csv`. List every changed or deleted file in line 1 of the report. Files that SV00 added to fill gaps are fine.
3. **Write `docs/sv/RUN_REPORT.md`:**
   - **Line 1:** everything not done, stubbed or substituted across the run (BLOCKED or SKIPPED steps, fallbacks, excluded chains), or "Nothing not done or substituted".
   - **Table:** step | status | key number or verdict | main files.
   - **"Decisions for the team":**
     - sign off the H12 window (SV02 proposal vs 538–547, and the SV01c effect);
     - adopt core v4 or keep v3 (SV05);
     - review the `visuals/` selection.
4. **Final message:** print `RUN_STATUS.md`, the `tl.py validate` output and line 1 of `RUN_REPORT.md`.

---

## Appendix A — Supervisor feedback (verbatim, team notes after the meeting on 2026-10-28)

```
[10/28/26 1:08 PM] ayau🫧 (taylor's version): we can check distribution of completeness, and if completeness correlates with any experimental conditions/parameters, or is it biology?
we can try geometrical approach since it's more sensitive for things like the length of the helix (could affect ligand)
consensus alignment of secondary structures
take helices without the problem of completeness, find a stable core
color outliers in clustering (but no pattern is also possible)
check spikes on figure in slide 5
check literature if only H12 is significant? so far seems so
[10/28/26 1:13 PM] ayau🫧 (taylor's version): also, for friday we can use that ppt but we have to extract extra data and leave only the main findings etc
[10/28/26 1:15 PM] ㅤ`Ye: I also remember the point about somehow being able to compare incomplete residues, for example, taking the middle atoms (based on start and end of the helix) and applying methods on these data. Huge limitations exist but we are playing around
```

## Appendix B — RCSB Data API query (test on 5TM1 first)

POST `https://data.rcsb.org/graphql`, body `{"query": "<below>", "variables": {"ids": ["5TM1", ...]}}`

```graphql
query ($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    exptl { method }
    exptl_crystal_grow { crystal_id pH temp method }
    diffrn { id crystal_id ambient_temp }
    symmetry { space_group_name_H_M }
    cell { length_a length_b length_c angle_alpha angle_beta angle_gamma }
    refine { ls_d_res_high ls_R_factor_R_work ls_R_factor_R_free }
    rcsb_accession_info { deposit_date initial_release_date }
    pdbx_audit_revision_history { ordinal revision_date major_revision minor_revision }
    polymer_entities {
      rcsb_polymer_entity_container_identifiers { entity_id auth_asym_ids uniprot_ids }
      entity_poly { rcsb_mutation_count pdbx_seq_one_letter_code_can }
      rcsb_polymer_entity_align {
        reference_database_name reference_database_accession
        aligned_regions { entity_beg_seq_id ref_beg_seq_id length }
      }
      polymer_entity_instances {
        rcsb_polymer_entity_instance_container_identifiers {
          auth_asym_id auth_to_entity_poly_seq_mapping
        }
      }
    }
    nonpolymer_entities {
      nonpolymer_comp { chem_comp { id name formula_weight } }
      nonpolymer_entity_instances {
        rcsb_nonpolymer_entity_instance_container_identifiers { auth_asym_id auth_seq_id }
        rcsb_nonpolymer_instance_validation_score { RSCC is_subject_of_investigation }
      }
    }
  }
}
```

## Appendix C — Meeting notes received 2026-11-04

```
Correlation between ligand size and pocket size
APO is the baseline and do relative according to it, volume normalized to 100%.
Go back to standard natural structure that is binding to estrogen receptor and calculate ligand volume and normalize to 100%.
And then you can calculate the percentage of filling.

Different affinities since of different assays (in vitro, ...) and create different IC50, ... There are a lot of different systems and you need to think how to compare based on different systems.

6 different action class types and based on their DOIs, show in the plot. Somehow try to show that.

Prof is saying about superagonists in vdr but scientists didn't find differences in H12 between simple agonists.

Can antagonist kill the receptor by heavily distabilizing them, the H12. Good to show the picture of H12 and the ligand and compare between several class types. This is very important. Is difference on structure linked to functional something.

One idea: get the SMILES, InChiKey from api and using chembl api to get info.

Prof asked to predict affinities based on folding but he isn't sure whether we can finish.
```

Clarifications (team, 2026-11-04):
- The teammate's action-class labels are not finished; the professor said the classes can be obtained from SMILES and the ChEMBL API. The 6 labels include agonist and SERM, and come from the papers (DOIs).
- About 7 apo structures are expected.
- Affinity prediction from folding is not part of this run.
