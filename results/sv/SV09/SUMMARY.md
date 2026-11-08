Not done / substituted: (1) "one figure per step" is read as naming one existing figure per step - no new figure was drawn (rule 1.1); SV00, SV01c, SV08, SV08b, SV10 and SV12 made no figure and say "Figure: none"; the "SV01 forest plot" is the pair of SV01b forests (one per outcome, as the SV01 figure text asks), and the SV02 residue-537 Ramachandran plot is a main figure as well, because the SV02 step text says so; (2) the text is kept as a tagged source (docs/sv/SV09_section.src.md, docs/sv/SV09_limitations.src.md: every number carries its file and column/key) and rendered by scripts/sv/SV09_docs.py, which appends it after the untouched earlier text and keeps plain copies (docs/sv/SV09_section.md, docs/sv/SV09_limitations.md) - supporting files not named in the step; (3) docs/RESULTS_SUMMARY.md is rendered from docs/templates/RESULTS_SUMMARY.md by scripts/collect_numbers_ext.py (template not edited: outside the step's write scope), so a re-render drops the new section; re-run scripts/sv/SV09_docs.py then SV09_number_check.py (noted in docs/EXPLORATION.md, "Suggested, not done (SV09)"); (4) number_check.csv has one extra leading column, doc, because the numbers of the new limitations were checked too; 87 of the 439 numbers are matched by a regular expression in a saved MD file (step SUMMARY/verdict files, docs/sv/literature*.md, docs/supervisor_feedback.md for dates and quoted requests) where no CSV/JSON holds them, and dates are compared as strings; the check verifies values, not wording - for residue ranges such as 331-340 the two end points are checked against the step's per-position verdict rows; (5) this SUMMARY.md was returned as text because the subagent's Write tool refused report files; the orchestrator saved it.

# SV09 - Integration

**Result:** new section "Supervisor requests 2026-10-28" appended at the end of docs/RESULTS_SUMMARY.md (earlier text byte-identical: the HEAD version is an exact prefix of the new file) and 12 new limitations appended to docs/limitations.md (same check). Number check: **439 numbers, 439 match, 0 mismatch** (results/sv/SV09/number_check.csv; 72 source files).

## What the section contains

- Intro: inputs, fixed choices (H12 538-547, core v3, S08 at 1.5 A), wording rule (H12 disorder / H12 mobility, never receptor activity or degradation; ligands by chemistry in SV00-SV11; SV12 classes = provisional ChEMBL-mechanism labels).
- Main figures for the next meeting (table, 10 files): SV01b forest plots (completeness, mobility), SV02 helix track, SV02 residue-537 Ramachandran plot, SV04 h12_map_partial (= SV03 H12 map with the SV04 overlay), SV06 outlier panels, SV07 spike profile, SV13 per-class H12 panels, SV13 forest plot, SV14 pocket vs ligand. Every other figure is listed as backup (14 files + the other SV02 Ramachandran plots).
- One paragraph per step (18): SV00, SV01a, SV01b, SV01c, SV02, SV03, SV04, SV05, SV06, SV07, SV08, SV10 (one sentence + table link), SV11 (Part A paragraph, Part B one sentence), then "Addition: notes received 2026-11-04": SV12, SV13, SV14, SV15, SV08b. Each quotes the request where one exists, gives the verdict and key numbers, the caveats from SUMMARY line 1 that change how a number reads, and its figure (main / backup / none).
- Caveats carried in:
  - pydssp fallback (SV02, SV03);
  - SV01a: cl7 near tie, and 14 incomplete chains for 16 predictors;
  - SV01b: lattice terms counted on modelled atoms, so partly mechanical for completeness (cl4 24.3 vs 9.1 symmetry-contact atoms, 12.3 of them on 546-550; 8.22 vs 6.66 on 538-545);
  - SV01c: verdict depends on the window; the entry level rests on 4 / 6 incomplete entries;
  - SV03: 4 of 37 hits on <= 5 entries; the 7-C-alpha axis error is 4.5 deg;
  - SV04: not independent of S08;
  - SV06: medoid substitution;
  - SV07: the sparse-sample rule is BLOCKED;
  - SV08: L536S gap;
  - SV12: teammate file absent; DOI route failed (1 of 2);
  - SV13: not testable (0 of 22 rows); SERM representative has no H12; matplotlib fallback;
  - SV14: pocket reference fallback (3 apo entries < 5); bimodal E2 pocket; box edge in 183 analysed chains;
  - SV15: first-run pairing correction (q 0.033 -> 0.057);
  - SV08b: volumes second-hand, without a method.
- Literature enters only from opened sources, as recorded in docs/sv/literature.md and docs/sv/literature_class.md (Li et al. 2011 for the unmodelled 331-340 / 462-463 loops; the SV08b answers with their counts).

## Key numbers in the section (all in number_check.csv)

- SV01:
  - SV01a mixed (chain cl7 Q+X 0.224 vs B 0.225; entry artifact).
  - SV01b artifact (chain cl7 0.414 vs 0.225, cl4 0.523 vs 0.079; mobility cl7 0.184 vs 0.076, cl4 0.210 vs 0.143).
  - SV01c mixed on 538-544 (cl7 0.130 vs 0.386); clash test 3 of 31 blocked at both positions, 0 at one only.
- SV02-SV07:
  - SV02: 8 consensus helices, H12 proposal 538-544, psi 537 +156.6 / -46.8 deg, overlap 0.
  - SV03: 12 of 16 descriptors separate the clusters; 37 within-cluster hits.
  - SV04: ARI 0.954 / 0.931.
  - SV05: core v4 128 positions, Jaccard 0.64, ARI 0.990 at 1.5 A.
  - SV06: 2 of 14 main-cluster outliers in cl7 vs 8.0 expected (q 0.0023).
  - SV07: 19 spikes, 14 genuine, 5 driven by <= 3 entries.
- SV11-SV15:
  - SV11: AF-2 groove pocket cl7 0.53 vs cl4 0.02; steroid entries 41/217 vs 5/124 (q 0.00011).
  - SV12: 51 of 362 rows labelled.
  - SV13: 0 of 22 testable.
  - SV14: rho 0.44 / 0.47; filling median 33.8 %.
  - SV15: 14 ligands with dpX; cl7 vs cl4 q 0.15 / 0.057.

## Number check (scripts/sv/SV09_number_check.py)

- Tokeniser: every number in the new section and in the new limitations, outside code spans; identifiers (H12, cl7, SV01b, Y537S, 5TM1, AF-2, H10/11, A^3) are not numbers. Each doc number is aligned with the rendered source, so a number without a tag fails; each tag's value is read from its file (CSV column + row filters, JSON key path, or MD regex) and compared at the displayed precision.
- RESULTS_SUMMARY section: 397 numbers (CSV 190, JSON 125, MD 82); limitations: 42 (CSV 6, JSON 31, MD 5). Match 439 / 439.
- Negative test (not saved): changing one value in a copy of the doc (0.414 -> 0.415) and adding an untagged number both fail the check.

## Verdict

Integration done: every DONE step (SV00, SV01a, SV01b, SV01c, SV02-SV08, SV10-SV15, SV08b) has one paragraph and one figure entry in docs/RESULTS_SUMMARY.md; all numbers match their files; new limitations added to docs/limitations.md.

## Files

- docs/RESULTS_SUMMARY.md (edited: section appended, 21.9 kB added; earlier 10.7 kB unchanged)
- docs/limitations.md (edited: 12 limitations appended, 5.6 kB added; earlier 1.2 kB unchanged)
- docs/sv/SV09_section.src.md (43.2 kB), docs/sv/SV09_section.md (21.9 kB)
- docs/sv/SV09_limitations.src.md (7.4 kB), docs/sv/SV09_limitations.md (5.6 kB)
- docs/EXPLORATION.md: appended "Suggested, not done (SV09, 2026-11-08)" (rule 1.1; existing text unchanged)
- results/sv/SV09/number_check.csv (98.8 kB, 439 rows)
- results/sv/SV09/env.txt (0.5 kB)
- results/sv/SV09/SUMMARY.md (this file)
- scripts/sv/SV09_docs.py (1.9 kB), scripts/sv/SV09_number_check.py (16.0 kB)

Run order: .venv/bin/python -B scripts/sv/SV09_docs.py; .venv/bin/python -B scripts/sv/SV09_number_check.py
