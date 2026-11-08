"""FINAL (RUN_SV1028): write docs/sv/RUN_REPORT.md from RUN_STATUS.md plus the run-level notes."""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
status = (ROOT / "docs/sv/RUN_STATUS.md").read_text(encoding="utf-8").splitlines()
rows = [l for l in status if l.startswith("| ") and not l.startswith("| step") and not l.startswith("|---")]
line1 = (
    "Not done / substituted: W17 logged MANUAL, not DONE (tl.py lets only a human log DONE for a manual task); "
    "packages installed from public PyPI by the user, not through the Socket proxy (unreachable without VPN); "
    "mkdssp unavailable -> pydssp fallback in SV02 ('H' = any helix, no 3_10/pi fractions), inherited by SV03, SV04, SV07; "
    "PyMOL unavailable -> SV13 Part C matplotlib C-alpha fallback; no SVG renderer -> VIS1 method (b) replaced by PNG gutters, 9 images not_cut; "
    "SV00 used the download function of scripts/fetch_data.py instead of the whole script (it writes outside rule 1.2), original download date "
    "taken from the Q03 log line, pH given only as text left NaN, no fresh RCSB search for newer entries; "
    "SV06 panel B uses the Y10 complete-H12 medoids (the S12 cl4 medoid has no modelled H12); SV04 adds an undeclared figure panel (b); "
    "SV07 'sparse-sample effect' verdict rule degenerate for cluster-mean spikes (BLOCKED part); SV01b step 8 not run (verdict not mixed), run in SV01c; "
    "SV12 teammate labels absent and route (c) failed its gate -> provisional ChEMBL-mechanism classes for 51 of 362 entry x ligand rows; "
    "SV13: all class comparisons not testable (coverage gate), SERM representative without modelled H12; SV15 class comparison not testable, pairing "
    "rule corrected to one dpX per assay after the first run (disclosed); SV14 pocket 100 % = wild-type E2 pocket (fallback: 3 apo entries < 5); "
    "SV08/SV08b: L536S mechanism and ERalpha volume methods have no openable primary source; 0 chains excluded as remediated; "
    "read-only check: 0 of 2576 snapshot files changed or deleted."
)
out = [line1, "", "# RUN_REPORT - RUN_SV1028 (supervisor requests 2026-10-28, addition 2026-11-04)", "",
       "Branch `sv1028-run` (from `tracker-tasks` @ 8828a30); nothing committed. Per-step details: results/sv/<ID>/SUMMARY.md.", "",
       "| step | status | key number or verdict | main files |", "|---|---|---|---|"] + rows + [
    "",
    "## Decisions for the team",
    "",
    "1. **Sign off the H12 window.** SV02 proposes 538-544 from the pydssp consensus (cl7 538-545, cl4 537-544); primary analyses used 538-547. "
    "On 538-544 the SV01 verdict changes from artifact to mixed for completeness at chain level only (the event becomes rare: 10 cl7 and 19 cl4 "
    "incomplete chains); mobility stays artifact (results/sv/SV01c/side_by_side.csv). A re-run with mkdssp would separate 3_10/pi turns before sign-off.",
    "2. **Adopt core v4 or keep v3.** v4 = 128 positions (adds H1, H9; drops 25 helix-end positions), Jaccard 0.64 with v3; RMSD Spearman 0.9975, "
    "3 chains change cluster at 1.5 A, cl4 sub-clusters regroup at 1.0 A (ARI 0.817) (results/sv/SV05/).",
    "3. **Review the `visuals/` selection** (77 files, visuals/INDEX.md; exclusions with reasons in results/sv/VIS2/selection.csv).",
    "4. **Replace the provisional ChEMBL classes with the teammate's labels when they are ready, then re-run SV13-SV15** "
    "(only agonist and SERM pass the 5-entry gate now; the SV13 model code was checked on pseudo-labels).",
    "5. **Confirm W17** (presentation took place): `python tracker/tl.py log W17 DONE --who human --ev docs/supervisor_feedback.md`.",
    "6. **Keep in mind:** docs/RESULTS_SUMMARY.md is rendered from a template by scripts/collect_numbers_ext.py; a re-render drops the new SV "
    "section - re-run scripts/sv/SV09_docs.py and SV09_number_check.py after it.",
]
(ROOT / "docs/sv/RUN_REPORT.md").write_text("\n".join(out) + "\n", encoding="utf-8")
print("RUN_REPORT rows:", len(rows))
