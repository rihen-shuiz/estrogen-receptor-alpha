# QC summary

Source: `results/entries.csv` (QC flags, Q05 and Q06) and `results/lbd_chains.csv` (chain eligibility, Q04).

Cutoffs: resolution fails above 2.5 A, R-free fails above 0.30, mean B fails above 50 A^2. A missing value is `unknown`, never `pass`. The combined flag `qc_pass` is `fail` if any flag fails, `pass` if all three pass, and `unknown` otherwise.

Sets:

- All entries: 526 deposited ESR1 and ESR2 X-ray entries (Q01, Q02).
- Selected ESR1 set: 351 entries with at least one ESR1 receptor chain marked `eligible_historical=True` (180 LBD CA and 60 core CA screen, unchanged from the original project). The 39 ESR2 validation entries are outside this set.

| Flag | Rule | All entries (526): pass / fail / unknown | Selected ESR1 set (351): pass / fail / unknown |
|---|---|---|---|
| Resolution | 2.5 A | 481 / 45 / 0 | 315 / 36 / 0 |
| R-free | 0.30 | 513 / 11 / 2 | 343 / 8 / 0 |
| Mean B | 50 A^2 | 424 / 102 / 0 | 267 / 84 / 0 |
| Combined | all three pass | 397 / 128 / 1 | 247 / 104 / 0 |

Failing entries stay in every file. They are flagged, not deleted.

The two R-free `unknown` entries and the one combined `unknown` entry are not in the selected ESR1 set.
