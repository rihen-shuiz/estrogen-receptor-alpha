# SV02 - data-derived H12 window (proposal only)

Status: **proposed, not adopted.** The H12 window for all primary analyses stays 538-547 (run file section 2);
a human signs off any change. SV01c tests the proposal as a sensitivity analysis.

Method: pydssp 0.9.1 (fallback for mkdssp, not installed); 'H' = any helix (3-, 4- and 5-turn helices together); [DEFAULT] consensus helix = >= 4 consecutive positions with H in >= 80 % of the observing chains, low-support positions (observed in < 50 % of chains) never a boundary.

**Data-derived H12 window (global consensus, 788 chains): 538-544** (7 positions; run before trimming 538-544; low-support positions inside: 0; min H fraction 0.933, min fraction observing 0.965).

- cl7 consensus H12: 538-545
- cl4 consensus H12: 537-544

## Comparison with the historical windows

| window | positions | shared with proposal | only in window | only in proposal | Jaccard | start shift | end shift |
|---|---|---|---|---|---|---|---|
| 538-547 (settings h12_range, primary) | 538-547 | 7 | [545, 546, 547] | - | 0.70 | +0 | -3 |
| 537-546 (helix table AC4, 7NFB struct_conf) | 537-546 | 7 | [537, 545, 546] | - | 0.70 | +1 | -2 |
| 536-545 | 536-545 | 7 | [536, 537, 545] | - | 0.70 | +2 | -1 |

The proposal equals none of the three historical windows.
It differs from the primary window 538-547, so SV01c re-runs SV01b on it (run file SV01c step 9).

## Per-position support, 530-550

| position | helix table | observing global | H global | H cl7 | H cl4 | obs cl7 | obs cl4 |
|---|---|---|---|---|---|---|---|
| 530 | H10/11 | 0.665 | 0.635 | 0.828 | 0.009 | 0.888 | 0.353 |
| 531 | H10/11 | 0.645 | 0.364 | 0.429 | 0.181 | 0.858 | 0.353 |
| 532 | loop | 0.675 | 0.051 | 0.016 | 0.146 | 0.849 | 0.438 |
| 533 | loop | 0.746 | 0.037 | 0.003 | 0.100 | 0.822 | 0.641 |
| 534 | loop | 0.821 | 0.000 | 0.000 | 0.000 | 0.858 | 0.772 |
| 535 | loop | 0.892 | 0.003 | 0.000 | 0.000 | 0.917 | 0.863 |
| 536 | loop | 0.945 | 0.077 | 0.000 | 0.167 | 0.964 | 0.927 |
| 537 | H12 | 0.964 | 0.408 | 0.000 | 0.965 | 0.980 | 0.948 |
| 538 | H12 | 0.964 | 0.992 | 0.998 | 0.984 | 0.978 | 0.951 |
| 539 | H12 | 0.966 | 0.996 | 1.000 | 0.990 | 0.978 | 0.954 |
| 540 | H12 | 0.970 | 0.995 | 0.995 | 0.994 | 0.982 | 0.957 |
| 541 | H12 | 0.968 | 0.997 | 0.995 | 1.000 | 0.982 | 0.954 |
| 542 | H12 | 0.971 | 0.993 | 0.993 | 0.997 | 0.984 | 0.957 |
| 543 | H12 | 0.970 | 0.984 | 0.986 | 0.987 | 0.982 | 0.957 |
| 544 | H12 | 0.966 | 0.933 | 0.986 | 0.882 | 0.980 | 0.951 |
| 545 | H12 | 0.935 | 0.693 | 0.924 | 0.364 | 0.978 | 0.894 |
| 546 | H12 | 0.853 | 0.246 | 0.325 | 0.099 | 0.975 | 0.705 |
| 547 | loop | 0.702 | 0.136 | 0.135 | 0.142 | 0.969 | 0.365 |
| 548 | loop | 0.607 | 0.065 | 0.052 | 0.121 | 0.867 | 0.277 |
| 549 | loop | 0.236 | 0.048 | 0.000 | 0.130 | 0.263 | 0.210 |
| 550 | loop | 0.122 | 0.000 | 0.000 | 0.000 | 0.088 | 0.173 |

Values from results/sv/SV02/sse_per_position.csv and consensus_{global,cl7,cl4}.csv.
