# TASK LANGUAGE (tasks v1)

File: `esr1.tasks`. One line = one statement. Plain text, read by humans and by `tl.py`.
The task file says WHAT must exist and HOW to check it. It never says what is done; that lives only in `esr1.log`.

## Grammar

```
file     := "LANG tasks v1" { comment | task }
comment  := "#" text                       (whole line only)
task     := header { field }
header   := "T" ID SIZE KIND '"' title '"'          (starts in column 0)
field    := two spaces, KEY, space, value           (indented)
```

Header words:

| word  | values | meaning |
|-------|--------|---------|
| ID    | letter + 2 digits, e.g. `S04` | unique. Letter = group (below) |
| SIZE  | `XS` `S` `M` | XS < 30 min, S <= 2 h, M <= half a day. No larger size exists: split instead |
| KIND  | `code` `run` `fig` `text` `decide` `manual` | `code` write a script; `run` execute something and keep the output; `fig` make a figure; `text` write prose from result files; `decide` / `manual` need a human, Claude must not mark them DONE |
| title | quoted | imperative, short |

Fields (each on its own line, in any order; `out`, `chk`, `orig` may repeat):

| key   | value | meaning |
|-------|-------|---------|
| help  | one sentence | plain-language purpose. Required |
| meta  | `wk=N who=X ph=P` | wk = project week from the description; who = role A/B/C/D/ALL; ph = phase. Required |
| needs | IDs | tasks that must be DONE or SKIP first |
| tool  | comma list | libraries / programs expected (hint) |
| out   | path or glob | what the task produces (used as log evidence) |
| chk   | TYPE args | one machine check. At least one per task. See below |
| orig  | text | where the same thing exists in the ORIGINAL project folder. Hint only: port it, then re-check |
| src   | text | where the requirement comes from in the project description |

Groups (ID letter): `Q` data + QC, `S` structure, `R` redundancy, `C` chemistry, `L` chemistry-structure linkage, `E` ESR1/ESR2 comparison, `X` extensions (ESR2, pockets, ChEMBL), `T` statistics, `F` PyMOL figures, `W` repo + manuscript.
Phases: `qc` `core` `ext` `fig` `stat` `wrap`. Rule from the project description: finish `qc` and `core` before `ext`; `needs` enforces it.
Roles: A data + QC, B structure + redundancy, C chemistry, D comparison + figures + manuscript.

## chk types

Paths are relative to the repo root; `**` matches any depth; when a glob matches several files the first (sorted) is used.

| chk | args | passes when |
|-----|------|-------------|
| `file` | path | file exists and is not empty |
| `glob` | pattern `>=N` | at least N files match (N default 1) |
| `lines` | path `>=N` | at least N non-blank lines |
| `csv` | path `rows>=N` `cols=a,b,c` | at least N data rows and all named columns present |
| `npy` | path `square>=N` | numpy array is 2-D, square, side >= N |
| `grep` | path-or-glob regex | regex matches somewhere in any matching file |
| `cmd` | shell command | exit code 0 (runs in repo root, 300 s limit; only with `--cmd`) |
| `man` | sentence | a human must confirm; never passes by itself |

A task is DONE only when every non-`man` chk passes, and, if it has a `man` chk, a human logged it.

## Example

```
T S04 S run "Pairwise Calpha RMSD matrix (ESR1) - save to disk"
  help N x N RMSD over common LBD Calpha positions after the core fit; save with np.save.
  meta wk=3 who=B ph=core
  needs S03
  tool numpy
  out results/rmsd_entry.npy
  chk npy results/**/rmsd*.npy square>=300
  orig esr_analysis/calculate_rmsd.py
```

## Rules for editing the task file

- Only a human edits `esr1.tasks`. Exception: Claude may correct an `out` or `chk` PATH when the correct output exists under another name, and must then log `NOTE` (see LOG_LANG.md) naming the task, old path and new path. Claude never loosens a threshold (`rows>=`, `>=N`, `cols=`) or deletes a `chk`.
- New work found during the project: add a new task (next free number in its group). Do not enlarge an existing one.
- Counts in `chk` (e.g. 351 selected entries) are floors chosen from earlier results; they are not claims about the data.
