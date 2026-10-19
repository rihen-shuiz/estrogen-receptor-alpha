# LOG LANGUAGE (log v1)

File: `esr1.log`. Append-only. One line = one event. The log is the only record of progress.
Write lines with `tl.py log` (it checks the task first). Never edit or delete an old line; to correct, append a new line.

## Grammar

```
file  := "LANG log v1" { line }
line  := TIME " | " ID " | " STATE { " | " KEY "=" VALUE }
TIME  := ISO 8601 with offset, minute precision, e.g. 2026-10-14T16:42+05:00
ID    := task id (S04) or "-" for a general note
```

VALUE is one line of text and never contains ` | `.

Keys: `who` (required, except on NOTE), `chk`, `ev`, `ref`, `note`.

| key  | meaning |
|------|---------|
| who  | `haiku`, `opus`, `sonnet`, `audit` (written by `tl.py audit`), or `human` |
| chk  | `passed/total` of the automatic checks, written by `tl.py` |
| ev   | evidence: comma-separated paths that were produced or inspected. Required on DONE |
| ref  | git commit hash if the work was committed |
| note | free text. Required on every state except DONE and REOPEN |

## States

| state | use it when | needs |
|-------|-------------|-------|
| `WIP` | work started and was interrupted | note: where it stopped |
| `PART` | some checks pass, some fail | note: which checks fail and why |
| `BLOCK` | cannot continue without something outside the task (missing input, no network, spec path wrong, needs a decision) | note: exactly what is needed |
| `MANUAL` | Claude finished its part and a human must confirm or do the rest | note: what the human must do |
| `DONE` | every automatic chk passes | `chk=n/n`, `ev=` |
| `SKIP` | deliberately not doing it (human decided) | note: reason and who decided |
| `REOPEN` | a DONE task is found wrong or stale | note: what is wrong |
| `NOTE` | general remark, id is `-` | note |

Current state of a task = its last log line (`REOPEN` returns it to TODO; no line = TODO).
`DONE` by `who=audit` means only "the files exist and pass the shape checks"; it is not a review of the content.

## Rules

1. DONE is written only by `tl.py log ID DONE --ev ...`, which runs the checks and refuses when any fail. Never write a DONE line by hand.
2. Tasks with a `man` check: Claude logs `MANUAL`. Only the human logs DONE (`--who human`).
3. Evidence must be real files that exist. No evidence, no DONE.
4. Every number written into a `text` task must come from a file in `ev`.
5. One task per work cycle; log it before starting the next.
6. Failure is information: a `BLOCK` or `PART` line with a precise note is a good result. A DONE that does not hold is the worst result.

## Examples

```
2026-10-14T16:42+05:00 | Q06 | DONE | who=haiku | chk=1/1 | ev=results/entries.csv
2026-10-14T16:51+05:00 | Q07 | BLOCK | who=haiku | note=needs Q04 output; results/lbd_chains.csv missing
2026-10-14T17:03+05:00 | S04 | PART | who=haiku | chk=0/1 | ev=results/rmsd_entry.npy | note=chk npy fails: matrix is 280x280, need square>=300; 20 chains missing from alignment
2026-10-14T17:10+05:00 | - | NOTE | who=haiku | note=spec-path-fix Q04: results/lbd_chains.csv -> results/lbd_chains_esr1.csv
2026-10-14T18:00+05:00 | Q10 | DONE | who=human | chk=1/1 | ev=results/ligand_action_labels.csv | note=sources checked
```

## Commands (`tl.py`)

```
python tracker/tl.py validate                  check both files
python tracker/tl.py status                    progress board, next ready tasks, tasks needing a human, blocked tasks
python tracker/tl.py show S04                  task card + its history
python tracker/tl.py check S04 [--cmd]         run the checks only
python tracker/tl.py audit [--write] [--cmd]   check every TODO task; --write logs DONE/PART as who=audit
python tracker/tl.py log S04 DONE --ev results/rmsd_entry.npy
python tracker/tl.py log S04 PART --note "what is missing"
python tracker/tl.py log Q10 MANUAL --note "labels drafted; human must verify sources"
```

Run from the repo root (checks use the current directory), or set `ESR_ROOT`.
