# ESR1 project tracker (paste into the repo's CLAUDE.md)

Files in `tracker/`: `esr1.tasks` (what must exist), `esr1.log` (what is done), `tl.py` (the only way to read state and write the log), `TASK_LANG.md`, `LOG_LANG.md`. Read both language files once before working.

## Session start
1. `python tracker/tl.py validate` must print 0 errors. Stop and report if not.
2. If `esr1.log` has no task lines yet, run `python tracker/tl.py audit --write`. This logs every task whose outputs already pass the checks (who=audit).
3. Spot-check 5 of the audit DONE lines: open the evidence file, confirm it is real and has the content the task title promises (not just the right shape). If one is wrong, `log ID REOPEN --note ...`.
4. `python tracker/tl.py status`.

## Work cycle (one task at a time)
1. Take the first task under "READY for Claude". Ignore "NEEDS HUMAN" tasks.
2. `python tracker/tl.py show ID`. The card is the whole spec.
3. Do the work. If the card has `orig`, the original project folder may already contain a working version: read it, port it, run it again here. Do not assume it is correct or that its output is current.
4. `python tracker/tl.py log ID DONE --ev <files you produced>`. If it refuses, fix the failing check, or log `PART` or `BLOCK` with a precise `--note`.
5. Next task. Run `status` again every few tasks.

## Hard rules
- Never edit `esr1.log` by hand. Never log DONE without `tl.py` accepting it.
- Never edit `esr1.tasks`, except correcting an `out`/`chk` path when the output exists under another name; then also append `tl.py log - NOTE --note "spec-path-fix ID: old -> new"`. Never lower a threshold or delete a check.
- `decide` and `manual` tasks: do the preparation, log `MANUAL` with what the human must do, move on. Do not log DONE for them.
- Do not run `git commit`, `git push`, or delete files, unless the task card says so.
- No network or missing input: log `BLOCK`, do not fake data. Never invent numbers, sources, ligand labels or p-values. A label without a real source is `unknown`.
- Do not recompute a saved matrix (`rmsd_*.npy`, `tanimoto_*.npy`); load it. Do not drop failing-QC entries from files; flag them. Do not relax QC cutoffs silently.
- Keep each change inside the task's `out` paths. No refactors of other code.

## Finish message after each cycle (max 3 lines)
Task id, state logged, evidence path. Nothing else.
