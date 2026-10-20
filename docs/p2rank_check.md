# P2Rank check (task X07)

- Machine: Windows 11, 16 CPUs. Java: OpenJDK 21.0.11 LTS (`java -version`).
- P2Rank 2.5.1, the official release from https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz (275,625,956 bytes). It is unpacked to `C:\Users\duise\tools\p2rank_2.5.1`, outside the repository.
- Test run: `prank.bat predict -f test_data\1fbl.pdb -o <scratch>`. It finished successfully in 17.2 s and wrote `1fbl.pdb_predictions.csv` and `1fbl.pdb_residues.csv`. The P2Rank JVM was capped at 2 GiB of memory.

## Running the tracker checks on Windows

On Windows, `tl.py` runs `chk cmd` lines through `cmd.exe`, which has no POSIX `command -v`. A shim, `C:\Users\duise\tools\bin\command.cmd`, maps `command -v NAME` to `where NAME`:

- it exits 0 only when NAME resolves on PATH;
- it exits 1 otherwise (checked with a name that does not exist).

To reproduce the checks:

```
set PATH=C:\Users\duise\tools\bin;C:\Users\duise\tools\p2rank_2.5.1;%PATH%
python tracker/tl.py check X07 --cmd
```

On Linux or macOS the shim is not needed. Put the P2Rank folder on PATH and run the same command.
