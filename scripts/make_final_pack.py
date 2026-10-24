"""Tracker Y41: final deliverable final_pack_<date>.zip with a MANIFEST.

Contents: docs/ (incl. gallery.html, interactive_embedding.html, figures/INDEX.md, templates/), figures/ (all PNG/SVG/
PDF/TXT/PML incl. figures/ext/), every CSV / JSON / small NPY under results/ (top level, results/ext/, results/rebuild/
report and tables, results/alignment/ metadata) plus the key matrices (rmsd_entry.npy, rmsd_old_core.npy,
rmsd_esr2.npy, results/ext/rmsd_*.npy), scripts/, esr_analysis/ code and config (not the bundled raw data), tracker/
(tasks, log, tl.py, language files), run_all.sh, requirements-lock.txt, requirements.txt, notebooks/.
Left out (re-downloadable or bulk): data/raw, data/clean, data/cache, results/rebuild/raw, results/rebuild/ccd,
results/p2rank per-chain files, esr_analysis/resources/data, .venv, __pycache__, old zip packs.
MANIFEST: docs/MANIFEST_final.csv (path, bytes, sha256) - written first, then zipped together with everything.
Run from the repo root: python scripts/make_final_pack.py [YYYYMMDD]
"""
import csv
import datetime as dt
import hashlib
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SKIP_DIRS = {"__pycache__", ".venv", ".git", "raw", "clean", "cache", "p2rank", "ccd", "pymol_panels"}


def wanted(p: Path) -> bool:
    rel = p.relative_to(ROOT).as_posix()
    if any(part in SKIP_DIRS for part in p.relative_to(ROOT).parts[:-1]):
        return False
    if rel.startswith("esr_analysis/resources/data") or rel.endswith(".zip") or rel.endswith(".pyc"):
        return False
    if rel.startswith("results/"):
        if p.suffix in (".csv", ".json", ".txt", ".md"):
            return p.stat().st_size < 30_000_000
        if p.suffix == ".npy":
            return p.stat().st_size < 8_000_000
        if p.suffix == ".npz":
            return "alignment" in rel
        return False
    return True


def collect():
    files = []
    for top in ("docs", "figures", "results", "scripts", "esr_analysis", "tracker", "notebooks"):
        for p in sorted((ROOT / top).rglob("*")):
            if p.is_file() and wanted(p):
                files.append(p)
    for name in ("run_all.sh", "requirements-lock.txt", "requirements.txt", "README.md", "CLAUDE.md", ".gitignore"):
        if (ROOT / name).exists():
            files.append(ROOT / name)
    files += sorted((ROOT / "figures/ext/pymol_panels").glob("*.png"))
    return files


def main():
    date = sys.argv[1] if len(sys.argv) > 1 else dt.date.today().strftime("%Y%m%d")
    files = [f for f in collect() if f.name != "MANIFEST_final.csv"]
    man = ROOT / "docs/MANIFEST_final.csv"
    with man.open("w", newline="", encoding="utf8") as fh:
        w = csv.writer(fh)
        w.writerow(["path", "bytes", "sha256"])
        for f in files:
            w.writerow([f.relative_to(ROOT).as_posix(), f.stat().st_size, hashlib.sha256(f.read_bytes()).hexdigest()])
    out = ROOT / f"final_pack_{date}.zip"
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for f in files + [man]:
            z.write(f, f.relative_to(ROOT).as_posix())
    print(out, f"{out.stat().st_size / 1e6:.1f} MB", len(files) + 1, "files")


if __name__ == "__main__":
    main()
