"""Tracker X08: P2Rank pocket prediction on every eligible ESR1 chain.

Input: data/clean/<entry>_<chain>.pdb (scripts/clean_chain_pdbs.py; receptor chain plus its ligand of interest as
HETATM). P2Rank builds pockets from protein atoms; the ligand atoms are used later (X09) for the overlap test.
Run: P2Rank `predict` on a dataset file listing every PDB (one JVM, THREADS threads), default P2Rank model,
visualizations off. Output: results/p2rank/<entry>_<chain>.pdb_predictions.csv and ..._residues.csv per chain, plus
run.log and params.txt from P2Rank, and results/p2rank/run_info.json (P2Rank version, command, chain count).
Run from the repo root: python scripts/p2rank_batch.py --prank /path/to/prank(.bat)
"""
import argparse
import json
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/p2rank"
THREADS = 16


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--prank", required=True)
    args = parser.parse_args()
    pdbs = sorted((ROOT / "data/clean").glob("*.pdb"))
    if not pdbs:
        raise SystemExit("no PDB files in data/clean; run scripts/clean_chain_pdbs.py first")
    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        dataset = Path(tmp) / "esr1_chains.ds"
        dataset.write_text("\n".join(str(p) for p in pdbs) + "\n")
        cmd = [args.prank, "predict", str(dataset), "-o", str(OUT), "-threads", str(THREADS), "-visualizations", "0"]
        result = subprocess.run(cmd, cwd=Path(args.prank).parent, capture_output=True, text=True)
    version = subprocess.run([args.prank, "-v"], cwd=Path(args.prank).parent, capture_output=True, text=True).stdout
    predictions = sorted(OUT.glob("*_predictions.csv"))
    info = {"p2rank_version": next((l.strip() for l in version.splitlines() if "P2Rank" in l), version.strip()[:80]),
            "command": " ".join(["prank", "predict", "<data/clean/*.pdb dataset>", "-threads", str(THREADS), "-visualizations", "0"]),
            "n_input_pdb": len(pdbs), "n_predictions": len(predictions), "exit_code": result.returncode}
    (OUT / "run_info.json").write_text(json.dumps(info, indent=2) + "\n")
    if result.returncode or len(predictions) != len(pdbs):
        raise SystemExit(f"P2Rank failed or incomplete: {info}\n{result.stderr[-2000:]}")
    print(f"X08: {len(predictions)} prediction files for {len(pdbs)} chains -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
