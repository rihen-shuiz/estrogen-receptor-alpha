#!/usr/bin/env bash
# Reproducible runner for Stage 12 (curated ligand action analysis, 20 figures, 10 PyMOL H11/H12 renders).
# Usage:
#   bash reproduce_curated.sh --dataset esr1_dataset_final.csv
#   bash reproduce_curated.sh --dataset esr1_dataset_final.csv --download-pymol
set -euo pipefail
cd "$(dirname "$0")"

DATASET="esr1_dataset_final.csv"
PYMOL_MODE="--render-only"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset)
      DATASET="$2"
      shift 2
      ;;
    --download-pymol)
      PYMOL_MODE=""
      shift
      ;;
    --skip-pymol)
      PYMOL_MODE="SKIP"
      shift
      ;;
    --render-only)
      PYMOL_MODE="--render-only"
      shift
      ;;
    *)
      echo "Unknown arg: $1" >&2
      exit 1
      ;;
  esac
done

PY="${PYTHON:-.venv/bin/python}"
if [[ ! -x "$PY" ]]; then PY="python3"; fi
export PYTHONWARNINGS=ignore MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/mpl}"

echo ">>> [1/6] Join & audit ($DATASET)"
"$PY" -B pipeline/curated_action/01_join_and_audit.py --dataset "$DATASET"

echo ">>> [2/6] Redundancy, set cover, medoids & QC sensitivity"
"$PY" -B pipeline/curated_action/02_redundancy_and_qc.py

echo ">>> [3/6] Chemical space (RDKit ECFP4 & MACCS)"
"$PY" -B pipeline/curated_action/03_chem_space.py

echo ">>> [4/6] Statistical effects & priority figures (F01-F12)"
"$PY" -B pipeline/curated_action/04a_effects_and_priority_figs.py

echo ">>> [5/6] Appendix figures (F13-F20) & manifest"
"$PY" -B pipeline/curated_action/04b_appendix_figs.py

if [[ "$PYMOL_MODE" != "SKIP" ]]; then
  echo ">>> [6/6] PyMOL 10-ligand H11/H12 visualizations ($PYMOL_MODE)"
  "$PY" -B pipeline/curated_action/05_pymol_10_ligands.py $PYMOL_MODE
fi

echo ">>> Stage 12 curated action pipeline finished cleanly."
