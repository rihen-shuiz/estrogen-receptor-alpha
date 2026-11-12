#!/usr/bin/env bash
# Rebuild every result table, figure and checked number from data/ in pipeline order.
#   bash reproduce.sh                 # uses python3 on PATH; PYTHON=/path/to/python bash reproduce.sh to choose
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
if [[ "${PY}" == */* ]]; then PY="$(cd "$(dirname "${PY}")" && pwd)/$(basename "${PY}")"; fi   # absolute path
export PYTHONWARNINGS=ignore MPLBACKEND=Agg
for stage in pipeline/[0-9][0-9]_*.py; do
  echo ">>> ${stage}"
  "$PY" -B "${stage}" > /dev/null
done
echo ">>> figures"
(cd pipeline/figures && for fig in fig_*.py; do "$PY" -B "${fig}" > /dev/null; done)
echo ">>> evidence: copy source tables, redraw evidence figures, re-check every number"
"$PY" -B evidence/rebuild.py --refresh | tail -1
