"""Orchestrator helper: update one row of docs/sv/RUN_STATUS.md (step | status | result | outputs)."""
import sys, re
from pathlib import Path
p = Path(__file__).resolve().parents[2] / "docs/sv/RUN_STATUS.md"
step, status, result, outputs = sys.argv[1:5]
lines = p.read_text(encoding="utf-8").splitlines()
row = f"| {step} | {status} | {result.replace('|', '/')} | {outputs.replace('|', '/')} |"
for i, l in enumerate(lines):
    if l.startswith(f"| {step} |"):
        lines[i] = row; break
else:
    lines.append(row)
p.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(row[:160])
