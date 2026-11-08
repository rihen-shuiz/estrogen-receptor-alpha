"""SV00: record the interpreter and package versions used by the SV00 scripts in results/sv/SV00/env.txt (rule 1.5).
Run: .venv/bin/python -I scripts/sv/SV00_env.py
"""
import datetime as dt
import platform
import sys
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKAGES = ["numpy", "pandas", "scipy", "scikit-learn", "gemmi", "rdkit", "biopython", "requests"]

lines = [f"# SV00 environment, {dt.datetime.now().astimezone().isoformat(timespec='minutes')}",
         f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}"]
for name in PACKAGES:
    try:
        lines.append(f"{name} {metadata.version(name)}")
    except metadata.PackageNotFoundError:
        lines.append(f"{name} not installed")
lines += ["nothing installed, upgraded or downgraded by SV00",
          "network: files.rcsb.org (mmCIF), data.rcsb.org/graphql (Data API), www.ebi.ac.uk/pdbe/api (SIFTS)"]
out = ROOT / "results/sv/SV00/env.txt"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("\n".join(lines) + "\n")
print("\n".join(lines))
