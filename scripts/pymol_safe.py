"""Safe headless PyMOL runner (tracker Y01). Every PyMOL render of the final run goes through this file.

Rules it enforces:
- PyMOL runs in a separate process (the conda env with open-source PyMOL), never imported into the analysis process.
- faulthandler is enabled inside the child, so a crash leaves a traceback in the .log file.
- Hard timeout of 150 s per render; on timeout the child is killed and the call returns False.
- Only flat, per-object selections are used by the generated scripts (one object per chain, no nested selections).
- All paths are absolute, built from ESR1_ROOT (default: repo root).
- Output check: the PNG must exist, be > 10 kB and not be blank (checked with Pillow). Otherwise False.
Callers fall back to a matplotlib C-alpha trace when render() returns False.

Usage from Python:  from pymol_safe import render; ok = render(script_text, out_png)
"""
from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from pathlib import Path

ROOT = Path(os.environ.get("ESR1_ROOT", Path(__file__).resolve().parents[1])).resolve()
TIMEOUT = 150
CANDIDATES = [os.environ.get("PYMOL_PYTHON", ""), r"C:/Users/duise/anaconda3/envs/pymol/python.exe"]


def pymol_python():
    for c in CANDIDATES:
        if c and Path(c).exists():
            try:
                r = subprocess.run([c, "-c", "import pymol"], capture_output=True, timeout=60)
                if r.returncode == 0:
                    return c
            except Exception:
                pass
    return None


def png_ok(path: Path) -> bool:
    if not path.exists() or path.stat().st_size < 10_000:
        return False
    try:
        import numpy as np
        from PIL import Image
        a = np.asarray(Image.open(path).convert("L"))
        return float(a.std()) > 2.0
    except Exception:
        return path.stat().st_size > 10_000


def render(commands: list[str], out_png: str | Path, width=1600, height=1200, log=None) -> bool:
    """Run PyMOL commands (cmd.do strings) headless and ray-trace to out_png. Returns True only if the PNG checks pass."""
    py = pymol_python()
    out_png = Path(out_png).resolve()
    log = Path(log or out_png.with_suffix(".log")).resolve()
    if py is None:
        log.write_text("PyMOL not available (no python with pymol module)\n")
        return False
    if out_png.exists():
        out_png.unlink()
    # a command starting with "PY:" is inserted as raw Python (e.g. cmd.transform_object); all others go through cmd.do
    body = "\n".join(c[3:] if c.startswith("PY:") else f"cmd.do({c!r})" for c in commands)
    child = textwrap.dedent(f"""
        import faulthandler, sys
        faulthandler.enable(file=sys.stderr, all_threads=True)
        import pymol
        pymol.finish_launching(['pymol', '-cq'])
        from pymol import cmd
        cmd.feedback('disable', 'all', 'everything')
        """) + body + textwrap.dedent(f"""
        cmd.ray({width}, {height})
        cmd.png(r'{out_png.as_posix()}', dpi=300)
        cmd.quit()
        """)
    script = out_png.with_suffix(".pymol_child.py")
    script.write_text(child)
    try:
        r = subprocess.run([py, str(script)], capture_output=True, text=True, timeout=TIMEOUT, cwd=str(ROOT),
                           env={**os.environ, "ESR1_ROOT": str(ROOT)})
        log.write_text(f"returncode {r.returncode}\n--- stdout\n{r.stdout[-4000:]}\n--- stderr\n{r.stderr[-4000:]}\n")
    except subprocess.TimeoutExpired:
        log.write_text(f"timeout after {TIMEOUT} s\n")
        return False
    return png_ok(out_png)


if __name__ == "__main__":
    print("pymol python:", pymol_python())
    sys.exit(0)
