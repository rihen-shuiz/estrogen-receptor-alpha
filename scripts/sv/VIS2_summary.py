"""VIS2 SUMMARY.md and env.txt, written from the saved outputs only (rule 1.7).

Inputs: results/sv/VIS2/selection.csv (build), results/sv/VIS2/check.csv (check), visuals/INDEX.md, visuals/*.png,
figures/panels/manifest.csv. Run after VIS2_build.py and VIS2_check.py:
    .venv/bin/python -B scripts/sv/VIS2_summary.py
"""
from __future__ import annotations

import platform
import re
import sys
import time
from importlib.metadata import version
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results/sv/VIS2"
VIS = ROOT / "visuals"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import VIS2_check as chk  # noqa: E402  (INDEX parser and claim regex, reused read-only)


def size(p: Path) -> str:
    n = p.stat().st_size
    return f"{n / 1024 / 1024:.1f} MB" if n >= 1024 * 1024 else f"{n / 1024:.1f} KB"


def category(reason: str) -> str:
    r = reason.lower()
    if r.startswith("empty panel"):
        return "empty panel (no data)"
    if r.startswith("legend-only"):
        return "legend-only file"
    if r.startswith("superseded") or "september" in r:
        return "superseded or September (351-entry) version"
    if r.startswith(("same claim", "same content", "duplicate", "covered by", "colour-by map")):
        return "duplicate (another graph kept for the claim)"
    return "tied to no stated claim (QC, diagnostic, context)"


def part_of(loc: str) -> str:
    if loc.startswith("SV"):
        return "Supervisor requests (SV)"
    return {"The data": "The data", "Storyline": "Storyline"}.get(loc, "section " + loc.lstrip("§"))


def main() -> int:
    sel = pd.read_csv(RES / "selection.csv")
    checks = pd.read_csv(RES / "check.csv")
    rows = [r for r in chk.parse_index() if "bad" not in r]
    files = sorted(p for p in VIS.iterdir() if p.is_file())
    pngs = [p for p in files if p.suffix == ".png"]
    manifest = pd.read_csv(ROOT / "figures/panels/manifest.csv")

    inc = sel[sel.decision == "included"]
    exc = sel[sel.decision == "excluded"].copy()
    exc["category"] = exc.reason.map(category)
    n_leg = sum(r["nn"].endswith("b") for r in rows)
    n_graph = len(rows) - n_leg
    by_method = sel.groupby("method").size().to_dict()
    inc_method = inc.groupby("method").size().to_dict()

    # story parts (first claim location of each graph row)
    parts: dict[str, list[str]] = {}
    for r in rows:
        if r["nn"].endswith("b"):
            continue
        loc = chk.CLAIM_RE.findall(r["claim"])[0][2].strip()
        parts.setdefault(part_of(loc), []).append(r["nn"])
    n_claims = int((checks.check == "claim").sum())
    docs = sorted({m[1].strip() for r in rows for m in chk.CLAIM_RE.findall(r["claim"])})

    # SV09 main figures -> visuals numbers
    pool_src = {}
    for _, m in manifest.iterrows():
        key = m["source_path(s)"].split(";")[0] if m["method"] in ("single", "not_cut") else m["panel_path"]
        pool_src[key] = m["source_path(s)"].split(";")[0]
    mains = checks[(checks.check == "sv09_main") & (checks.item != "table found")].item.tolist()
    main_lines = []
    for f in mains:
        nns = [r["file"] for r in rows if r["source"] == f or pool_src.get(r["source"]) == f]
        n_pan = int(((manifest["source_path(s)"].str.split(";").str[0] == f) & manifest.panel_path.notna()).sum())
        how = f"{len(nns)} of {n_pan} VIS1 panels" if n_pan else "original"
        main_lines.append(f"  - {f}: {how} -> {', '.join(nns)}")

    sv13_empty = exc[exc.category == "empty panel (no data)"].key.str.extract(r"__(\w+)\.png$")[0].tolist()
    putty = exc[exc.key.str.contains("Y30_pymol_bputty__0")]
    substitutes = [r for r in rows if r["source"] in ("figures/dendrogram.png", "figures/diversity_curve.png")]
    n_fail = int((checks.result == "FAIL").sum())

    line1 = (
        "Not done / substituted: (1) 'Always include the five main SV figures (SV09)': the SV09 table lists "
        f"{len(mains)} files (with the SV13 per-class H12 panels, the SV13 forest plot and the SV14 main figure of the "
        f"2026-11-04 addition); all {len(mains)} are in visuals/, multi-graph ones as their VIS1 panels, except "
        f"{len(sv13_empty)} empty SV13 class panels ({', '.join(sv13_empty)}: 'no chain carries this class', text "
        "only), which are not copied. (2) Two claims use a keep-all (788-chain) image from figures/ that RESULTS_SUMMARY "
        "does not cite, because it prints the stated number: "
        + "; ".join(f"{r['file']} from {r['source']}" for r in substitutes)
        + " (the 7 clusters at 1.5 A, with the cited Y07 heatmap kept for the 3.49 A claim; 90% of the 1.5 A clusters "
        "at N=488, in place of the cited Y26_rarefaction panel 01, which shows the curve without the number). "
        f"(3) The Y30 B-putty pair ({len(putty)} panels, cited in section 4 'PyMOL putty') is excluded: the same claim "
        "(H12 / core-v3 B ratio 1.53 vs 1.33) is shown for all chains by Y21_helix_b; so only the Y30 medoids legend is "
        "needed (suffix b). (4) The claim check reads the documents with line breaks as single spaces and \\| as | "
        "(RESULTS_SUMMARY wraps lines; GFM tables need the escape) and checks each quote inside its cited section or "
        "SV paragraph; INDEX.md itself is not counted as a visual. (5) Supporting files not named in the run file: "
        "results/sv/VIS2/selection.csv, results/sv/VIS2/check.csv, scripts/sv/VIS2_build.py. Nothing installed."
    )

    L = [line1, "", "# VIS2 - Collect the valuable visuals into one folder", "", "**Verdict**", ""]
    L.append(f"- Pool (VIS1 manifest): {len(sel)} items - "
             + ", ".join(f"{k} {v}" for k, v in sorted(by_method.items())) + ".")
    L.append(f"- Copied into visuals/: {len(rows)} files = {n_graph} graphs + {n_leg} legend (suffix b); by pool "
             "method " + ", ".join(f"{k} {v}" for k, v in sorted(inc_method.items())) + ". Excluded: "
             f"{len(exc)}.")
    for cat, n in exc.category.value_counts().items():
        L.append(f"  - {cat}: {n}")
    L.append(f"- Story order (RESULTS_SUMMARY.md): "
             + "; ".join(f"{p} {v[0]}-{v[-1]} ({len(v)})" for p, v in parts.items()) + ".")
    L.append(f"- INDEX claims: {n_claims} verbatim quotes in {len(rows)} rows, from " + ", ".join(docs) + ".")
    by = checks.groupby(["check", "result"]).size().unstack(fill_value=0)
    L.append("- Check (scripts/sv/VIS2_check.py, results/sv/VIS2/check.csv): "
             + "; ".join(f"{c} {int(by.loc[c].get('pass', 0))} pass / {int(by.loc[c].get('FAIL', 0))} fail"
                         for c in by.index)
             + (". All pass." if n_fail == 0 else f". {n_fail} FAILED."))
    sha = checks[checks.check == "sha1"].detail.iloc[0]
    L.append(f"- Distinct SHA-1 among the {len(pngs)} PNG files: {sha}; every copy equals its source SHA-1 "
             "(check 'source').")
    L.append("")
    L.append("## Main SV figures (SV09 table) in visuals/")
    L.append("")
    L.extend(main_lines)
    L.append("")
    L.append("## Selection rule as applied (scripts/sv/VIS2_build.py)")
    L.append("")
    L.append("- One graph per claim stated in docs/RESULTS_SUMMARY.md or an SV verdict file: for a claim cited with a "
             "multi-graph figure, the panel that shows the stated number; where two graphs carry the same claim, the "
             "one that states it (reason per pool item in selection.csv, e.g. 'same claim as NN_...').")
    L.append("- Keep-all (788-chain) graphs over the September (351-entry) core-task and provisional images "
             "(figures/0*-11*, figures/provisional_*, figures/qc_outliers, results/figures/, results/archive_v2_strict/); "
             "graphs at the declared 1.5 A cut (multi-cut graphs only where the claim is the 1.5 A value).")
    L.append("- not_cut originals are copied as they are (one file, several graphs): "
             + ", ".join(r["file"] for r in rows if (manifest.loc[manifest["source_path(s)"].str.split(';').str[0]
                                                               == r["source"], "method"] == "not_cut").any()) + ".")
    L.append("")
    L.append("## Files")
    L.append("")
    out_files = [VIS / "INDEX.md", RES / "selection.csv", RES / "check.csv", RES / "SUMMARY.md", RES / "env.txt",
                 ROOT / "scripts/sv/VIS2_build.py", ROOT / "scripts/sv/VIS2_check.py",
                 ROOT / "scripts/sv/VIS2_summary.py"]
    for f in out_files:
        if f.exists() or f.name in ("SUMMARY.md", "env.txt"):
            L.append(f"- {f.relative_to(ROOT)}" + (f" ({size(f)})" if f.exists() else ""))
    tot = sum(p.stat().st_size for p in pngs)
    L.append(f"- visuals/*.png: {len(pngs)} files ({tot / 1024 / 1024:.1f} MB)")
    if "## Suggested, not done (VIS2" in (ROOT / "docs/EXPLORATION.md").read_text():
        L.append("- docs/EXPLORATION.md: appended 'Suggested, not done (VIS2, 2026-11-08)' (rule 1.1)")
    L.append("")
    L.append("Status proposed: " + ("DONE" if n_fail == 0 else "BLOCKED (check failures)"))

    (RES / "SUMMARY.md").write_text("\n".join(L) + "\n")
    env = [f"# VIS2 environment, {time.strftime('%Y-%m-%dT%H:%M%z')}",
           f"python {platform.python_version()} at {sys.executable}",
           f"platform {platform.platform()}",
           f"pandas {version('pandas')}",
           f"Pillow {version('Pillow')}",
           "copy: shutil.copy2 (files copied, never moved); hashes: hashlib.sha1",
           "network: none",
           "nothing installed, upgraded or downgraded by VIS2"]
    (RES / "env.txt").write_text("\n".join(env) + "\n")
    # second pass so SUMMARY.md lists its own final size
    text = (RES / "SUMMARY.md").read_text()
    text = re.sub(r"- results/sv/VIS2/SUMMARY.md( \([^)]*\))?", f"- results/sv/VIS2/SUMMARY.md ({size(RES / 'SUMMARY.md')})", text)
    text = re.sub(r"- results/sv/VIS2/env.txt( \([^)]*\))?", f"- results/sv/VIS2/env.txt ({size(RES / 'env.txt')})", text)
    (RES / "SUMMARY.md").write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
