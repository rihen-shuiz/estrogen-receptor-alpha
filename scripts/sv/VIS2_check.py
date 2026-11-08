"""VIS2 check (RUN_SV1028, step 5): verifies visuals/INDEX.md against the documents and the files.

Checks (one row each in results/sv/VIS2/check.csv):
  claim      every quoted claim of every INDEX row exists in the cited document, inside the cited location
             (RESULTS_SUMMARY: 'The data', '§N', 'Storyline' = that '## ' section; 'SVxx' = that step paragraph;
             other documents: whole file). Whitespace is collapsed on both sides (the document wraps lines) and \\|
             in the table reads as |.
  listed     every file in visuals/ (INDEX.md itself excluded) is listed in INDEX, and every listed file exists
  opens      every listed file opens as an image (PIL verify + full decode)
  sha1       no two files in visuals/ share a SHA-1
Extra integrity checks: the source path of every row is a VIS1 pool item (figures/panels/manifest.csv) and the copy
has the same SHA-1 as its source; numbering is 01..N without gaps, a 'b' row follows its graph; every main-figure
file of the SV09 table in RESULTS_SUMMARY is represented by an original or one of its panels.

Run: .venv/bin/python -B scripts/sv/VIS2_check.py   (exit code 1 if any check fails)
"""
from __future__ import annotations

import csv
import hashlib
import re
import sys
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
VIS = ROOT / "visuals"
INDEX = VIS / "INDEX.md"
MANIFEST = ROOT / "figures/panels/manifest.csv"
RS = ROOT / "docs/RESULTS_SUMMARY.md"
OUT = ROOT / "results/sv/VIS2/check.csv"

CLAIM_RE = re.compile(r'"((?:[^"\\]|\\.)*)" \(([^,()]+), ([^()]+)\)')
SPLIT_RE = re.compile(r"(?<!\\)\|")


def norm(text: str) -> str:
    return " ".join(text.split())


def sha1(path: Path) -> str:
    return hashlib.sha1(path.read_bytes()).hexdigest()


def location_text(doc: Path, loc: str) -> str | None:
    lines = doc.read_text().splitlines()
    if doc.resolve() != RS.resolve():
        return "\n".join(lines)
    if loc == "The data":
        start = lambda s: s.strip() == "## The data"
        stop = lambda s: s.startswith("## ")
    elif loc == "Storyline":
        start = lambda s: s.startswith("## Storyline")
        stop = lambda s: s.startswith("## ")
    elif re.fullmatch(r"§\d+", loc):
        num = loc[1:]
        start = lambda s: s.startswith(f"## {num}.")
        stop = lambda s: s.startswith("## ")
    elif re.fullmatch(r"SV\d+[a-z]?", loc):
        start = lambda s: s.startswith(f"**{loc} - ")
        stop = lambda s: s.startswith("**SV") or s.startswith("#")
    else:
        return None
    out, inside = [], False
    for line in lines:
        if inside and stop(line):
            break
        if not inside and start(line):
            inside = True
        if inside:
            out.append(line)
    return "\n".join(out) if out else None


def parse_index():
    rows = []
    for line in INDEX.read_text().splitlines():
        if not line.startswith("| ") or line.startswith("| NN ") or line.startswith("|---"):
            continue
        cells = [c.strip() for c in SPLIT_RE.split(line.strip())[1:-1]]
        if len(cells) != 5:
            rows.append({"bad": line})
            continue
        nn, fcell, what, claim, src = cells
        rows.append({"nn": nn, "file": fcell.strip("`"), "what": what, "claim": claim,
                     "source": re.findall(r"`([^`]+)`", src)[0] if "`" in src else ""})
    return rows


def main() -> int:
    results = []

    def add(check, item, ok, detail=""):
        results.append({"check": check, "item": item, "result": "pass" if ok else "FAIL", "detail": detail})

    rows = parse_index()
    for r in rows:
        if "bad" in r:
            add("parse", r["bad"][:60], False, "row does not have 5 cells")
    rows = [r for r in rows if "bad" not in r]

    # claims
    for r in rows:
        claims = CLAIM_RE.findall(r["claim"])
        if not claims:
            add("claim", r["file"], False, "no quoted claim found")
        for quote, doc, loc in claims:
            q = norm(quote.replace("\\|", "|"))
            dpath = ROOT / doc.strip()
            if not dpath.exists():
                add("claim", r["file"], False, f"document missing: {doc}")
                continue
            text = location_text(dpath, loc.strip())
            if text is None:
                add("claim", r["file"], False, f"location not found: {doc}, {loc}")
                continue
            ok = q in norm(text)
            add("claim", r["file"], ok, f"{doc}, {loc}: {q[:70]}" + ("" if ok else " -- NOT FOUND"))

    # listed / opens
    listed = [r["file"] for r in rows]
    on_disk = sorted(p.name for p in VIS.iterdir() if p.is_file() and p.name != "INDEX.md")
    for name in on_disk:
        add("listed", name, name in listed, "" if name in listed else "file not in INDEX")
    for name in listed:
        path = VIS / name
        if not path.exists():
            add("listed", name, False, "listed in INDEX but missing")
            continue
        try:
            with Image.open(path) as im:
                im.verify()
            with Image.open(path) as im:
                im.load()
                size = im.size
            add("opens", name, True, f"{size[0]}x{size[1]}")
        except Exception as exc:  # noqa: BLE001 - any decode error is a failure
            add("opens", name, False, repr(exc))

    # sha1 unique
    hashes = {}
    for name in on_disk:
        hashes.setdefault(sha1(VIS / name), []).append(name)
    dups = {h: n for h, n in hashes.items() if len(n) > 1}
    add("sha1", f"{len(on_disk)} files", not dups, "; ".join(",".join(n) for n in dups.values()) or
        f"{len(hashes)} distinct SHA-1")

    # source integrity
    pool, pool_src = set(), {}
    with MANIFEST.open(newline="") as fh:
        for m in csv.DictReader(fh):
            key = m["source_path(s)"].split(";")[0] if m["method"] in ("single", "not_cut") else m["panel_path"]
            pool.add(key)
            pool_src[key] = m["source_path(s)"].split(";")[0]
    for r in rows:
        src = ROOT / r["source"]
        ok = r["source"] in pool and src.exists() and (VIS / r["file"]).exists() and \
            sha1(src) == sha1(VIS / r["file"])
        add("source", r["file"], ok, r["source"] if ok else f"not a pool item or SHA-1 differs: {r['source']}")

    # numbering
    expect, prev_main = 1, None
    for r in rows:
        nn = r["nn"]
        if nn.endswith("b"):
            ok = prev_main is not None and nn[:-1] == prev_main and r["file"].startswith(nn + "_")
            add("numbering", r["file"], ok, "legend follows its graph" if ok else "b row out of place")
        else:
            ok = nn == f"{expect:02d}" and r["file"].startswith(nn + "_")
            add("numbering", r["file"], ok, "" if ok else f"expected {expect:02d}")
            prev_main = nn
            expect += 1

    # SV09 main figures represented
    rs = RS.read_text()
    m = re.search(r"### Main figures for the next meeting\n(.*?)\n### ", rs, re.S)
    mains = sorted(set(re.findall(r"`(figures/sv/[^`]+\.png)`", m.group(1)))) if m else []
    add("sv09_main", "table found", bool(mains), f"{len(mains)} files")
    used = {r["source"] for r in rows} | {pool_src.get(r["source"], "") for r in rows}
    for f in mains:
        add("sv09_main", f, f in used, "" if f in used else "main figure not in visuals/")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["check", "item", "result", "detail"])
        w.writeheader()
        w.writerows(results)
    n_fail = sum(r["result"] == "FAIL" for r in results)
    by = {}
    for r in results:
        by.setdefault(r["check"], [0, 0])[r["result"] == "FAIL"] += 1
    for c, (ok, bad) in by.items():
        print(f"{c}: {ok} pass, {bad} fail")
    for r in results:
        if r["result"] == "FAIL":
            print("FAIL", r["check"], r["item"], r["detail"])
    print("ALL PASS" if n_fail == 0 else f"{n_fail} FAILED")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
