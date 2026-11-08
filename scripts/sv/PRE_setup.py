"""PRE (RUN_SV1028): read-only snapshot, RUN_STATUS skeleton, per-step briefs. Stdlib only."""
import csv, hashlib, os, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "docs/sv/RUN_SV1028.md"
text = RUN.read_text(encoding="utf-8")

def section(start_pat, end_pat):
    m = re.search(start_pat, text, re.M)
    if not m:
        sys.exit(f"section not found: {start_pat}")
    rest = text[m.start():]
    nl = rest.index("\n") + 1  # search for the end after the heading line
    e = re.search(end_pat, rest[nl:], re.M)
    return rest[: e.start() + nl if e else len(rest)].rstrip().rstrip("-").rstrip()

def snapshot():
    rows = []
    for top in ["results", "figures", "data/clean", "esr_analysis"]:
        base = ROOT / top
        if not base.exists():
            continue
        for p in sorted(base.rglob("*")):
            if p.is_file():
                h = hashlib.sha1(p.read_bytes()).hexdigest()
                rows.append((p.relative_to(ROOT).as_posix(), h, p.stat().st_size))
    with open(ROOT / "docs/sv/readonly_snapshot.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["path", "sha1", "bytes"]); w.writerows(rows)
    return len(rows)

STEPS = ["PRE", "SV00", "SV01a", "SV06", "SV10", "SV11", "SV02", "SV01b", "SV03", "SV04", "SV05",
         "SV07", "SV01c", "SV08", "PRE addition", "SV12", "SV14", "SV15", "SV13", "SV08b", "SV09",
         "VIS1", "VIS2", "FINAL"]

def status():
    p = ROOT / "docs/sv/RUN_STATUS.md"
    if p.exists():
        return "exists"
    lines = ["# RUN_STATUS — RUN_SV1028", "", "| step | status | result | main outputs |", "|---|---|---|---|"]
    lines += [f"| {s} | TODO | | |" for s in STEPS]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return "created"

def briefs():
    rules = section(r"^## 1\. Rules for every step", r"^## 2\.")
    order = section(r"^## 2\. Order and dependencies", r"^## 3\.")
    data = section(r"^## 3\. Data access", r"^## 4\.")
    addition = section(r"^### Addition 2026-11-04", r"^### PRE addition")
    integ = section(r"^### Integration of the addition", r"^### SV09")
    appB = section(r"^## Appendix B", r"^## Appendix C")
    step_heads = {
        "SV00": r"^### SV00 ", "SV01": r"^### SV01 ", "SV02": r"^### SV02 ", "SV03": r"^### SV03 ",
        "SV04": r"^### SV04 ", "SV05": r"^### SV05 ", "SV06": r"^### SV06 ", "SV07": r"^### SV07 ",
        "SV08": r"^### SV08 ", "SV10": r"^### SV10 ", "SV11": r"^### SV11 ", "SV12": r"^### SV12 ",
        "SV14": r"^### SV14 ", "SV15": r"^### SV15 ", "SV13": r"^### SV13 ", "SV08b": r"^### SV08b ",
        "SV09": r"^### SV09 ", "VIS1": r"^### VIS1 ", "VIS2": r"^### VIS2 ",
    }
    steps = {k: section(v, r"^(### |## )") for k, v in step_heads.items()}
    need_data = {"SV00", "SV01a", "SV01b", "SV11", "SV12", "SV15"}
    ids = ["SV00", "SV01a", "SV01b", "SV01c", "SV02", "SV03", "SV04", "SV05", "SV06", "SV07", "SV08",
           "SV10", "SV11", "SV12", "SV14", "SV15", "SV13", "SV08b", "SV09", "VIS1", "VIS2"]
    for i in ids:
        key = i[:4] if i.startswith("SV01") else i
        parts = [f"# Brief for step {i} (RUN_SV1028)", "",
                 "Extracted verbatim from docs/sv/RUN_SV1028.md: §1, §2" + (", §3" if i in need_data else "")
                 + " and the step text.", "", rules, "", order]
        if i in need_data:
            parts += ["", data]
        if i == "SV00":
            parts += ["", appB]
        if i in {"SV12", "SV14", "SV15", "SV13", "SV08b"}:
            parts += ["", addition]
        if i in {"SV09", "VIS2"}:
            parts += ["", integ]
        parts += ["", steps[key]]
        if i.startswith("SV01"):
            parts += ["", f"**This subagent runs only part {i} of SV01** (outputs in results/sv/{i}/)."]
        (ROOT / f"docs/sv/briefs/{i}.md").write_text("\n".join(parts) + "\n", encoding="utf-8")
    return len(ids)

if __name__ == "__main__":
    print("snapshot files:", snapshot())
    print("RUN_STATUS:", status())
    print("briefs:", briefs())
