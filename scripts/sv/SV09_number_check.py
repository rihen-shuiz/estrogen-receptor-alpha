"""SV09 number check: every number in the new docs text must match the file it comes from.

The SV09 text is written as a tagged source (docs/sv/SV09_section.src.md, docs/sv/SV09_limitations.src.md).
Each number carries a tag

    ⟦display ‖ alias ‖ locator ‖ filters ‖ transform⟧

- alias:     a short name for a saved file (ALIASES below).
- locator:   CSV -> column name, or #count (rows matching the filters), #nrows, #ncols;
             JSON -> path of keys / list indices separated by '/';
             MD / TXT -> regular expression whose group 1 is the number (filters = occurrence index, default 0).
- filters:   CSV only, 'col=value&col~substring' (numeric values compared as numbers).
- transform: comma-separated, applied in order: x100, abs, neg, len, min, max, first.

scripts/sv/SV09_docs.py renders the sources (tags -> display) and appends them to docs/RESULTS_SUMMARY.md and
docs/limitations.md. This script then
1. extracts every number from the new section of each doc (code spans, identifiers such as H12, cl7, SV01b,
   Y537S, 5TM1, H10/11 and A^3 exponents excluded),
2. aligns them with the numbers of the rendered source, so each doc number is tied to its tag (a number with no
   tag is a failure),
3. reads the tagged value from its file and compares it with the displayed value at the displayed precision,
4. writes results/sv/SV09/number_check.csv (doc, number, context, source file, column/key, value found, match).

Run: .venv/bin/python -B scripts/sv/SV09_number_check.py   (exit code 1 if any number fails)
"""
from __future__ import annotations

import csv
import difflib
import json
import re
import sys
from functools import lru_cache
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

SOURCES = {
    "RESULTS_SUMMARY": (ROOT / "docs/sv/SV09_section.src.md", ROOT / "docs/RESULTS_SUMMARY.md",
                        "## Supervisor requests 2026-10-28"),
    "limitations": (ROOT / "docs/sv/SV09_limitations.src.md", ROOT / "docs/limitations.md",
                    "# Limitations added by RUN_SV1028"),
}
OUT = ROOT / "results/sv/SV09/number_check.csv"

ALIASES = {
    "FB": "docs/supervisor_feedback.md",
    "S00g": "results/sv/SV00/gate.json",
    "S00c": "results/sv/SV00/counts.csv",
    "S00m": "data/sv/manifest.csv",
    "S00miss": "results/sv/SV00/missing_by_column.csv",
    "S00t": "results/sv/SV00/chain_table.csv",
    "S00sum": "results/sv/SV00/SUMMARY.md",
    "S01a_h": "results/sv/SV01a/h12_status_counts.csv",
    "S01a_pe": "results/sv/SV01a/partial_ends.csv",
    "S01a_cc": "results/sv/SV01a/concordance_summary.csv",
    "S01a_vt": "results/sv/SV01a/verdict_table.csv",
    "S01a_mb": "results/sv/SV01a/model_blocks.csv",
    "S01b_vt": "results/sv/SV01b/verdict_table.csv",
    "S01b_or": "results/sv/SV01b/odds_ratios.csv",
    "S01b_lg": "results/sv/SV01b/lattice_contacts_by_group.csv",
    "S01b_mb": "results/sv/SV01b/model_blocks.csv",
    "S01b_sum": "results/sv/SV01b/SUMMARY.md",
    "S01c_n": "results/sv/SV01c/numbers.json",
    "S01c_sb": "results/sv/SV01c/side_by_side.csv",
    "S01c_st": "results/sv/SV01c/status_transitions.csv",
    "S01c_cl": "results/sv/SV01c/clash_summary.csv",
    "S02_hb": "results/sv/SV02/helix_boundaries.csv",
    "S02_nc": "results/sv/SV02/numbers_consensus.json",
    "S02_nr": "results/sv/SV02/numbers_rama.json",
    "S02_rs": "results/sv/SV02/rama_summary.csv",
    "S02_hw": "results/sv/SV02/h12_window.md",
    "S03_nt": "results/sv/SV03/numbers_tests.json",
    "S03_w": "results/sv/SV03/cl7_vs_cl4_wide.csv",
    "S03_ax": "results/sv/SV03/axis_method_check.csv",
    "S03_sum": "results/sv/SV03/SUMMARY.md",
    "S04_ari": "results/sv/SV04/ari.csv",
    "S04_ds": "results/sv/SV04/distance_summary.csv",
    "S04_nc": "results/sv/SV04/numbers_cluster.json",
    "S04_ve": "results/sv/SV04/validation_error.csv",
    "S04_sum": "results/sv/SV04/SUMMARY.md",
    "S05_n": "results/sv/SV05/numbers_core_v4.json",
    "S05_ari": "results/sv/SV05/ari.csv",
    "S05_sum": "results/sv/SV05/SUMMARY.md",
    "S06_n": "results/sv/SV06/numbers.json",
    "S06_h": "results/sv/SV06/hypergeom_tests.csv",
    "S07_n": "results/sv/SV07/numbers.json",
    "S07_v": "results/sv/SV07/verdict.csv",
    "S07_s": "results/sv/SV07/spikes.csv",
    "S07_ct": "results/sv/SV07/covariate_tests.csv",
    "S07_sum": "results/sv/SV07/SUMMARY.md",
    "S08_sum": "results/sv/SV08/SUMMARY.md",
    "S08_lit": "docs/sv/literature.md",
    "S08b_sum": "results/sv/SV08b/SUMMARY.md",
    "S08b_lit": "docs/sv/literature_class.md",
    "S10_n": "results/sv/SV10/numbers.json",
    "S11_n": "results/sv/SV11/numbers_partA.json",
    "S11_st": "results/sv/SV11/site_tests.csv",
    "S11_ster": "results/sv/SV11/steroid_tests.csv",
    "S11_sum": "results/sv/SV11/SUMMARY.md",
    "S12_n": "results/sv/SV12/numbers.json",
    "S12_cov": "results/sv/SV12/class_coverage.csv",
    "S12_sum": "results/sv/SV12/SUMMARY.md",
    "S13_n": "results/sv/SV13/numbers.json",
    "S13_g": "results/sv/SV13/gate.csv",
    "S13_r": "results/sv/SV13/representatives.csv",
    "S13_sum": "results/sv/SV13/SUMMARY.md",
    "S14_n": "results/sv/SV14/numbers.json",
    "S14_c": "results/sv/SV14/correlations.csv",
    "S14_t": "results/sv/SV14/tests.csv",
    "S14_ref": "results/sv/SV14/references.json",
    "S14_fg": "results/sv/SV14/filling_by_group.csv",
    "S14_sum": "results/sv/SV14/SUMMARY.md",
    "S15_n": "results/sv/SV15/numbers.json",
    "S15_t": "results/sv/SV15/tests.csv",
    "S15_t1": "results/sv/SV15/tests_first_run_type_units.csv",
    "S15_a": "results/sv/SV15/affinity_per_ligand.csv",
    "S15_v": "results/sv/SV15/verdict.md",
}

TAG = re.compile(r"⟦([^⟧]*)⟧")
NUM = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
# masked before tokenising: code spans, links, helix pairs such as H4/5 and H10/11
MASK = re.compile(r"`[^`\n]*`|https?://\S+|\b[A-Za-z]+\d+(?:/\d+)+\b")
IDENT_BEFORE = set("_|^.")
IDENT_AFTER = set("_|^")


# --------------------------------------------------------------------------------------------- tags / rendering
def parse_tag(body: str) -> dict:
    parts = [p.strip() for p in body.split("‖")]
    parts += [""] * (5 - len(parts))
    display, alias, locator, filters, transform = parts[:5]
    return {"display": display, "alias": alias, "locator": locator, "filters": filters, "transform": transform}


def render(src_text: str) -> tuple[str, list[tuple[int, int, dict]]]:
    """Replace every tag by its display value; return the text and (start, end, tag) spans in it."""
    out, spans, pos = [], [], 0
    cursor = 0
    for m in TAG.finditer(src_text):
        out.append(src_text[pos:m.start()])
        cursor += m.start() - pos
        tag = parse_tag(m.group(1))
        out.append(tag["display"])
        spans.append((cursor, cursor + len(tag["display"]), tag))
        cursor += len(tag["display"])
        pos = m.end()
    out.append(src_text[pos:])
    return "".join(out), spans


# --------------------------------------------------------------------------------------------- tokeniser
def tokenize(text: str) -> list[tuple[int, int, str]]:
    masked = MASK.sub(lambda m: " " * len(m.group()), text)
    toks = []
    for m in NUM.finditer(masked):
        s, e, tok = m.start(), m.end(), m.group()
        prev = masked[s - 1] if s > 0 else " "
        if tok[0] in "+-" and (prev.isalnum() or prev == ")"):
            if prev.isalpha():
                continue  # hyphenated identifier such as AF-2
            # a dash between two numbers (538-547, 0.5-5.4): not a sign
            tok, s, prev = tok[1:], s + 1, "-"
        if prev.isalnum() or prev in IDENT_BEFORE:
            continue
        nxt = masked[e] if e < len(masked) else " "
        if nxt.isalpha() or nxt in IDENT_AFTER:
            continue
        toks.append((s, e, tok))
    return toks


# --------------------------------------------------------------------------------------------- value lookup
@lru_cache(maxsize=None)
def _csv(path: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / path, dtype=str, keep_default_na=False)


@lru_cache(maxsize=None)
def _json(path: str):
    return json.loads((ROOT / path).read_text(encoding="utf8"))


@lru_cache(maxsize=None)
def _text(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf8")


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _eq(cell: str, want: str) -> bool:
    a, b = _num(cell), _num(want)
    if a is not None and b is not None:
        return abs(a - b) < 1e-12
    return cell == want


def lookup(tag: dict):
    alias = tag["alias"]
    if alias not in ALIASES:
        raise KeyError(f"unknown alias {alias!r}")
    path = ALIASES[alias]
    loc, flt, tr = tag["locator"], tag["filters"], [t for t in tag["transform"].split(",") if t]
    suffix = Path(path).suffix
    if suffix == ".csv":
        df = _csv(path)
        sub = df
        if flt:
            for f in flt.split("&"):
                if "~" in f and ("=" not in f or f.index("~") < f.index("=")):
                    col, pat = f.split("~", 1)
                    sub = sub[sub[col].str.contains(pat, regex=False)]
                else:
                    col, val = f.split("=", 1)
                    sub = sub[sub[col].map(lambda c: _eq(c, val))]
        if loc == "#count":
            value = len(sub)
        elif loc == "#nrows":
            value = len(sub)
        elif loc == "#ncols":
            value = df.shape[1]
        else:
            if len(sub) == 0:
                raise LookupError("no row matches the filters")
            if len(sub) > 1 and "first" not in tr:
                raise LookupError(f"{len(sub)} rows match the filters")
            value = sub.iloc[0][loc]
    elif suffix == ".json":
        cur = _json(path)
        for seg in loc.split("/"):
            cur = cur[int(seg)] if isinstance(cur, list) else cur[seg]
        value = cur
    else:
        hits = list(re.finditer(loc, _text(path)))
        idx = int(flt) if flt else 0
        if len(hits) <= idx:
            raise LookupError("regex not found")
        value = hits[idx].group(1)
    for t in tr:
        if t == "first":
            continue
        if t == "len":
            value = len(value)
        elif t == "min":
            value = min(value)
        elif t == "max":
            value = max(value)
        elif t == "x100":
            value = float(value) * 100
        elif t == "abs":
            value = abs(float(value))
        elif t == "neg":
            value = -float(value)
        else:
            raise ValueError(f"unknown transform {t}")
    return path, value


def matches(display: str, value) -> bool:
    d = display.lstrip("+")
    dv = _num(d)
    if dv is None:  # non-numeric display such as a date: literal comparison
        return str(value) == display
    v = _num(value)
    if v is None:
        return False
    if display.startswith("+") and v < 0:
        return False
    if "e" in d.lower():
        mant = d.lower().split("e")[0].lstrip("-")
        nd = len(mant.split(".")[1]) if "." in mant else 0
        return float(f"{v:.{nd}e}") == dv
    nd = len(d.split(".")[1]) if "." in d else 0
    return abs(v - dv) <= 0.5 * 10 ** (-nd) + 1e-9 * max(1.0, abs(v))


def fmt_value(value) -> str:
    if isinstance(value, float):
        return f"{value:.10g}"
    return str(value)


# --------------------------------------------------------------------------------------------- check
def section_of(doc_text: str, marker: str) -> str | None:
    i = doc_text.find(marker)
    return None if i < 0 else doc_text[i:]


def check_doc(name: str, src: Path, doc: Path, marker: str) -> list[dict]:
    rendered, spans = render(src.read_text(encoding="utf8"))
    section = section_of(doc.read_text(encoding="utf8"), marker)
    rows = []
    if section is None:
        return [{"doc": name, "number": "", "context": f"marker {marker!r} not found in {doc.name}",
                 "source_file": "", "column_key": "", "value_found": "", "match": False}]
    exp = tokenize(rendered)
    exp_tags = []
    for s, e, tok in exp:
        owner = next((t for (a, b, t) in spans if a <= s and e <= b), None)
        exp_tags.append(owner)
    for a, b, t in spans:  # a tag whose display the tokeniser does not count would escape the check
        if not any(a <= s and e <= b for s, e, _ in exp):
            rows.append({"doc": name, "number": t["display"], "context": "tag display not counted as a number",
                         "source_file": ALIASES.get(t["alias"], t["alias"]), "column_key": t["locator"],
                         "value_found": "", "match": False})
    got = tokenize(section)
    if section.rstrip() != rendered.rstrip():
        rows.append({"doc": name, "number": "", "context": "doc section differs from the rendered source",
                     "source_file": str(src.relative_to(ROOT)), "column_key": "", "value_found": "",
                     "match": False})
    sm = difflib.SequenceMatcher(a=[t for _, _, t in got], b=[t for _, _, t in exp], autojunk=False)
    pairs = {}
    for blk in sm.get_matching_blocks():
        for k in range(blk.size):
            pairs[blk.a + k] = blk.b + k
    for i, (s, e, tok) in enumerate(got):
        ctx = section[max(0, s - 60):min(len(section), e + 60)].replace("\n", " ")
        row = {"doc": name, "number": tok, "context": ctx.strip()}
        j = pairs.get(i)
        tag = exp_tags[j] if j is not None else None
        if tag is None:
            row.update(source_file="UNTAGGED", column_key="", value_found="", match=False)
            rows.append(row)
            continue
        key = tag["locator"]
        if tag["filters"]:
            key += f" [{tag['filters']}]"
        if tag["transform"]:
            key += f" ({tag['transform']})"
        try:
            path, value = lookup(tag)
            ok = matches(tag["display"], value) if NUM.fullmatch(tag["display"].lstrip("+")) else \
                str(value) == tag["display"]
            if tok.lstrip("+") not in tag["display"].lstrip("+"):
                ok = False
            row.update(source_file=path, column_key=key, value_found=fmt_value(value), match=bool(ok))
        except Exception as exc:  # noqa: BLE001 - report the failure in the table
            row.update(source_file=ALIASES.get(tag["alias"], tag["alias"]), column_key=key,
                       value_found=f"ERROR: {exc}", match=False)
        rows.append(row)
    matched_exp = set(pairs.values())
    for j, (s, e, tok) in enumerate(exp):
        if j not in matched_exp:
            rows.append({"doc": name, "number": tok, "context": "in the source but not in the doc",
                         "source_file": "", "column_key": "", "value_found": "", "match": False})
    return rows


def main() -> int:
    rows = []
    for name, (src, doc, marker) in SOURCES.items():
        rows.extend(check_doc(name, src, doc, marker))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["doc", "number", "context", "source_file", "column_key", "value_found", "match"]
    with OUT.open("w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    n = len(rows)
    bad = [r for r in rows if not r["match"]]
    by_doc = {}
    for r in rows:
        by_doc.setdefault(r["doc"], [0, 0])
        by_doc[r["doc"]][0] += 1
        by_doc[r["doc"]][1] += int(bool(r["match"]))
    print(f"numbers checked: {n}; matched: {n - len(bad)}; failed: {len(bad)}")
    for d, (a, b) in by_doc.items():
        print(f"  {d}: {b}/{a} matched")
    print(f"sources used: {len({r['source_file'] for r in rows if r['source_file']})}")
    for r in bad[:80]:
        print("FAIL", r["doc"], repr(r["number"]), "|", r["value_found"], "|", r["column_key"], "|",
              r["context"][:110])
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
