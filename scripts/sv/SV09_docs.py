"""SV09: render the tagged SV09 texts and append them to the two docs the step edits.

- docs/sv/SV09_section.src.md     -> docs/sv/SV09_section.md     -> end of docs/RESULTS_SUMMARY.md
- docs/sv/SV09_limitations.src.md -> docs/sv/SV09_limitations.md -> end of docs/limitations.md

Tags (see scripts/sv/SV09_number_check.py) are replaced by their display value. Earlier text of both docs is kept
byte for byte: the new text is appended after it; if the SV09 heading is already present (a re-run), only the text
from that heading to the end of the file is replaced. docs/RESULTS_SUMMARY.md is rendered from
docs/templates/RESULTS_SUMMARY.md by scripts/collect_numbers_ext.py; a re-render drops this section, which can then
be re-appended by running this script again.

Run: .venv/bin/python -B scripts/sv/SV09_docs.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from SV09_number_check import ROOT, SOURCES, render  # noqa: E402

PLAIN = {
    "RESULTS_SUMMARY": ROOT / "docs/sv/SV09_section.md",
    "limitations": ROOT / "docs/sv/SV09_limitations.md",
}


def main() -> None:
    for name, (src, doc, marker) in SOURCES.items():
        text, spans = render(src.read_text(encoding="utf8"))
        text = text.rstrip("\n") + "\n"
        if not text.startswith(marker):
            raise SystemExit(f"{src} must start with {marker!r}")
        PLAIN[name].write_text(text, encoding="utf8")
        old = doc.read_text(encoding="utf8")
        i = old.find(marker)
        head = old if i < 0 else old[:i]
        head = head.rstrip("\n") + "\n\n"
        doc.write_text(head + text, encoding="utf8")
        print(f"{doc.relative_to(ROOT)}: kept {len(head)} chars, appended {len(text)} chars "
              f"({len(spans)} tagged numbers); copy {PLAIN[name].relative_to(ROOT)}")


if __name__ == "__main__":
    main()
