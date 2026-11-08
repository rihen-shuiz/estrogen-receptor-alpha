"""SV12 env.txt: interpreter, package versions, network endpoints and the ChEMBL release (status resource, cached)."""
from __future__ import annotations

import json
import platform
import sys
import time
from importlib.metadata import version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv/common"))
import cached_http as ch  # noqa: E402


def main():
    st, path = ch.fetch_json(ch.chembl_url("status"), ROOT / "data/sv/chembl_cache", "status")
    lines = [f"# SV12 environment, {time.strftime('%Y-%m-%dT%H:%M%z')}",
             f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}"]
    for p in ("pandas", "numpy", "requests", "gemmi"):
        lines.append(f"{p} {version(p)}")
    lines += [f"ChEMBL release: {st.get('chembl_db_version')} (release date {st.get('chembl_release_date')}; "
              f"status resource cached {path.relative_to(ROOT)})",
              "network: data.rcsb.org REST chemcomp (306 CCD records, cache data/sv/rcsb_chemcomp/), "
              "data.rcsb.org GraphQL rcsb_primary_citation (5TM1 test + 1 batch, cache data/sv/api_cache/sv12_*), "
              "www.ebi.ac.uk/chembl/api/data (cache data/sv/chembl_cache/, index.jsonl lists every URL)",
              "client: scripts/sv/common/cached_http.py (<= 4 requests/s, retries with exponential backoff on "
              "network errors, 429 and 5xx; every page via page_meta.next)",
              "nothing installed, upgraded or downgraded by SV12",
              "teammate file data/sv/action_class_teammate.csv: "
              + ("present" if (ROOT / "data/sv/action_class_teammate.csv").exists() else "absent at run time")]
    (ROOT / "results/sv/SV12/env.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(json.dumps(st))


if __name__ == "__main__":
    main()
