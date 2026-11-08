"""SV00: one-off test of the Appendix B RCSB Data API query on 5TM1 (section 3: test once on 5TM1 first).
Output: data/sv/api_cache/test_5TM1.json (raw response).
Run: .venv/bin/python -I scripts/sv/SV00_api_test.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from SV00_metadata import API_QUERY, post_graphql  # noqa: E402

res = post_graphql(["5TM1"])
(ROOT / "data/sv/api_cache/test_5TM1.json").write_text(json.dumps(res, indent=1))
print("errors:", res.get("errors"))
e = res["data"]["entries"][0]
print({k: (v if not isinstance(v, list) else v[:2]) for k, v in e.items() if k not in ("polymer_entities", "nonpolymer_entities")})
print(json.dumps(e["polymer_entities"][0], indent=0)[:1500])
print(json.dumps(e["nonpolymer_entities"][:2], indent=0)[:800])
