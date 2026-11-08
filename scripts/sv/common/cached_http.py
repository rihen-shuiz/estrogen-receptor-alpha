"""Polite cached JSON GET/POST for RCSB and ChEMBL (RUN_SV1028 §3).

- At most 4 requests per second (min interval 0.25 s, below the §3 limit of 5/s), retries with exponential backoff
  on network errors, HTTP 429 and 5xx.
- Every response body is cached as raw JSON text in the given cache directory, under a file named
  <prefix>_<sha1 of the full request>.json; an index.jsonl next to it records file -> URL (and POST body).
- A cached response is reused instead of a new request (re-runs are offline once the cache is complete).
- chembl_get_all() pages through a ChEMBL list resource with page_meta.next and caches every page.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from urllib.parse import urlencode

import requests

MIN_INTERVAL_S = 0.25
MAX_TRIES = 6
TIMEOUT_S = 120
CHEMBL_HOST = "https://www.ebi.ac.uk"
CHEMBL_BASE = CHEMBL_HOST + "/chembl/api/data/"

_last = [0.0]
_session = requests.Session()
_session.headers.update({"User-Agent": "esr1-sv1028/1.0 (academic project; polite cached client)",
                         "Accept": "application/json"})
STATS = {"network": 0, "cache": 0, "retries": 0}


def _wait():
    dt = time.monotonic() - _last[0]
    if dt < MIN_INTERVAL_S:
        time.sleep(MIN_INTERVAL_S - dt)
    _last[0] = time.monotonic()


def _key(url: str, body: str | None) -> str:
    return hashlib.sha1((url + "\n" + (body or "")).encode()).hexdigest()


def fetch_json(url: str, cache_dir: Path, prefix: str, body: dict | None = None, allow_404: bool = False):
    """GET (or POST with a JSON body) url; return (parsed JSON or None for an allowed 404, cache path)."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    body_s = json.dumps(body, sort_keys=True) if body is not None else None
    path = cache_dir / f"{prefix}_{_key(url, body_s)}.json"
    if path.exists():
        STATS["cache"] += 1
        txt = path.read_text()
        return (None if txt == "__404__" else json.loads(txt)), path
    err = None
    for i in range(MAX_TRIES):
        _wait()
        try:
            if body is None:
                r = _session.get(url, timeout=TIMEOUT_S)
            else:
                r = _session.post(url, json=body, timeout=TIMEOUT_S)
            STATS["network"] += 1
            if r.status_code == 404 and allow_404:
                path.write_text("__404__")
                _index(cache_dir, path, url, body_s, 404)
                return None, path
            if r.status_code == 429 or r.status_code >= 500:
                err = f"HTTP {r.status_code}"
            else:
                r.raise_for_status()
                data = r.json()
                path.write_text(r.text)
                _index(cache_dir, path, url, body_s, r.status_code)
                return data, path
        except (requests.ConnectionError, requests.Timeout, json.JSONDecodeError) as e:
            err = repr(e)
        STATS["retries"] += 1
        time.sleep(min(60, 2 ** i))
    raise RuntimeError(f"failed after {MAX_TRIES} tries: {url} ({err})")


def _index(cache_dir: Path, path: Path, url: str, body_s: str | None, status: int):
    with open(cache_dir / "index.jsonl", "a") as fh:
        fh.write(json.dumps({"file": path.name, "url": url, "body": body_s, "status": status,
                             "fetched": time.strftime("%Y-%m-%dT%H:%M:%S%z")}) + "\n")


def chembl_url(resource: str, params: dict | None = None) -> str:
    q = urlencode(params or {}, safe=",")
    return f"{CHEMBL_BASE}{resource}.json" + (f"?{q}" if q else "")


def chembl_get_all(resource: str, params: dict, cache_dir: Path, prefix: str | None = None) -> list[dict]:
    """All records of a ChEMBL list resource, following page_meta.next; every page cached."""
    prefix = prefix or resource
    p = dict(params)
    p.setdefault("limit", 1000)
    url = chembl_url(resource, p)
    out = []
    key = None
    while url:
        data, _ = fetch_json(url, cache_dir, prefix)
        if key is None:
            key = [k for k in data if k != "page_meta"][0]
        out.extend(data[key])
        nxt = data["page_meta"].get("next")
        url = (CHEMBL_HOST + nxt) if nxt and nxt.startswith("/") else nxt
    return out
