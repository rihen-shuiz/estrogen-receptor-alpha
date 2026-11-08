"""SV01c (RUN_SV1028) step 8 on the proposed H12 window: medoid H12 clash test, run because the SV01c verdict is "mixed".

Runs scripts/sv/SV01b_clash.py unchanged (its main() and rules: Kabsch fit of the medoid C-alpha onto the target's
core v3 C-alpha, moved heavy atoms of the medoid H12, lattice neighbours = other chains' ASU copies and every
symmetry copy with the chain itself and its own bound peptide copy left out, clash = heavy atoms within 3.0 A,
'blocked' = >= 1 backbone atom of the moved H12 in clash), with the window parameters fixed before the SV01c models
were fitted (scripts/sv/SV01c_window.py docstring):
- targets = every chain (all clusters) whose H12 is partial or not modelled on the proposed window 538-544
  (results/sv/SV01c/chain_outcomes.csv, h12_status_538_544), instead of SV01a's 538-547 status;
- moved atoms = the medoid residues 538-544 (SV01b_clash.H12 = proposed window) of 5TM1|A|1 (cl7) and 7UJ7|C|1 (cl4);
- gate = results/sv/SV01c/numbers.json run_step8 (SV01c step verdict "mixed");
- lattice neighbour set and the chain's own peptide copy = results/sv/SV01b/lattice_contacts_per_chain.csv (step 6,
  reused, not recomputed).
SV01b_clash reads its two status/lattice inputs by fixed paths, so the pandas it uses is given a read_csv that serves
those two paths from the SV01c / SV01b files above; nothing else in it changes, and nothing is written to SV01a/SV01b.
Its output column n_obs_538_547 is renamed to the window it holds (n_obs_538_544).
Distance field: common/lattice_contacts.lattice_hits reports m.pos.dist(point), the distance to the neighbour mark's
stored image, which is not the contact distance when the neighbour is found across a lattice translation (gemmi's
find_atoms selects hits by the distance to the nearest image). For this run lattice_hits is replaced by the same
function with the distance taken as ns.dist(point, m.pos) (nearest image); the set of hits is checked to be identical
to the shared function's, so clash counts and 'blocked' are those of the shared code and only min_distance_A changes.
The shared module is not edited (clash_checks.json: distance_fix).
Outputs: results/sv/SV01c/clash_test.csv, clash_summary.csv, clash_checks.json; the step-8 result replaces the
last line of results/sv/SV01c/verdict.md (so run SV01c_window.py, then this script).
Run: .venv/bin/python -B scripts/sv/SV01c_clash.py (after SV01c_window.py)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
import SV01b_clash as C  # noqa: E402
from common import lattice_contacts as lc  # noqa: E402

OUT = ROOT / "results/sv/SV01c"
SV01B = ROOT / "results/sv/SV01b"
NUMS = json.loads((OUT / "numbers.json").read_text())
WIN = tuple(NUMS["window_proposed"])
WTAG = f"{WIN[0]}_{WIN[1]}"
SV01A_OUTCOMES = (ROOT / "results/sv/SV01a/chain_outcomes.csv").resolve()
LATTICE_AT_OUT = (OUT / "lattice_contacts_per_chain.csv").resolve()


class _PandasForClash:
    """pandas as seen by SV01b_clash: read_csv serves the H12 status on the proposed window and the SV01b lattice
    table; every other attribute is pandas itself."""

    def __getattr__(self, name):
        return getattr(pd, name)

    @staticmethod
    def read_csv(path, *args, **kw):
        p = Path(path).resolve()
        if p == SV01A_OUTCOMES:
            co = pd.read_csv(OUT / "chain_outcomes.csv")
            return co.rename(columns={f"h12_status_{WTAG}": "h12_status", f"n_obs_{WTAG}": "n_obs_window"})[
                ["chain_key", "h12_status", "cluster_group", "n_obs_window"]]
        if p == LATTICE_AT_OUT:
            return pd.read_csv(SV01B / "lattice_contacts_per_chain.csv", *args, **kw)
        return pd.read_csv(path, *args, **kw)


_ORIG_HITS = lc.lattice_hits
_FIX = dict(n_calls=0, n_hits=0, n_hits_stored_image_dist_gt_cutoff=0, max_nearest_image_dist_A=0.0,
            max_stored_image_dist_A=0.0, hit_sets_identical=True)


def _lattice_hits_nearest_image(model, ns, point, self_auth, kinds, exclude=None, cutoff=lc.CUTOFF_A):
    """lattice_contacts.lattice_hits with distance = ns.dist(point, mark) (nearest lattice image)."""
    out = []
    for m in ns.find_atoms(point, "\0", radius=cutoff):
        cra = m.to_cra(model)
        key = lc._copy_key(m, cra, point, cutoff)
        if key[2] == "asu" and key[0] == self_auth:
            continue
        if exclude is not None and key == exclude:
            continue
        d, d_stored = float(ns.dist(point, m.pos)), float(m.pos.dist(point))
        _FIX["n_hits"] += 1
        _FIX["n_hits_stored_image_dist_gt_cutoff"] += int(d_stored > cutoff)
        _FIX["max_nearest_image_dist_A"] = max(_FIX["max_nearest_image_dist_A"], d)
        _FIX["max_stored_image_dist_A"] = max(_FIX["max_stored_image_dist_A"], d_stored)
        out.append(key + (kinds.get(key[0], "other"), d))
    orig = _ORIG_HITS(model, ns, point, self_auth, kinds, exclude=exclude, cutoff=cutoff)
    _FIX["n_calls"] += 1
    if [h[:4] for h in orig] != [h[:4] for h in out]:
        _FIX["hit_sets_identical"] = False
    return out


def main():
    assert NUMS["run_step8"], "SV01c verdict not mixed: step 8 does not run"
    C.OUT, C.H12, C.pd = OUT, WIN, _PandasForClash()
    lc.lattice_hits = _lattice_hits_nearest_image
    try:
        C.main()
    finally:
        lc.lattice_hits = _ORIG_HITS
    assert _FIX["hit_sets_identical"] and _FIX["max_nearest_image_dist_A"] <= C.CLASH_A + 1e-6, _FIX
    ct = pd.read_csv(OUT / "clash_test.csv")
    ct.rename(columns={"n_obs_538_547": f"n_obs_{WTAG}"}).to_csv(OUT / "clash_test.csv", index=False)
    chk = json.loads((OUT / "clash_checks.json").read_text())
    chk.update(window=list(WIN), targets_from=f"results/sv/SV01c/chain_outcomes.csv h12_status_{WTAG}",
               lattice_from="results/sv/SV01b/lattice_contacts_per_chain.csv", distance_fix=_FIX)
    (OUT / "clash_checks.json").write_text(json.dumps(chk, indent=1))
    # step-8 result into the last line of verdict.md (written by SV01c_window.py; <= 10 lines)
    cs = pd.read_csv(OUT / "clash_summary.csv")
    k = cs[cs.cluster_group == "all"].set_index("measure").k
    vlines = (OUT / "verdict.md").read_text().splitlines()
    assert vlines[-1].endswith("Step 8 runs (verdict mixed).")
    vlines[-1] = vlines[-1].replace(
        "Step 8 runs (verdict mixed).",
        f"Step 8 (verdict mixed): of {chk['n_targets']} chains partial / not modelled on {WIN[0]}-{WIN[1]}, medoid H12 "
        f"blocked at both positions {int(k['both positions blocked'])}, only cl7 {int(k['only cl7 position blocked'])}, "
        f"only cl4 {int(k['only cl4 position blocked'])}, neither {int(k['neither position blocked'])} (clash_summary.csv).")
    (OUT / "verdict.md").write_text("\n".join(vlines) + "\n")


if __name__ == "__main__":
    main()
