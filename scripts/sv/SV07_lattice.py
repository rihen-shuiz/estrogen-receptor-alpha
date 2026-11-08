"""SV07 (RUN_SV1028) step 2, lattice and partner contacts at the spike positions (SV01b module, unchanged).

Spike positions: the step-1 rule of scripts/sv/SV07_spikes.py (spike_positions) applied to the Y10 card output
results/ext/cluster_displacement_profile.csv.
Contacts: scripts/sv/common/lattice_contacts.py (SV01b step 6 definitions: heavy atoms <= 4.0 A; ASU copy of another
chain = 'asu' (partner chains of the asymmetric unit), any other symmetry operation or lattice translation = 'sym'
(lattice); the chain's own bound peptide copy (SV00 rule) reported separately as own_peptide_atoms), for every kept
chain of results/sv/SV00/chain_table.csv at those positions.
Checks: residues found vs the SV00 C-alpha observation mask at the spike positions.
Outputs: results/sv/SV07/lattice_contacts_spikes_per_residue.csv, lattice_contacts_spikes_per_chain.csv,
lattice_checks.json.
Run: .venv/bin/python -B scripts/sv/SV07_lattice.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from common import lattice_contacts as lc  # noqa: E402
from SV07_spikes import spike_positions  # noqa: E402

OUT = ROOT / "results/sv/SV07"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    y10 = pd.read_csv(ROOT / "results/ext/cluster_displacement_profile.csv")
    pos = spike_positions(y10).position.astype(int).tolist()
    t = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", low_memory=False)
    lookup = lc.residue_lookup(t.chain_key)
    res_parts, chain_parts = [], []
    t0 = time.time()
    for i, (ent, g) in enumerate(t.groupby("entry_id", sort=False)):
        r, c = lc.entry_contacts(ent, list(zip(g.chain_key, g.auth_asym_id)), pos, lookup=lookup)
        res_parts.append(r)
        chain_parts.append(c)
        if i % 50 == 0:
            print(f"{i} entries, {time.time() - t0:.0f} s", flush=True)
    per_res = pd.concat(res_parts, ignore_index=True)
    per_chain = pd.concat(chain_parts, ignore_index=True).set_index("chain_key").loc[t.chain_key].reset_index()
    per_res.to_csv(OUT / "lattice_contacts_spikes_per_residue.csv", index=False)
    per_chain.to_csv(OUT / "lattice_contacts_spikes_per_chain.csv", index=False)
    mask = np.load(ROOT / "results/sv/SV00/observed_mask.npy")
    cols = pd.read_csv(ROOT / "results/sv/SV00/observed_mask_cols.csv").esr1_position.values
    w = np.isin(cols, pos)
    obs = {(k, int(p)) for i, k in enumerate(t.chain_key) for p in cols[w][mask[i, w]]}
    found = set(zip(per_res.chain_key, per_res.esr1_position.astype(int)))
    checks = dict(spike_positions=pos, n_chains=len(per_chain), n_entries=int(per_chain.entry_id.nunique()),
                  chains_not_in_model=int((~per_chain.chain_in_model).sum()),
                  residues_found=len(found), residues_ca_observed_mask=len(obs),
                  found_not_ca_observed=sorted(f"{k}:{p}" for k, p in found - obs),
                  ca_observed_not_found=sorted(f"{k}:{p}" for k, p in obs - found),
                  seconds=round(time.time() - t0, 1))
    (OUT / "lattice_checks.json").write_text(json.dumps(checks, indent=1))
    print(json.dumps({k: v for k, v in checks.items() if not isinstance(v, list) or len(v) < 30}, indent=1))


if __name__ == "__main__":
    main()
