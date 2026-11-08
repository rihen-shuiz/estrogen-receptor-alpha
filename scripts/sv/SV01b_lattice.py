"""SV01b step 6: lattice contacts of ESR1 residues 528-550 for every kept chain (module scripts/sv/common/lattice_contacts.py).

Input: results/sv/SV00/chain_table.csv (788 kept chains; auth_asym_id), restored mmCIFs data/raw/ (SV00), regenerated
Q04 tables data/sv/lbd_extract/ (SV00). Definitions are in the module docstring (contact <= 4.0 A heavy atoms, ASU copy
vs symmetry copy, the chain's own bound peptide copy reported separately and not counted as lattice).
Checks: own bound peptide vs SV00 chain_peptide.csv (auth chain and source); residues found vs the SV00 observation
mask (C-alpha observed) at 528-550.
Outputs: results/sv/SV01b/lattice_contacts_per_residue.csv, lattice_contacts_per_chain.csv, lattice_checks.json.
Run: .venv/bin/python -B scripts/sv/SV01b_lattice.py
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

OUT = ROOT / "results/sv/SV01b"
OUT.mkdir(parents=True, exist_ok=True)
WINDOW = (528, 550)


def main():
    t = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", low_memory=False)
    lookup = lc.residue_lookup(t.chain_key)
    res_parts, chain_parts = [], []
    t0 = time.time()
    for i, (ent, g) in enumerate(t.groupby("entry_id", sort=False)):
        r, c = lc.entry_contacts(ent, list(zip(g.chain_key, g.auth_asym_id)), range(WINDOW[0], WINDOW[1] + 1),
                                 lookup=lookup)
        res_parts.append(r)
        chain_parts.append(c)
        if i % 50 == 0:
            print(f"{i} entries, {time.time() - t0:.0f} s", flush=True)
    per_res = pd.concat(res_parts, ignore_index=True)
    per_chain = pd.concat(chain_parts, ignore_index=True)
    summ = lc.summarise(per_res, *WINDOW, prefix="")
    per_chain = per_chain.merge(summ, on="chain_key", how="left")
    for c in summ.columns:
        if c != "chain_key":
            per_chain[c] = per_chain[c].fillna(0).astype(int)
    per_chain = per_chain.set_index("chain_key").loc[t.chain_key].reset_index()
    per_res.to_csv(OUT / "lattice_contacts_per_residue.csv", index=False)
    per_chain.to_csv(OUT / "lattice_contacts_per_chain.csv", index=False)

    # checks
    pep = pd.read_csv(ROOT / "results/sv/SV00/chain_peptide.csv", keep_default_na=False)
    pep = pep.set_index("chain_key").loc[t.chain_key]
    same_pep = ((pep.peptide_auth_asym_id.values == per_chain.own_peptide_auth_asym_id.values)
                & (pep.peptide_source.values == per_chain.own_peptide_source.values))
    mask = np.load(ROOT / "results/sv/SV00/observed_mask.npy")
    cols = pd.read_csv(ROOT / "results/sv/SV00/observed_mask_cols.csv").esr1_position.values
    w = (cols >= WINDOW[0]) & (cols <= WINDOW[1])
    obs = {(k, int(p)) for i, k in enumerate(t.chain_key) for p in cols[w][mask[i, w]]}
    found = set(zip(per_res.chain_key, per_res.esr1_position.astype(int)))
    checks = dict(
        n_chains=len(per_chain), n_entries=int(per_chain.entry_id.nunique()),
        chains_not_in_model=int((~per_chain.chain_in_model).sum()),
        crystal_frame=per_chain.drop_duplicates("entry_id").crystal_frame.value_counts().to_dict(),
        own_peptide_agrees_with_SV00=int(same_pep.sum()),
        own_peptide_disagrees=per_chain.chain_key[~same_pep].tolist(),
        own_peptide_source=per_chain.own_peptide_source.value_counts().to_dict(),
        residues_528_550_found=len(found), residues_528_550_ca_observed_mask=len(obs),
        found_not_ca_observed=sorted(f"{k}:{p}" for k, p in found - obs),
        ca_observed_not_found=sorted(f"{k}:{p}" for k, p in obs - found),
        chains_with_asu_contact=int((per_chain.asu_atoms > 0).sum()),
        chains_with_sym_contact=int((per_chain.sym_atoms > 0).sum()),
        seconds=round(time.time() - t0, 1))
    (OUT / "lattice_checks.json").write_text(json.dumps(checks, indent=1))
    print(json.dumps({k: v for k, v in checks.items() if not isinstance(v, list) or len(v) < 20}, indent=1))


if __name__ == "__main__":
    main()
