"""Rebuild tracker Q04: per-chain LBD extraction (ESR1 residues 305-550) for all 526 ESR1+ESR2 entries.

Inputs: mmCIFs (data/raw, else results/rebuild/raw from scripts/rebuild_fetch_raw.py), ID lists
data/pdb_ids_esr1.txt + data/pdb_ids_esr2.txt (in that order), UniProt sequences and settings from
esr_analysis/resources.
Method: esr_analysis.extract_lbd.inventory_entry for every entry (block 04 logic: identity from deposited
UniProt accession, mapping by deposited segments + alignment, ESR2 numbered onto ESR1 through
map_residues.receptor_equivalence, eligibility = >=180 LBD CA and >=60 core CA). Entry rows get
source_path, source_sha256, query_memberships='esr1' (as the block 04 CLI writes them).
lbd_chains.csv = receptor_chains rows with identity esr1/esr2, prefixed with entry_id (=pdb_id) and
chain_id (=label_asym_id).
Outputs: results/rebuild/lbd_extract/{entries,receptor_chains,residue_mapping}.csv, results/rebuild/lbd_chains.csv
Run from the repo root: PYTHONPATH=. python scripts/rebuild_lbd_extract.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from esr_analysis._task_io import RESOURCE_ROOT, sequences
from esr_analysis.extract_lbd import inventory_entry
from esr_analysis.map_residues import receptor_equivalence
from esr_analysis.read_structure import read_cif
from esr_analysis.run_io import sha256
from rebuild_common import OUT, SETTINGS, raw_path, read_ids, rel


def main():
    seq = sequences(RESOURCE_ROOT, SETTINGS)
    equivalence = receptor_equivalence(seq)
    entries, chains, mapping = [], [], []
    ids = read_ids()
    for i, pdb in enumerate(ids, 1):
        path = raw_path(pdb)
        data, atoms = read_cif(path)
        entry, local_chains, local_mapping = inventory_entry(pdb, data, atoms, seq, equivalence, SETTINGS)
        entry.update(source_path=rel(path), source_sha256=sha256(path), query_memberships="esr1")
        entries.append(entry)
        chains.extend(local_chains)
        mapping.extend(local_mapping)
        if i % 50 == 0:
            print(f"  {i}/{len(ids)}", flush=True)
    out = OUT / "lbd_extract"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(entries).to_csv(out / "entries.csv", index=False)
    chain_table = pd.DataFrame(chains)
    chain_table.to_csv(out / "receptor_chains.csv", index=False)
    pd.DataFrame(mapping).to_csv(out / "residue_mapping.csv", index=False)
    lbd = chain_table[chain_table.identity.isin(["esr1", "esr2"])].copy()
    lbd.insert(0, "chain_id", lbd.label_asym_id)
    lbd.insert(0, "entry_id", lbd.pdb_id)
    lbd.to_csv(OUT / "lbd_chains.csv", index=False)
    print(f"Q04 rebuild: {len(entries)} entries, {len(chain_table)} chain/models, {len(lbd)} receptor chains")


if __name__ == "__main__":
    main()
