"""SV00 step 3-4: regenerate the Q04 LBD extraction from the restored mmCIFs and extract the kept-chain backbone.

Code path (unchanged, imported): esr_analysis.read_structure.read_cif (altlocs: highest occupancy, then blank/A) and
esr_analysis.extract_lbd.inventory_entry, called exactly as scripts/rebuild_lbd_extract.py calls them (Q04 rebuild,
tracker Y31), for all 526 frozen IDs (ESR1 list then ESR2 list). Only the output folder differs: data/sv/lbd_extract/
instead of results/ (read-only, rule 1.2).
Backbone: for each of the 788 kept chains (results/rmsd_entry_index.csv, matrix order) the N, CA, C, O atoms of the
altloc-selected atom table, usable occupancy as in extract_lbd.residue_observations (occupancy > 0 or unknown),
deposited (crystal) frame, numbered with the regenerated residue mapping (label_seq_id -> ESR1/UniProt position).
Outputs: data/sv/lbd_extract/{entries,receptor_chains,residue_mapping}.csv, data/sv/lbd_chains.csv,
data/sv/backbone/backbone.npz (coords [chain, position 305-550, atom N/CA/C/O, xyz], NaN = not observed; B, occupancy),
data/sv/backbone/backbone_atoms.csv.gz (long table, one row per atom, with label/auth numbering).
Run: .venv/bin/python -I scripts/sv/SV00_extract.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from esr_analysis._task_io import RESOURCE_ROOT, sequences  # noqa: E402
from esr_analysis.extract_lbd import inventory_entry  # noqa: E402
from esr_analysis.map_residues import receptor_equivalence  # noqa: E402
from esr_analysis.read_structure import read_cif  # noqa: E402
from esr_analysis.run_io import sha256  # noqa: E402
from fetch_data import read_ids  # noqa: E402

# the settings file the Q04 rebuild (scripts/rebuild_common.SETTINGS) reads
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
OUT = ROOT / "data/sv"
ATOMS = ["N", "CA", "C", "O"]
LBD = np.arange(SETTINGS["esr1_lbd_start"], SETTINGS["esr1_lbd_end"] + 1)


def backbone_rows(chain_key, atoms, mapping_rows):
    pdb, label, model = chain_key.split("|")
    to_e1 = {int(r["label_seq_id"]): r["esr1_position"] for r in mapping_rows
             if r["chain_key"] == chain_key and r["esr1_position"] is not None and np.isfinite(r["esr1_position"])}
    a = atoms[(atoms.label_asym_id == label) & (atoms.pdbx_PDB_model_num == model) & atoms.label_atom_id.isin(ATOMS)]
    a = a[a.occupancy.isna() | (a.occupancy > 0)]
    a = a[a.label_seq_id.astype(str).str.isdigit()]
    rows = []
    for at in a.itertuples():
        seq = int(at.label_seq_id)
        e1 = to_e1.get(seq)
        rows.append({"chain_key": chain_key, "entry_id": pdb, "label_seq_id": seq,
                     "esr1_position": int(e1) if e1 is not None else np.nan,
                     "auth_seq_id": at.auth_seq_id, "ins_code": "" if at.pdbx_PDB_ins_code in ("?", ".") else at.pdbx_PDB_ins_code,
                     "resname": at.label_comp_id, "atom": at.label_atom_id,
                     "altloc": "" if at.label_alt_id in ("?", ".") else at.label_alt_id,
                     "x": at.Cartn_x, "y": at.Cartn_y, "z": at.Cartn_z, "occupancy": at.occupancy, "B": at.B_iso_or_equiv})
    return rows


def main():
    seq = sequences(RESOURCE_ROOT, SETTINGS)
    equivalence = receptor_equivalence(seq)
    kept = pd.read_csv(ROOT / "results/rmsd_entry_index.csv").chain_key.tolist()
    kept_by_entry = {}
    for key in kept:
        kept_by_entry.setdefault(key.split("|")[0], []).append(key)
    entries, chains, mapping, bb = [], [], [], []
    ids = read_ids()
    for i, pdb in enumerate(ids, 1):
        path = ROOT / "data/raw" / f"{pdb}.cif"
        data, atoms = read_cif(path)
        entry, local_chains, local_mapping = inventory_entry(pdb, data, atoms, seq, equivalence, SETTINGS)
        entry.update(source_path=path.relative_to(ROOT).as_posix(), source_sha256=sha256(path), query_memberships="esr1")
        entries.append(entry)
        chains.extend(local_chains)
        mapping.extend(local_mapping)
        for key in kept_by_entry.get(pdb, []):
            bb.extend(backbone_rows(key, atoms, local_mapping))
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

    bbd = pd.DataFrame(bb)
    (OUT / "backbone").mkdir(parents=True, exist_ok=True)
    bbd.to_csv(OUT / "backbone/backbone_atoms.csv.gz", index=False)
    coords = np.full((len(kept), len(LBD), len(ATOMS), 3), np.nan)
    bfac = np.full((len(kept), len(LBD), len(ATOMS)), np.nan)
    occ = np.full((len(kept), len(LBD), len(ATOMS)), np.nan)
    row_of = {k: i for i, k in enumerate(kept)}
    inlbd = bbd[bbd.esr1_position.between(LBD[0], LBD[-1])]
    dup = inlbd.duplicated(["chain_key", "esr1_position", "atom"]).sum()
    if dup:
        raise ValueError(f"{dup} duplicated backbone atoms after altloc selection")
    r = inlbd.chain_key.map(row_of).to_numpy()
    p = (inlbd.esr1_position.astype(int) - LBD[0]).to_numpy()
    a = inlbd.atom.map({x: j for j, x in enumerate(ATOMS)}).to_numpy()
    coords[r, p, a] = inlbd[["x", "y", "z"]].to_numpy(float)
    bfac[r, p, a] = inlbd.B.to_numpy(float)
    occ[r, p, a] = inlbd.occupancy.to_numpy(float)
    np.savez_compressed(OUT / "backbone/backbone.npz", chain_keys=np.asarray(kept), positions=LBD,
                        atoms=np.asarray(ATOMS), coords=coords, b_factor=bfac, occupancy=occ,
                        frame=np.asarray("deposited crystal frame (no superposition)"))
    print(f"SV00 extract: {len(entries)} entries, {len(chain_table)} chain/models, {len(lbd)} receptor chains; "
          f"backbone {len(kept)} chains, {int(np.isfinite(coords[..., 0]).sum())} LBD backbone atoms")


if __name__ == "__main__":
    main()
