"""Ligand of interest per eligible ESR1 chain (tracker task Q18).

Input: results/lbd_chains.csv (chains), results/lbd_extract/residue_mapping.csv
(canonical ESR1 position per observed residue) and the deposited mmCIF files in data/raw/.
Output: results/chain_ligand_of_interest.csv, one row per chain and ligand. Every chain is
kept. A chain without a pocket ligand gets an empty ligand_id and status apo_or_unclear.
ligand_key is entry|label_asym_id|auth_seq_id+insertion code of the ligand residue.

Rules, per chain (heavy atoms only):
1. Flag: the entry's pdbx_entity_instance_feature rows mark components as 'SUBJECT OF
   INVESTIGATION'. The flag gives only the component ID. A flagged non-polymer with a heavy
   atom within 4.5 A of a pocket residue of this chain is this chain's ligand (source
   subject_of_investigation). Each such non-polymer gets one row.
2. Otherwise: a non-polymer with heavy atoms within 4.5 A of at least 3 pocket residues of
   this chain is the ligand (source pocket_contact).
Waters, ions, crystallization additives and buffers never count as a ligand (EXCLUDED).

Run from anywhere: python scripts/ligand_of_interest.py
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from esr_analysis.read_structure import category_rows, read_cif  # noqa: E402

# Pocket residues in canonical ESR1 numbering. Source: task Q18 in tracker/esr1.tasks
# (project description). Not re-checked against the literature here.
POCKET_POSITIONS = (343, 346, 347, 350, 353, 383, 384, 387, 388, 391, 394, 404, 421, 424, 425, 428, 521, 524, 525)
CONTACT_CUTOFF_A = 4.5
MIN_POCKET_CONTACTS = 3
SUBJECT_FLAG = "SUBJECT OF INVESTIGATION"
WATERS = {"HOH", "DOD", "WAT", "H2O"}
IONS = {"NA", "K", "MG", "CA", "ZN", "CL", "BR", "IOD", "CD", "MN", "FE", "FE2", "CU", "CU1", "CO",
        "NI", "LI", "RB", "CS", "SR", "BA", "HG", "PT", "AU", "AG", "TL", "PB", "EU", "TB", "GD",
        "LA", "YB", "IR", "OS", "RU", "RH", "PD", "AL", "GA", "IN", "NO3"}
# Includes the Q16 is_additive list (results/ligands_ranked.csv) plus cryoprotectants and buffers.
ADDITIVES = {"GOL", "EDO", "PEG", "PG4", "PE4", "PE8", "P6G", "PGE", "1PE", "2PE", "MPD", "SO4", "PO4",
             "DMS", "DMSO", "ACY", "ACT", "FMT", "BME", "DTT", "MES", "TRS", "EPE", "CIT", "TAR", "TLA",
             "MLI", "IMD", "SCN", "BCT", "CAC"}
EXCLUDED = WATERS | IONS | ADDITIVES

CHAINS_CSV = ROOT / "results/lbd_chains.csv"
MAPPING_CSV = ROOT / "results/lbd_extract/residue_mapping.csv"
OUT_CSV = ROOT / "results/chain_ligand_of_interest.csv"
SUBJECT, POCKET = "subject_of_investigation", "pocket_contact"


def heavy(atoms):
    return atoms[~atoms.type_symbol.str.upper().isin(["H", "D"])]


def pocket_map():
    """chain_key -> {label_seq_id: (canonical position, one-letter residue)} for observed pocket residues."""
    cols = ["chain_key", "label_seq_id", "canonical_position", "observed_aa", "has_residue", "identity"]
    rm = pd.read_csv(MAPPING_CSV, dtype=str, keep_default_na=False, usecols=cols)
    positions = {str(p) for p in POCKET_POSITIONS}
    rm = rm[(rm.identity == "esr1") & (rm.has_residue == "True") & rm.canonical_position.isin(positions)]
    out = {}
    for row in rm.itertuples(index=False):
        out.setdefault(row.chain_key, {})[row.label_seq_id] = (int(row.canonical_position), row.observed_aa)
    return out


def contacted_positions(ligand_xyz, pocket_xyz, pocket_positions):
    """Canonical pocket positions with a heavy atom within the cutoff of any ligand heavy atom."""
    if len(pocket_xyz) == 0:
        return []
    distance = np.linalg.norm(ligand_xyz[:, None, :] - pocket_xyz[None, :, :], axis=2)
    near = distance.min(axis=0) <= CONTACT_CUTOFF_A
    return sorted({pocket_positions[j] for j in np.flatnonzero(near)})


def rows_for_entry(entry_id, chains, pockets):
    data, atoms = read_cif(ROOT / "data/raw" / f"{entry_id}.cif")
    atoms = heavy(atoms)
    nonpoly = atoms[atoms.label_seq_id.isin([".", "?"]) & ~atoms.label_comp_id.isin(EXCLUDED)]
    instances = {}
    for key, group in nonpoly.groupby(["label_asym_id", "label_comp_id", "auth_seq_id", "pdbx_PDB_ins_code"], sort=True):
        instances[key] = group[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float)
    flagged = {r["comp_id"] for r in category_rows(data, "_pdbx_entity_instance_feature")
               if r.get("feature_type") == SUBJECT_FLAG} - EXCLUDED

    rows = []
    for chain in chains:
        pocket = pockets.get(chain["chain_key"], {})
        chain_atoms = atoms[(atoms.label_asym_id == chain["chain_id"]) & (atoms.pdbx_PDB_model_num == chain["model"])
                            & atoms.label_seq_id.isin(list(pocket))]
        pocket_xyz = chain_atoms[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float)
        pocket_positions = [pocket[s][0] for s in chain_atoms.label_seq_id]
        residue_letter = {pocket[s][0]: pocket[s][1] for s in chain_atoms.label_seq_id}
        contacts = {key: contacted_positions(xyz, pocket_xyz, pocket_positions) for key, xyz in instances.items()}
        picks = [(key, found, SUBJECT) for key, found in contacts.items() if key[1] in flagged and found]
        if not picks:
            picks = [(key, found, POCKET) for key, found in contacts.items() if len(found) >= MIN_POCKET_CONTACTS]
        base = {"chain_key": chain["chain_key"], "entry_id": entry_id}
        if not picks:
            rows.append({**base, "ligand_id": "", "ligand_key": "", "source": "", "n_pocket_residues": 0,
                         "pocket_residues": "", "status": "apo_or_unclear"})
            continue
        status = "multi" if len(picks) > 1 else "single"
        for (asym, comp, auth_seq, ins), found, source in sorted(picks):
            ins = "" if ins in ("?", ".") else ins
            rows.append({**base, "ligand_id": comp, "ligand_key": f"{entry_id}|{asym}|{auth_seq}{ins}",
                         "source": source, "n_pocket_residues": len(found),
                         "pocket_residues": ";".join(f"{residue_letter[p]}{p}" for p in found),
                         "status": status})
    return rows


def main():
    chains = pd.read_csv(CHAINS_CSV, dtype=str, keep_default_na=False)
    eligible = chains[(chains.identity == "esr1") & (chains.eligible_historical == "True")]
    pockets = pocket_map()
    rows = []
    for entry_id, group in eligible.groupby("entry_id", sort=True):
        rows += rows_for_entry(entry_id, group.to_dict("records"), pockets)
    out = pd.DataFrame(rows, columns=["chain_key", "entry_id", "ligand_id", "ligand_key", "source",
                                      "n_pocket_residues", "pocket_residues", "status"])
    out = out.sort_values(["chain_key", "ligand_id"], kind="stable")
    out.to_csv(OUT_CSV, index=False)
    print(f"eligible ESR1 chains: {eligible.chain_key.nunique()}  entries: {eligible.entry_id.nunique()}")
    print(f"rows: {len(out)}  chains: {out.chain_key.nunique()}")
    print(out.status.value_counts().to_string())
    print(out.source.value_counts().to_string())
    print(f"wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
