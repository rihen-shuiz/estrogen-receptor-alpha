"""What did the deposition declare and what has coordinates? Read mmCIF
categories and atom records, preserving label/auth IDs, model and insertion
codes. Select alternate atoms by highest occupancy, then blank/A/name ties.
These are deposited records, not a biological interpretation of chain counts.
"""
import numpy as np
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def category_rows(data, category):
    prefix = category + "."
    columns = {}
    for key, value in data.items():
        if key.startswith(prefix):
            columns[key[len(prefix):]] = value if isinstance(value, list) else [value]
    if not columns:
        return []
    count = len(next(iter(columns.values())))
    if any(len(values) != count for values in columns.values()):
        raise ValueError(f"Unequal column lengths in mmCIF {category}")
    return [{name: values[i] for name, values in columns.items()} for i in range(count)]


def atom_records(data, select_alternates=True):
    """Preserve all conformers for contacts; select one atom for CA mapping."""
    atoms = pd.DataFrame(category_rows(data, "_atom_site"))
    if atoms.empty:
        raise ValueError("No coordinate records in mmCIF")
    for column in ["Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv"]:
        atoms[column] = pd.to_numeric(atoms[column], errors="coerce")
    if not np.isfinite(atoms[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy()).all():
        raise ValueError("Nonfinite atom coordinates in mmCIF")
    atoms["pdbx_PDB_model_num"] = atoms.get("pdbx_PDB_model_num", "1").astype(str)
    if not select_alternates:
        return atoms
    atoms["alt_rank"] = atoms["label_alt_id"].replace({".": "0", "?": "0", "A": "1"})
    atoms = atoms.sort_values(["occupancy", "alt_rank"], ascending=[False, True], kind="stable")
    # Nonpolymers have no label_seq_id: preserve their author residue/insertion identity.
    atom_identity = ["pdbx_PDB_model_num", "label_asym_id", "auth_seq_id", "pdbx_PDB_ins_code", "label_atom_id"]
    atoms = atoms.drop_duplicates(atom_identity).sort_index().drop(columns="alt_rank")
    return atoms


def read_cif(path):
    data = MMCIF2Dict(str(path))
    return data, atom_records(data)


def first_value(data, names):
    for name in names:
        for value in data.get(name, []):
            if value not in (".", "?", ""):
                return value, name
    return None, "missing"
