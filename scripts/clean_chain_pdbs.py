"""Receptor chain plus its ligand of interest: one clean PDB per eligible ESR1 chain (tracker task S17).

Input: results/lbd_chains.csv (eligible ESR1 chains), results/chain_ligand_of_interest.csv
(ligand_key per chain, from task Q18) and the deposited mmCIF files in data/raw/.
Output: data/clean/<entry>_<chain>.pdb for every eligible chain. A file holds the chain's
polymer atoms and its ligand instance(s) as HETATM, with deposited atom names, author chain
and residue numbers, coordinates, occupancy and B factors. Waters, ions, additives and
partner chains are not written. A chain without a ligand of interest gets only its protein.
PDB residue names have 3 characters: a longer CCD ID uses its last 3 characters, recorded in
data/clean/resname_map.csv.

Run from anywhere: python scripts/clean_chain_pdbs.py
"""
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from esr_analysis.read_structure import read_cif  # noqa: E402

CHAINS_CSV = ROOT / "results/lbd_chains.csv"
LIGANDS_CSV = ROOT / "results/chain_ligand_of_interest.csv"
OUT_DIR = ROOT / "data/clean"
MISSING = ("?", ".")


def ligand_resname(ligand_id):
    return ligand_id if len(ligand_id) <= 3 else ligand_id[-3:]


def residue_ins(value):
    return "" if value in MISSING else value


def pdb_atom(serial, row, record, chain_char, resname, ins):
    """One fixed-column PDB ATOM/HETATM record (PDB format v3.3 column positions)."""
    element = row.type_symbol.strip()
    name = row.auth_atom_id
    if len(name) >= 4:
        atom4 = name[:4]
    elif len(element) == 2:
        atom4 = name.ljust(4)
    else:
        atom4 = (" " + name).ljust(4)
    occupancy = 1.0 if pd.isna(row.occupancy) else row.occupancy
    bfactor = 0.0 if pd.isna(row.B_iso_or_equiv) else row.B_iso_or_equiv
    return (f"{record:<6}{serial:>5} {atom4} {resname:>3} {chain_char}{int(row.auth_seq_id):>4}{ins:1}   "
            f"{row.Cartn_x:8.3f}{row.Cartn_y:8.3f}{row.Cartn_z:8.3f}{occupancy:6.2f}{bfactor:6.2f}"
            f"          {element:>2}")


def main():
    chains = pd.read_csv(CHAINS_CSV, dtype=str, keep_default_na=False)
    eligible = chains[(chains.identity == "esr1") & (chains.eligible_historical == "True")]
    ligands = pd.read_csv(LIGANDS_CSV, dtype=str, keep_default_na=False)
    ligands = ligands[ligands.ligand_id != ""]
    by_chain = {key: group.to_dict("records") for key, group in ligands.groupby("chain_key")}
    resnames = {ligand_id: ligand_resname(ligand_id) for ligand_id in sorted(set(ligands.ligand_id))}
    if len(set(resnames.values())) != len(resnames):
        raise ValueError("ligand residue names collide after shortening to 3 characters")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"ligand_id": list(resnames), "pdb_resname": list(resnames.values())}).to_csv(
        OUT_DIR / "resname_map.csv", index=False)

    written = with_ligand = ligand_atoms = 0
    for entry_id, group in eligible.groupby("entry_id", sort=True):
        _, atoms = read_cif(ROOT / "data/raw" / f"{entry_id}.cif")
        ins_key = atoms.pdbx_PDB_ins_code.map(residue_ins)
        res_key = atoms.auth_seq_id + ins_key
        for chain in group.to_dict("records"):
            protein = atoms[(atoms.label_asym_id == chain["label_asym_id"])
                            & (atoms.pdbx_PDB_model_num == chain["model"])
                            & ~atoms.label_seq_id.isin(MISSING)]
            if protein.empty:
                raise ValueError(f"{chain['chain_key']}: no polymer atoms in the mmCIF")
            chain_char = chain["auth_asym_id"]
            lines, serial, last = [], 1, None
            for row in protein.itertuples(index=False):
                ins = residue_ins(row.pdbx_PDB_ins_code)
                resname = row.auth_comp_id
                lines.append(pdb_atom(serial, row, row.group_PDB, chain_char, resname, ins))
                serial += 1
                last = (resname, ins, row.auth_seq_id)
            resname, ins, seq = last
            lines.append(f"TER   {serial:>5}      {resname:>3} {chain_char}{int(seq):>4}{ins:1}")
            serial += 1
            for ligand in by_chain.get(chain["chain_key"], []):
                _, asym, auth_ins = ligand["ligand_key"].split("|")
                lig = atoms[(atoms.label_asym_id == asym) & (atoms.label_comp_id == ligand["ligand_id"])
                            & (res_key == auth_ins)]
                if lig.empty:
                    raise ValueError(f"{ligand['ligand_key']}: ligand atoms not found in mmCIF")
                for row in lig.itertuples(index=False):
                    lines.append(pdb_atom(serial, row, "HETATM", chain_char,
                                          resnames[ligand["ligand_id"]], residue_ins(row.pdbx_PDB_ins_code)))
                    serial += 1
                ligand_atoms += len(lig)
            lines.append("END")
            (OUT_DIR / f"{entry_id}_{chain['chain_id']}.pdb").write_text("\n".join(lines) + "\n")
            written += 1
            with_ligand += chain["chain_key"] in by_chain

    print(f"wrote {written} PDB files to {OUT_DIR.relative_to(ROOT)}; chains with ligand: {with_ligand}; "
          f"ligand atoms written: {ligand_atoms}")


if __name__ == "__main__":
    main()
