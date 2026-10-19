"""Diagnostic for Q18: what sits near each apo_or_unclear chain (read-only, writes nothing).

For every entry with an apo_or_unclear chain in results/chain_ligand_of_interest.csv, prints the
SUBJECT OF INVESTIGATION components and every non-excluded non-polymer with the number of
pocket residues it contacts, per eligible ESR1 chain. Use it to judge the apo rows by hand.

Run from anywhere: python scripts/check_apo_chains.py
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import ligand_of_interest as loi  # noqa: E402
from esr_analysis.read_structure import category_rows, read_cif  # noqa: E402


def main():
    out = pd.read_csv(loi.OUT_CSV, dtype=str, keep_default_na=False)
    chains = pd.read_csv(loi.CHAINS_CSV, dtype=str, keep_default_na=False)
    pockets = loi.pocket_map()
    for entry_id in sorted(set(out.entry_id[out.status == "apo_or_unclear"])):
        data, atoms = read_cif(ROOT / "data/raw" / f"{entry_id}.cif")
        atoms = loi.heavy(atoms)
        flags = sorted({r["comp_id"] for r in category_rows(data, "_pdbx_entity_instance_feature")
                        if r.get("feature_type") == loi.SUBJECT_FLAG})
        nonpoly = atoms[atoms.label_seq_id.isin([".", "?"]) & ~atoms.label_comp_id.isin(loi.EXCLUDED)]
        comps = sorted(set(nonpoly.label_comp_id))
        print(f"== {entry_id} flagged={flags} non-excluded non-polymers={comps}")
        selected = chains[(chains.entry_id == entry_id) & (chains.identity == "esr1")
                          & (chains.eligible_historical == "True")]
        for chain in selected.to_dict("records"):
            pocket = pockets.get(chain["chain_key"], {})
            pocket_atoms = atoms[(atoms.label_asym_id == chain["chain_id"])
                                 & atoms.label_seq_id.isin(list(pocket))]
            pocket_xyz = pocket_atoms[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float)
            pocket_positions = [pocket[s][0] for s in pocket_atoms.label_seq_id]
            found = []
            for (asym, comp, auth_seq, ins), group in nonpoly.groupby(
                    ["label_asym_id", "label_comp_id", "auth_seq_id", "pdbx_PDB_ins_code"]):
                contacts = loi.contacted_positions(group[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float),
                                                   pocket_xyz, pocket_positions)
                found.append(f"{comp}[{asym}:{auth_seq}] contacts={len(contacts)}")
            print(f"   {chain['chain_key']} pocket_atoms={len(pocket_atoms)} | "
                  + ("; ".join(found) if found else "no non-excluded non-polymer"))


if __name__ == "__main__":
    main()
