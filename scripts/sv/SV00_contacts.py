"""SV00 step 6 (contact part): dimer partner at H11 and the peptide bound to each receptor chain.

Structures: the restored mmCIFs (data/raw), read with gemmi as scripts/ext_h12_contacts.py (Y44) does: waters,
ligands and hydrogens removed, model 1, NeighborSearch with crystal symmetry (5 A cells; non-standard SCALEn ->
standardize_crystal_frame, a rigid change). All conformers are kept for contacts (read_structure.atom_records policy).
Chain roles from the regenerated Q04 table (data/sv/lbd_extract/receptor_chains.csv, model 1, keyed by auth_asym_id):
receptor = identity esr1/esr2; peptide = any other polypeptide of <= peptide_max_length residues (settings, the
existing peptide definition of Y02/Y44); the rest = other.
Copies: an atom found at symmetry image 0 whose deposited position is itself within the cut-off is an ASU copy;
every other hit (other symmetry operation or a lattice translation, including the chain's own image) is a symmetry
copy. Copies are grouped by (chain, image index, asu/symmetry).
Dimer partner (step 6): for each kept chain, the receptor copy (not the chain itself in the ASU) with the most heavy
atoms of the chain's H11 (ESR1 497-530; author numbers from the regenerated residue mapping) within 4.0 A; ties ->
ASU before symmetry, then chain name. Source asu / symmetry / none.
Bound peptide (for 'peptide' per chain and 'partner's peptide'): for every receptor chain of the 351 entries, the
peptide copy with the most receptor heavy atoms within 4.0 A (same tie rule); none if no peptide atom is that close.
Outputs: results/sv/SV00/dimer_partner.csv (one row per kept chain), results/sv/SV00/chain_peptide.csv (one row per
receptor chain of the 351 entries).
Run: .venv/bin/python -I scripts/sv/SV00_contacts.py
"""
import json
import sys
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV00"
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
H11 = (497, 530)
CUTOFF_A = 4.0


def structure(path):
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_ligands_and_waters()
    st.remove_hydrogens()
    model = st[0]
    try:
        ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
        frame = "standard"
    except RuntimeError:
        st.standardize_crystal_frame()
        model = st[0]
        ns = gemmi.NeighborSearch(model, st.cell, 5).populate()
        frame = "standardized"
    return model, ns, frame


def contacts(model, ns, chain_name, residue_filter, kinds, want):
    """{(partner chain, image_idx, source): set of query atoms} for partners whose kind is in `want`."""
    found = {}
    chain = next((ch for ch in model if ch.name == chain_name), None)
    if chain is None:
        return found, False
    for res in chain:
        if residue_filter is not None and res.seqid.num not in residue_filter:
            continue
        for at in res:
            for m in ns.find_atoms(at.pos, "\0", radius=CUTOFF_A):
                cra = m.to_cra(model)
                if kinds.get(cra.chain.name) not in want:
                    continue
                asu = m.image_idx == 0 and cra.atom.pos.dist(at.pos) <= CUTOFF_A + 1e-6
                if asu and cra.chain.name == chain_name:
                    continue
                key = (cra.chain.name, int(m.image_idx), "asu" if asu else "symmetry")
                found.setdefault(key, set()).add((res.seqid.num, res.seqid.icode, at.name, at.altloc))
    return found, True


def best(found):
    if not found:
        return None, 0, 0, 0
    ranked = sorted(found.items(), key=lambda kv: (-len(kv[1]), kv[0][2] != "asu", kv[0][0], kv[0][1]))
    second = len(ranked[1][1]) if len(ranked) > 1 else 0
    return ranked[0][0], len(ranked[0][1]), second, len(ranked)


def main():
    idx = pd.read_csv(ROOT / "results/rmsd_entry_index.csv")
    rc = pd.read_csv(ROOT / "data/sv/lbd_extract/receptor_chains.csv", keep_default_na=False, na_values=[""],
                     low_memory=False)
    rc = rc[rc.model.astype(str) == "1"]
    mp = pd.read_csv(ROOT / "data/sv/lbd_extract/residue_mapping.csv", low_memory=False,
                     usecols=["chain_key", "esr1_position", "has_residue", "auth_seq_id"])
    mp = mp[mp.chain_key.isin(idx.chain_key) & mp.has_residue.astype(bool) & mp.esr1_position.between(*H11)]
    h11 = mp.groupby("chain_key").auth_seq_id.apply(lambda s: set(int(float(x)) for x in s.dropna()))
    pmax = SETTINGS["peptide_max_length"]
    prow, drow = [], []
    for ent, g in idx.groupby("entry_id"):
        r = rc[rc.pdb_id == ent]
        kinds, label_of, desc, length, ident = {}, {}, {}, {}, {}
        for x in r.itertuples():
            if x.identity in ("esr1", "esr2"):
                kinds[x.auth_asym_id] = "receptor"
                label_of[x.auth_asym_id] = x.label_asym_id
                ident[x.auth_asym_id] = x.identity
            elif str(x.polymer_type).startswith("polypeptide") and x.sequence_length <= pmax:
                kinds[x.auth_asym_id] = "peptide"
            else:
                kinds.setdefault(x.auth_asym_id, "other")
            desc[x.auth_asym_id] = x.description
            length[x.auth_asym_id] = x.sequence_length
        model, ns, frame = structure(ROOT / "data/raw" / f"{ent}.cif")
        for auth in sorted(a for a, k in kinds.items() if k == "receptor"):
            found, ok = contacts(model, ns, auth, None, kinds, {"peptide"})
            key, n, second, ncand = best(found)
            prow.append({"entry_id": ent, "auth_asym_id": auth, "label_asym_id": label_of[auth],
                         "chain_key": f"{ent}|{label_of[auth]}|1", "chain_in_model": ok,
                         "peptide_bound": key is not None, "peptide_auth_asym_id": key[0] if key else "",
                         "peptide_source": key[2] if key else "none", "peptide_contact_atoms": n,
                         "peptide_second_best_atoms": second, "n_peptide_copies_in_contact": ncand,
                         "peptide_length": length.get(key[0]) if key else np.nan,
                         "peptide_description": desc.get(key[0], "") if key else "", "crystal_frame": frame})
        auth_of = {v: k for k, v in label_of.items()}
        for c in g.itertuples():
            label = c.chain_key.split("|")[1]
            auth = auth_of.get(label)
            res = h11.get(c.chain_key, set())
            found, ok = contacts(model, ns, auth, res, kinds, {"receptor"}) if auth else ({}, False)
            key, n, second, ncand = best(found)
            drow.append({"chain_key": c.chain_key, "entry_id": ent, "auth_asym_id": auth, "n_h11_residues": len(res),
                         "dimer_partner_auth_asym_id": key[0] if key else "",
                         "dimer_partner_chain_key": f"{ent}|{label_of[key[0]]}|1" if key else "",
                         "dimer_partner_identity": ident.get(key[0], "") if key else "",
                         "dimer_partner_source": key[2] if key else "none",
                         "dimer_partner_image_idx": key[1] if key else np.nan,
                         "dimer_partner_self_image": bool(key and key[0] == auth),
                         "dimer_partner_h11_atoms": n, "dimer_second_best_h11_atoms": second,
                         "n_receptor_copies_at_h11": ncand, "crystal_frame": frame})
    pep = pd.DataFrame(prow)
    dim = pd.DataFrame(drow)
    pep.to_csv(OUT / "chain_peptide.csv", index=False)
    dim.to_csv(OUT / "dimer_partner.csv", index=False)
    print("dimer partner source:", dim.dimer_partner_source.value_counts().to_dict(),
          "self image:", int(dim.dimer_partner_self_image.sum()))
    kept_pep = pep[pep.chain_key.isin(idx.chain_key)]
    print("kept chains with a bound peptide:", int(kept_pep.peptide_bound.sum()), "of", len(kept_pep),
          kept_pep.peptide_source.value_counts().to_dict())


if __name__ == "__main__":
    main()
