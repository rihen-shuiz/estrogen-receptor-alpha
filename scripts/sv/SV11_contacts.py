"""SV11 Part A, step 2 (structure part): receptor residues in contact with the bound coactivator peptide and with the
SV00 dimer partner, for every kept chain, computed on the SV00 structures.

Structures: the restored mmCIFs in data/raw (SV00), read exactly as scripts/sv/SV00_contacts.py does (its structure()
is imported): waters, ligands and hydrogens removed, model 1, gemmi NeighborSearch with crystal symmetry (5 A cells;
non-standard SCALEn -> standardize_crystal_frame). All conformers are kept, as in SV00.
Contact [DEFAULT, SV11 rule 2]: any receptor heavy atom within 4.0 A of any heavy atom of the partner copy.
Copies are grouped as in SV00: (partner chain, symmetry image index, asu/symmetry); an atom found at image 0 whose
deposited position is itself within the cut-off is an ASU copy, every other hit is a symmetry copy.
Chain roles from the SV00 regenerated Q04 table (data/sv/lbd_extract/receptor_chains.csv, model 1):
  receptor = identity esr1/esr2; peptide [SV11 rule 2] = a non-receptor polymer chain of <= 30 residues
  (entity sequence_length); everything else = other.
AF-2 coactivator groove: the peptide copy with the most receptor heavy atoms within 4.0 A of the chain (SV00 tie rule:
  ASU before symmetry, then chain name, then image index); chains with no peptide atom within 4.0 A have none.
Dimer interface: the SV00 dimer partner copy of the chain (results/sv/SV00/dimer_partner.csv: partner chain, image
  index, source); every residue of the chain in contact with that copy (not only H11).
Numbering: UniProt = author number + uniprot_author_offset (SV00 chain_table / chain_offsets.csv).
Checks written to the output: the bound peptide equals SV00 peptide_bound (chain, source, contact atoms) and the H11
  (UniProt 497-530) atom count of the dimer copy equals SV00 dimer_partner_h11_atoms.
Output: results/sv/SV11/chain_site_contacts.csv (two rows per kept chain: AF-2 groove, dimer interface).
Run: .venv/bin/python -B -I scripts/sv/SV11_contacts.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from SV00_contacts import structure  # noqa: E402  (identical structure handling to SV00)

OUT = ROOT / "results/sv/SV11"
CUTOFF_A = 4.0
PEPTIDE_MAX_RESIDUES = 30
H11 = (497, 530)
GROOVE, DIMER = "AF-2 groove", "dimer interface"


def all_contacts(model, ns, chain_name, kinds):
    """{(partner chain, image_idx, asu/symmetry): {(resnum, icode, atom, altloc)}} over every atom of the chain,
    partners of kind peptide or receptor only (the chain itself in the ASU is skipped)."""
    found = {}
    chain = next((ch for ch in model if ch.name == chain_name), None)
    if chain is None:
        return found, False
    for res in chain:
        for at in res:
            for m in ns.find_atoms(at.pos, "\0", radius=CUTOFF_A):
                cra = m.to_cra(model)
                if kinds.get(cra.chain.name) not in ("peptide", "receptor"):
                    continue
                asu = m.image_idx == 0 and cra.atom.pos.dist(at.pos) <= CUTOFF_A + 1e-6
                if asu and cra.chain.name == chain_name:
                    continue
                key = (cra.chain.name, int(m.image_idx), "asu" if asu else "symmetry")
                found.setdefault(key, set()).add((res.seqid.num, res.seqid.icode, at.name, at.altloc))
    return found, True


def ranked(found, keys):
    """SV00 tie rule: most atoms, ASU before symmetry, chain name, image index."""
    return sorted(keys, key=lambda k: (-len(found[k]), k[2] != "asu", k[0], k[1]))


def residues_uniprot(atoms, offset):
    res = sorted({(num, ic) for num, ic, _, _ in atoms})
    n_icode = sum(1 for _, ic in res if ic.strip())
    return sorted({int(num + offset) for num, _ in res}), n_icode


def main():
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", keep_default_na=False, na_values=[""])
    dp = pd.read_csv(ROOT / "results/sv/SV00/dimer_partner.csv", keep_default_na=False, na_values=[""])
    cp = pd.read_csv(ROOT / "results/sv/SV00/chain_peptide.csv", keep_default_na=False, na_values=[""])
    rc = pd.read_csv(ROOT / "data/sv/lbd_extract/receptor_chains.csv", keep_default_na=False, na_values=[""],
                     low_memory=False)
    rc = rc[rc.model.astype(str) == "1"]
    dp = dp.set_index("chain_key")
    cp = cp.set_index("chain_key")
    rows = []
    for ent, g in ct.groupby("entry_id", sort=True):
        kinds, length, ptype = {}, {}, {}
        for x in rc[rc.pdb_id == ent].itertuples():
            if x.identity in ("esr1", "esr2"):
                kinds[x.auth_asym_id] = "receptor"
            elif x.sequence_length <= PEPTIDE_MAX_RESIDUES:
                kinds[x.auth_asym_id] = "peptide"
            else:
                kinds.setdefault(x.auth_asym_id, "other")
            length[x.auth_asym_id] = x.sequence_length
            ptype[x.auth_asym_id] = x.polymer_type
        model, ns, frame = structure(ROOT / "data/raw" / f"{ent}.cif")
        for c in g.itertuples():
            offset = float(c.uniprot_author_offset)
            found, ok = all_contacts(model, ns, c.auth_asym_id, kinds)
            base = dict(chain_key=c.chain_key, entry_id=ent, auth_asym_id=c.auth_asym_id, chain_in_model=ok,
                        uniprot_author_offset=offset, crystal_frame=frame)
            # AF-2 groove: the bound peptide copy
            pkeys = ranked(found, [k for k in found if kinds.get(k[0]) == "peptide"])
            sv0 = cp.loc[c.chain_key]
            if pkeys:
                k = pkeys[0]
                res, nic = residues_uniprot(found[k], offset)
                row = dict(base, site=GROOVE, has_partner=True, partner_auth_asym_id=k[0], partner_image_idx=k[1],
                           partner_source=k[2], partner_length=length.get(k[0]), partner_polymer_type=ptype.get(k[0]),
                           n_contact_atoms=len(found[k]),
                           n_second_best_atoms=len(found[pkeys[1]]) if len(pkeys) > 1 else 0,
                           n_partner_copies_in_contact=len(pkeys), n_contact_residues=len(res),
                           n_contact_residues_icode=nic, contact_residues_uniprot=";".join(map(str, res)))
                row["sv00_match"] = bool(sv0.peptide_bound and sv0.peptide_auth_asym_id == k[0]
                                         and sv0.peptide_source == k[2]
                                         and int(sv0.peptide_contact_atoms) == len(found[k]))
            else:
                row = dict(base, site=GROOVE, has_partner=False, partner_auth_asym_id="", partner_image_idx=np.nan,
                           partner_source="none", partner_length=np.nan, partner_polymer_type="", n_contact_atoms=0,
                           n_second_best_atoms=0, n_partner_copies_in_contact=0, n_contact_residues=0,
                           n_contact_residues_icode=0, contact_residues_uniprot="",
                           sv00_match=not bool(sv0.peptide_bound))
            rows.append(row)
            # dimer interface: the SV00 partner copy
            d = dp.loc[c.chain_key]
            if d.dimer_partner_source in ("asu", "symmetry"):
                k = (d.dimer_partner_auth_asym_id, int(d.dimer_partner_image_idx), d.dimer_partner_source)
                atoms = found.get(k, set())
                res, nic = residues_uniprot(atoms, offset)
                h11_atoms = sum(1 for num, _, _, _ in atoms if H11[0] <= num + offset <= H11[1])
                rows.append(dict(base, site=DIMER, has_partner=True, partner_auth_asym_id=k[0], partner_image_idx=k[1],
                                 partner_source=k[2], partner_length=length.get(k[0]),
                                 partner_polymer_type=ptype.get(k[0]), n_contact_atoms=len(atoms),
                                 n_second_best_atoms=np.nan, n_partner_copies_in_contact=np.nan,
                                 n_contact_residues=len(res), n_contact_residues_icode=nic,
                                 contact_residues_uniprot=";".join(map(str, res)),
                                 partner_self_image=bool(d.dimer_partner_self_image), partner_h11_atoms=h11_atoms,
                                 sv00_match=h11_atoms == int(d.dimer_partner_h11_atoms)))
            else:
                rows.append(dict(base, site=DIMER, has_partner=False, partner_auth_asym_id="",
                                 partner_image_idx=np.nan, partner_source="none", n_contact_atoms=0,
                                 n_contact_residues=0, n_contact_residues_icode=0, contact_residues_uniprot="",
                                 sv00_match=True))
    out = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT / "chain_site_contacts.csv", index=False)
    for site, s in out.groupby("site"):
        print(site, "chains with partner:", int(s.has_partner.sum()), "of", len(s), "sv00 match:",
              int(s.sv00_match.sum()), "icode residues:", int(s.n_contact_residues_icode.sum()),
              "sources:", s.partner_source.value_counts().to_dict())


if __name__ == "__main__":
    main()
