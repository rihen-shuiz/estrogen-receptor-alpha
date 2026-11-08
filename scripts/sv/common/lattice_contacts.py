"""Lattice contacts of chosen ESR1 residues with other chains (RUN_SV1028 SV01b step 6; SV07 reuses this module).

Definitions (fixed in SV01b before any model was fitted; the structure handling is the one of SV00_contacts.py / Y44):
- Structure: the restored mmCIF data/raw/{entry}.cif read with gemmi, model 1, waters, ligands and hydrogens removed,
  so "other chains" are polymer chains (receptor copies, peptides, other proteins). NeighborSearch with crystal
  symmetry (5 A cells); a non-standard SCALEn -> standardize_crystal_frame (a rigid change, SV00 rule).
- All conformers are kept for contacts (esr_analysis.read_structure.atom_records policy); a query atom is
  (residue, atom name) and is counted once whatever its altlocs.
- Query residues: the chain's residues mapped to the requested UniProt positions through the regenerated Q04
  residue mapping (data/sv/lbd_extract/residue_mapping.csv: auth_seq_id + insertion code of each esr1_position).
- Copies (SV00 rule): a neighbour atom found at symmetry image 0 whose deposited position is itself within the
  cut-off is an ASU copy; every other hit (another symmetry operation or a lattice translation, including the
  chain's own image) is a symmetry copy. Copies are keyed (auth chain, image index, 'asu' | 'symmetry').
- Contact: heavy-atom distance <= 4.0 A (CUTOFF_A).
- Chain kinds per auth chain, from data/sv/lbd_extract/receptor_chains.csv (model 1): receptor (identity esr1 or
  esr2), peptide (any other polypeptide of <= peptide_max_length residues, settings = 50), other.
- The chain's own bound peptide copy (SV00 definition: the peptide copy with the most heavy atoms of the receptor
  chain within 4.0 A; ties -> ASU before symmetry, then chain name, then image) is reported separately as
  own_peptide_* and is NOT a lattice contact: it is the 'peptide' predictor of SV01 block B.
- Lattice contact counts per residue = number of distinct query heavy atoms with >= 1 neighbour of that class
  within the cut-off: asu_atoms (ASU copies of other chains), sym_atoms (symmetry copies), each split by the
  neighbour's kind (receptor / peptide / other).

Main entry point: entry_contacts(entry_id, chains, positions) -> (per-residue DataFrame, per-chain DataFrame).
lattice_hits(...) returns the lattice neighbours of an arbitrary point (used by the SV01b medoid-H12 clash test).
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data/raw"
EXTRACT = ROOT / "data/sv/lbd_extract"
CUTOFF_A = 4.0
PEPTIDE_MAX = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())[
    "peptide_max_length"]
BACKBONE = ("N", "CA", "C", "O")


# ------------------------------------------------------------------------------------------------ inputs
@lru_cache(maxsize=1)
def _receptor_chains() -> pd.DataFrame:
    rc = pd.read_csv(EXTRACT / "receptor_chains.csv", keep_default_na=False, na_values=[""], low_memory=False)
    return rc[rc.model.astype(str) == "1"]


def chain_kinds(entry_id: str) -> dict:
    """{auth_asym_id: 'receptor' | 'peptide' | 'other'} for one entry (SV00_contacts rule)."""
    kinds = {}
    for x in _receptor_chains()[lambda r: r.pdb_id == entry_id].itertuples():
        if x.identity in ("esr1", "esr2"):
            kinds[x.auth_asym_id] = "receptor"
        elif str(x.polymer_type).startswith("polypeptide") and x.sequence_length <= PEPTIDE_MAX:
            kinds[x.auth_asym_id] = "peptide"
        else:
            kinds.setdefault(x.auth_asym_id, "other")
    return kinds


def residue_lookup(chain_keys) -> dict:
    """{chain_key: {(auth_seq_id, icode): esr1_position}} for residues with coordinates (icode ' ' when blank)."""
    keys = set(chain_keys)
    mp = pd.read_csv(EXTRACT / "residue_mapping.csv", low_memory=False,
                     usecols=["chain_key", "esr1_position", "has_residue", "auth_seq_id", "insertion_code"])
    mp = mp[mp.chain_key.isin(keys) & mp.has_residue.astype(bool) & mp.esr1_position.notna()]
    out = {k: {} for k in keys}
    for r in mp.itertuples():
        ic = r.insertion_code if isinstance(r.insertion_code, str) and r.insertion_code.strip() not in ("", "?", ".") \
            else " "
        out[r.chain_key][(int(float(r.auth_seq_id)), ic)] = int(r.esr1_position)
    return out


def load_structure(entry_id: str):
    """(model, NeighborSearch, frame) for data/raw/{entry}.cif; frame 'standard' or 'standardized'."""
    st = gemmi.read_structure(str(RAW / f"{entry_id}.cif"))
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


def get_chain(model, auth: str):
    return next((ch for ch in model if ch.name == auth), None)


# ------------------------------------------------------------------------------------------------ copies
def _copy_key(m, cra, at_pos, cutoff):
    asu = m.image_idx == 0 and cra.atom.pos.dist(at_pos) <= cutoff + 1e-6
    return cra.chain.name, int(m.image_idx), "asu" if asu else "symmetry"


def bound_peptide_copy(model, ns, auth: str, kinds: dict, cutoff: float = CUTOFF_A):
    """SV00 bound peptide of a receptor chain: the peptide copy with the most receptor heavy atoms within the
    cut-off (ties: ASU first, chain name, image). Returns ((chain, image, source), n_atoms) or (None, 0)."""
    chain = get_chain(model, auth)
    if chain is None:
        return None, 0
    found = {}
    for res in chain:
        for at in res:
            for m in ns.find_atoms(at.pos, "\0", radius=cutoff):
                cra = m.to_cra(model)
                if kinds.get(cra.chain.name) != "peptide":
                    continue
                key = _copy_key(m, cra, at.pos, cutoff)
                if key[2] == "asu" and key[0] == auth:
                    continue
                found.setdefault(key, set()).add((res.seqid.num, res.seqid.icode, at.name, at.altloc))
    if not found:
        return None, 0
    ranked = sorted(found.items(), key=lambda kv: (-len(kv[1]), kv[0][2] != "asu", kv[0][0], kv[0][1]))
    return ranked[0][0], len(ranked[0][1])


# ------------------------------------------------------------------------------------------------ contacts
def residue_contacts(model, ns, auth: str, lookup: dict, positions, kinds: dict, own_peptide=None,
                     cutoff: float = CUTOFF_A) -> list[dict]:
    """One row per query residue (UniProt position in `positions`) with lattice-contact atom counts."""
    chain = get_chain(model, auth)
    pos_set = set(int(p) for p in positions)
    rows = []
    if chain is None:
        return rows
    for res in chain:
        p = lookup.get((res.seqid.num, res.seqid.icode))
        if p is None or p not in pos_set:
            continue
        atoms = {}
        for at in res:
            atoms.setdefault(at.name, []).append(at)
        cls = {k: set() for k in ("asu", "sym", "asu_receptor", "asu_peptide", "asu_other", "sym_receptor",
                                  "sym_peptide", "sym_other", "own_peptide")}
        partners = {"asu": set(), "sym": set()}
        for name, confs in atoms.items():
            for at in confs:
                for m in ns.find_atoms(at.pos, "\0", radius=cutoff):
                    cra = m.to_cra(model)
                    key = _copy_key(m, cra, at.pos, cutoff)
                    if key[2] == "asu" and key[0] == auth:
                        continue  # the chain itself
                    if own_peptide is not None and key == own_peptide:
                        cls["own_peptide"].add(name)
                        continue
                    src = "asu" if key[2] == "asu" else "sym"
                    kind = kinds.get(key[0], "other")
                    cls[src].add(name)
                    cls[f"{src}_{kind}"].add(name)
                    partners[src].add(key[0] if src == "asu" else f"{key[0]}@{key[1]}")
        row = dict(esr1_position=p, auth_seq_id=res.seqid.num, icode=res.seqid.icode.strip(), residue=res.name,
                   n_heavy_atoms=len(atoms))
        row.update({f"{k}_atoms": len(v) for k, v in cls.items()})
        row.update(asu_partners=";".join(sorted(partners["asu"])), sym_partners=";".join(sorted(partners["sym"])))
        rows.append(row)
    return rows


def entry_contacts(entry_id: str, chains, positions, lookup: dict | None = None, cutoff: float = CUTOFF_A):
    """Lattice contacts for the chains [(chain_key, auth_asym_id), ...] of one entry at the UniProt `positions`.
    Returns (per-residue DataFrame, per-chain DataFrame with the own bound peptide copy and the crystal frame)."""
    chains = list(chains)
    lookup = lookup if lookup is not None else residue_lookup([c for c, _ in chains])
    kinds = chain_kinds(entry_id)
    model, ns, frame = load_structure(entry_id)
    res_rows, chain_rows = [], []
    for ck, auth in chains:
        pep, n_pep = bound_peptide_copy(model, ns, auth, kinds, cutoff)
        rows = residue_contacts(model, ns, auth, lookup.get(ck, {}), positions, kinds, own_peptide=pep,
                                cutoff=cutoff)
        for r in rows:
            r.update(chain_key=ck, entry_id=entry_id)
        res_rows += rows
        chain_rows.append(dict(chain_key=ck, entry_id=entry_id, auth_asym_id=auth, chain_in_model=get_chain(
            model, auth) is not None, own_peptide_auth_asym_id=pep[0] if pep else "",
            own_peptide_image_idx=pep[1] if pep else np.nan, own_peptide_source=pep[2] if pep else "none",
            own_peptide_contact_atoms=n_pep, crystal_frame=frame))
    cols = ["chain_key", "entry_id"]
    rdf = pd.DataFrame(res_rows)
    if len(rdf):
        rdf = rdf[cols + [c for c in rdf.columns if c not in cols]]
    return rdf, pd.DataFrame(chain_rows)


def summarise(per_residue: pd.DataFrame, lo: int, hi: int, prefix: str) -> pd.DataFrame:
    """Per-chain sums of the per-residue atom counts over UniProt lo-hi (columns prefixed)."""
    d = per_residue[per_residue.esr1_position.between(lo, hi)]
    cols = [c for c in d.columns if c.endswith("_atoms") and c != "n_heavy_atoms"]
    s = d.groupby("chain_key")[cols + ["n_heavy_atoms"]].sum()
    s["n_residues"] = d.groupby("chain_key").size()
    s["n_residues_asu_contact"] = d[d.asu_atoms > 0].groupby("chain_key").size()
    s["n_residues_sym_contact"] = d[d.sym_atoms > 0].groupby("chain_key").size()
    s = s.fillna(0).astype(int)
    s.columns = [f"{prefix}{c}" for c in s.columns]
    return s.reset_index()


# ------------------------------------------------------------------------------------------------ point queries
def lattice_hits(model, ns, point: gemmi.Position, self_auth: str, kinds: dict, exclude=None,
                 cutoff: float = CUTOFF_A) -> list[tuple]:
    """Lattice neighbours of an arbitrary point (e.g. a superposed medoid atom): [(auth, image, source, kind,
    distance)], the chain itself (ASU copy) and the copy `exclude` (its own bound peptide) left out."""
    out = []
    for m in ns.find_atoms(point, "\0", radius=cutoff):
        cra = m.to_cra(model)
        key = _copy_key(m, cra, point, cutoff)
        if key[2] == "asu" and key[0] == self_auth:
            continue
        if exclude is not None and key == exclude:
            continue
        # distance to the nearest lattice image (find_atoms selects by it; m.pos is the stored image)
        out.append(key + (kinds.get(key[0], "other"), float(ns.dist(point, m.pos))))
    return out


def selected_atoms(residue) -> dict:
    """{atom name: gemmi.Atom}: one conformer per atom, highest occupancy then blank/A (read_structure policy)."""
    best = {}
    for at in residue:
        rank = (-at.occ, 0 if at.altloc in ("\0", " ", "") else 1 if at.altloc == "A" else 2)
        if at.name not in best or rank < best[at.name][0]:
            best[at.name] = (rank, at)
    return {k: v[1] for k, v in best.items()}
