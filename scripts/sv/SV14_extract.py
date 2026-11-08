"""SV14 stage 1: receptor-chain and ligand-of-interest heavy atoms in the core v3 frame, plus the apo check.

Inputs (rule 1.10):
  - SV00 card: results/sv/SV00/chain_table.csv (788 kept chains), data/sv/lbd_extract/receptor_chains.csv
    (label_asym_id, model) and residue_mapping.csv (label_seq_id -> UniProt position), restored mmCIFs data/raw/.
  - Q18 card: results/chain_ligand_of_interest.csv (ligand_key = entry|label_asym_id|auth_seq_id+ins per instance).
  - SV11 card: results/sv/SV11/site_residues.csv (ligand-pocket site, in_site True: 18 residues).
  - SV03 card: data/sv/SV03/core_v3_transforms.csv via scripts/sv/common/core_frame.py.
  - SV12 cache: data/sv/rcsb_chemcomp/*.json (rcsb_chem_comp_info.atom_count_heavy per CCD).

Atom selection ([DEFAULT] step 1): the receptor chain's own polymer atoms only (label_asym_id of the chain, model 1,
label_seq_id present), heavy atoms (type_symbol not H/D), altlocs by highest occupancy then blank/A via
esr_analysis.read_structure.read_cif (the SV00 code path). data/clean/*.pdb (what P2Rank saw) hold the same atom
selection but keep every altloc atom (no altloc column) and any hydrogens, so they are not used for volumes.
Ligand of interest: the chain's Q18 instances, heavy atoms, altloc-selected per atom; when instances carry
different alternate-location ids (alternate binding modes modelled as separate residues), the alt id with the
highest occupancy is kept (ties: alphabetical), instances with a blank alt id are always kept.

Apo check ([DEFAULT] step 5) for chains without a ligand of interest: any residue other than the chain's own polymer
residues and water (HOH/DOD/WAT), with >= 6 heavy atoms, having a heavy atom within 4.0 A of a heavy atom of the
chain's ligand-pocket residues (asymmetric unit, model 1, altloc-selected atoms).

Outputs: data/sv/SV14/atoms_core.npz, data/sv/SV14/chains.csv, data/sv/SV14/ligand_instances.csv,
data/sv/SV14/apo_check.csv.
Run: .venv/bin/python -B scripts/sv/SV14_extract.py
"""
from __future__ import annotations

import glob
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts/sv"))
from esr_analysis.read_structure import read_cif  # noqa: E402
from common.core_frame import load_transforms, to_core  # noqa: E402

OUT = ROOT / "data/sv/SV14"
MISSING = ("?", ".")
WATER = {"HOH", "DOD", "WAT"}
BONDI = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98}
OTHER_R = 1.80
APO_CUTOFF_A = 4.0
APO_MIN_HEAVY = 6


def radius(sym: str) -> float:
    return BONDI.get(str(sym).strip().upper(), OTHER_R)


def heavy(df: pd.DataFrame) -> pd.DataFrame:
    return df[~df.type_symbol.str.strip().str.upper().isin(["H", "D"])]


def site_positions() -> set[int]:
    s = pd.read_csv(ROOT / "results/sv/SV11/site_residues.csv")
    s = s[(s.site == "ligand pocket") & (s.in_site.astype(str) == "True")]
    return set(s.uniprot_position.astype(int))


def ccd_heavy() -> dict:
    out = {}
    for f in glob.glob(str(ROOT / "data/sv/rcsb_chemcomp/chemcomp_*.json")):
        d = json.loads(Path(f).read_text())
        if "chem_comp" in d:
            out[d["chem_comp"]["id"]] = int(d["rcsb_chem_comp_info"]["atom_count_heavy"])
    return out


def choose_instances(inst: list[dict]) -> tuple[list[dict], str]:
    """Blank-alt instances always; among alternate ids keep the one with the highest occupancy (ties alphabetical)."""
    blank = [i for i in inst if i["alt_ids"] in ("", ".")]
    alts = [i for i in inst if i["alt_ids"] not in ("", ".")]
    if not alts:
        return blank, "all instances (no alternate ids)"
    occ = {}
    for i in alts:
        occ.setdefault(i["alt_ids"], []).append(i["mean_occupancy"])
    best = sorted(occ, key=lambda a: (-max(occ[a]), a))[0]
    kept = blank + [i for i in alts if i["alt_ids"] == best]
    return kept, f"alt id {best} kept of {sorted(occ)} (highest occupancy)"


def process_entry(args):
    entry_id, chains, lig_rows, resmap, site, transforms = args
    _, atoms = read_cif(ROOT / "data/raw" / f"{entry_id}.cif")
    atoms = atoms[atoms.pdbx_PDB_model_num == "1"]
    atoms = heavy(atoms).copy()
    ins = atoms.pdbx_PDB_ins_code.where(~atoms.pdbx_PDB_ins_code.isin(MISSING), "")
    atoms["res_key"] = atoms.auth_seq_id.astype(str) + ins
    xyz_all = atoms[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float)
    out = []
    for ch in chains:
        ck = ch["chain_key"]
        own = (atoms.label_asym_id == ch["label_asym_id"]) & ~atoms.label_seq_id.isin(MISSING)
        rec = atoms[own]
        mp = resmap.get(ck, {})
        pos = rec.label_seq_id.astype(int).map(mp).fillna(-1).astype(int).to_numpy()
        rxyz = to_core(rec[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float), transforms[ck])
        rrad = rec.type_symbol.map(radius).to_numpy(float)
        # ligand of interest
        inst = []
        for lr in lig_rows.get(ck, []):
            _, asym, rk = lr["ligand_key"].split("|")
            la = atoms[(atoms.label_asym_id == asym) & (atoms.label_comp_id == lr["ligand_id"]) & (atoms.res_key == rk)]
            alt = sorted(set(la.label_alt_id) - set(MISSING))
            inst.append({"chain_key": ck, "ligand_id": lr["ligand_id"], "ligand_key": lr["ligand_key"],
                         "n_heavy_modelled": int(len(la)), "alt_ids": ";".join(alt),
                         "mean_occupancy": float(la.occupancy.mean()) if len(la) else np.nan,
                         "index": la.index.to_numpy()})
        kept, rule = choose_instances(inst) if inst else ([], "no ligand of interest")
        kept_keys = {i["ligand_key"] for i in kept}
        for i in inst:
            i["kept"] = i["ligand_key"] in kept_keys
            i["instance_rule"] = rule
        if kept:
            idx = np.concatenate([i["index"] for i in kept])
            la = atoms.loc[idx]
            lxyz = to_core(la[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float), transforms[ck])
            lrad = la.type_symbol.map(radius).to_numpy(float)
        else:
            lxyz, lrad = np.zeros((0, 3)), np.zeros(0)
        # apo check (chains without a ligand of interest)
        apo = None
        if not inst:
            site_mask = np.isin(pos, list(site))
            sxyz = rec[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float)[site_mask]
            other = ~own.to_numpy() & ~atoms.label_comp_id.isin(WATER).to_numpy()
            oat = atoms[other]
            oxyz = xyz_all[other]
            rid = oat.label_asym_id + "|" + oat.label_comp_id + "|" + oat.res_key + "|" + oat.label_seq_id
            nheavy = rid.map(rid.value_counts())
            near = np.zeros(len(oat), bool)
            if len(oat) and len(sxyz):
                dmin, _ = cKDTree(sxyz).query(oxyz, k=1, distance_upper_bound=APO_CUTOFF_A + 1e-9)
                near = dmin <= APO_CUTOFF_A
            found = sorted(set(rid[near & (nheavy.to_numpy() >= APO_MIN_HEAVY)]))
            near_small = sorted(set(rid[near & (nheavy.to_numpy() < APO_MIN_HEAVY)]))
            apo = {"chain_key": ck, "entry_id": entry_id, "n_site_heavy_atoms": int(len(sxyz)),
                   "n_site_residues_observed": int(len(set(pos[site_mask]))),
                   "residues_ge6_within_4A": ";".join(found), "n_residues_ge6_within_4A": len(found),
                   "residues_lt6_within_4A": ";".join(near_small), "passes_apo_check": len(found) == 0}
        out.append({"chain_key": ck, "rec_xyz": rxyz.astype(np.float32), "rec_rad": rrad.astype(np.float32),
                    "rec_pos": pos.astype(np.int16), "lig_xyz": lxyz.astype(np.float32),
                    "lig_rad": lrad.astype(np.float32),
                    "inst": [{k: v for k, v in i.items() if k != "index"} for i in inst], "apo": apo})
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv")
    rc = pd.read_csv(ROOT / "data/sv/lbd_extract/receptor_chains.csv", dtype=str)
    ct = ct.merge(rc[["chain_key", "label_asym_id", "model"]], on="chain_key", how="left")
    assert ct.label_asym_id.notna().all() and (ct.model == "1").all()
    rm = pd.read_csv(ROOT / "data/sv/lbd_extract/residue_mapping.csv",
                     usecols=["chain_key", "label_seq_id", "esr1_position"])
    rm = rm[rm.chain_key.isin(set(ct.chain_key))]
    resmap = {k: dict(zip(g.label_seq_id.astype(int), g.esr1_position.astype(int))) for k, g in rm.groupby("chain_key")}
    q18 = pd.read_csv(ROOT / "results/chain_ligand_of_interest.csv", dtype=str, keep_default_na=False)
    q18 = q18[q18.ligand_id != ""]
    lig_rows = {k: g.to_dict("records") for k, g in q18.groupby("chain_key")}
    site = site_positions()
    transforms = load_transforms()
    assert set(ct.chain_key) <= set(transforms)
    jobs = []
    for e, g in ct.groupby("entry_id", sort=True):
        keys = set(g.chain_key)
        jobs.append((e, g[["chain_key", "label_asym_id"]].to_dict("records"),
                     {k: v for k, v in lig_rows.items() if k in keys}, {k: resmap[k] for k in keys}, site,
                     {k: transforms[k] for k in keys}))
    with Pool(8) as pool:
        res = [r for chunk in pool.imap(process_entry, jobs, chunksize=2) for r in chunk]
    order = {k: i for i, k in enumerate(ct.chain_key)}
    res.sort(key=lambda r: order[r["chain_key"]])
    ccd = ccd_heavy()

    def cat(name):
        return np.concatenate([r[name] for r in res])

    rec_n = np.array([len(r["rec_rad"]) for r in res])
    lig_n = np.array([len(r["lig_rad"]) for r in res])
    np.savez_compressed(OUT / "atoms_core.npz", chain_keys=np.array([r["chain_key"] for r in res]),
                        rec_count=rec_n, rec_xyz=cat("rec_xyz"), rec_rad=cat("rec_rad"), rec_pos=cat("rec_pos"),
                        lig_count=lig_n, lig_xyz=cat("lig_xyz"), lig_rad=cat("lig_rad"))
    inst = pd.DataFrame([i for r in res for i in r["inst"]])
    inst["ccd_heavy_atoms"] = inst.ligand_id.map(ccd)
    inst["incomplete"] = inst.n_heavy_modelled < inst.ccd_heavy_atoms
    inst.to_csv(OUT / "ligand_instances.csv", index=False)
    apo = pd.DataFrame([r["apo"] for r in res if r["apo"] is not None])
    apo.to_csv(OUT / "apo_check.csv", index=False)
    kept = inst[inst.kept]
    agg = kept.groupby("chain_key").agg(
        ligand_ids_used=("ligand_id", lambda s: ";".join(sorted(set(s)))),
        n_instances_used=("ligand_key", "size"), n_instances_q18=("ligand_key", "size"),
        lig_heavy_modelled=("n_heavy_modelled", "sum"), lig_heavy_ccd=("ccd_heavy_atoms", "sum"),
        ligand_incomplete=("incomplete", "any"), instance_rule=("instance_rule", "first"))
    agg["n_instances_q18"] = inst.groupby("chain_key").size().reindex(agg.index)
    ch = pd.DataFrame({"chain_key": [r["chain_key"] for r in res], "n_rec_heavy": rec_n, "n_lig_heavy": lig_n,
                       "n_site_residues_observed": [len(set(r["rec_pos"][np.isin(r["rec_pos"], list(site))]))
                                                    for r in res]})
    ch = ch.merge(agg.reset_index(), on="chain_key", how="left")
    ch.to_csv(OUT / "chains.csv", index=False)
    print(f"chains {len(ch)}; receptor heavy atoms {rec_n.sum()}; ligand heavy atoms {lig_n.sum()}; "
          f"instances {len(inst)} (kept {int(inst.kept.sum())}); CCD missing {int(inst.ccd_heavy_atoms.isna().sum())}; "
          f"incomplete chains {int(ch.ligand_incomplete.fillna(False).sum())}; apo-check chains {len(apo)} "
          f"(pass {int(apo.passes_apo_check.sum()) if len(apo) else 0})")


if __name__ == "__main__":
    main()
