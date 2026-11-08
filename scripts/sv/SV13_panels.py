"""SV13 (RUN_SV1028) Part C data: one saved view, class representatives, C-alpha traces and ligand atoms in the core
v3 frame, projected onto the view (the figure is drawn from these CSVs by scripts/sv/SV13_figures.py).

PyMOL is not available (PRE addition, docs/sv/env_pre_1104.txt), so the step's fallback is used: C-alpha traces
projected with matplotlib onto the saved view (step 12). The representative's "cartoon" is drawn as a thick C-alpha
trace; its ligand of interest as sticks (bonds from interatomic distances).

Inputs (rule 1.10):
  - core v3 frame: data/sv/SV03/core_v3_transforms.csv via scripts/sv/common/core_frame.py (to_core) and the SV00 gate
    C-alpha in that frame (load_ca_core: positions 305-550, NaN = not observed; unmodelled residues stay absent);
    core v3 ranges esr_analysis/resources/config/alignment_core_v3.json.
  - classes: results/sv/SV12/chain_class.csv (primary_set rows; connectivity_only rows carry no class).
  - S04/S05 cards: results/rmsd_entry.npy + results/rmsd_entry_index.csv (existing chain x chain C-alpha RMSD).
  - SV00 chain table (resolution_A, s08 cluster); SV01a chain_outcomes.csv (H12 complete 538-547).
  - SV01b mobility_per_residue.csv (per-residue C-alpha B z vs the chain's core v3 C-alpha B; 528-547).
  - Q18 card: results/chain_ligand_of_interest.csv (ligand instances); restored mmCIFs data/raw/<entry>.cif read with
    esr_analysis.read_structure.read_cif (the SV00/SV14 code path); instance rule = SV14 choose_instances
    (cross-checked against data/sv/SV14/atoms_core.npz).

Rules (step 8-10; declared where the run file is open):
  - View: principal axes of the C-alpha of 5TM1|A|1 (cl7 medoid, results/medoids.csv) at the core v3 positions plus
    its H12 538-547, in the core v3 frame (the PyMOL 'orient' convention): screen x = PC1, y = PC2, z = x cross y
    (towards the viewer); signs chosen so that the 5TM1|A|1 H12 centroid lies at +x and +y. Orthographic
    projection [sx, sy, depth] = R (x_core - centre). Axis limits = the box of everything drawn in any panel plus
    2 A, identical for every panel. Saved in results/sv/SV13/view.txt.
  - Class chains = primary-set chains of the class (all clusters). Representative [DEFAULT]: lowest median RMSD
    (results/rmsd_entry.npy) to the other chains of the class; ties by resolution (lower first), then chain_key.
  - Traces: residues 528-550 of every class chain (observed C-alpha only); the representative's full observed
    305-550 trace for context; its H12 538-547 coloured by the SV01b per-residue B z.
  - Ligand bonds: heavy-atom pairs closer than r_cov(i) + r_cov(j) + 0.45 A (covalent radii C 0.76, N 0.71,
    O 0.66, F 0.57, S 1.05, Cl 1.02, Br 1.20, I 1.39, other 0.80).
Outputs: results/sv/SV13/view.txt, representatives.csv; data/sv/SV13/panel_traces.csv, panel_ligand_atoms.csv,
panel_ligand_bonds.csv, panel_meta.json.
Run: .venv/bin/python -B scripts/sv/SV13_panels.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts/sv"))
from esr_analysis.read_structure import read_cif  # noqa: E402
from common.core_frame import load_ca_core, load_transforms, to_core  # noqa: E402
from SV14_extract import choose_instances  # noqa: E402  (SV14 instance rule, imported not copied)

OUT = ROOT / "results/sv/SV13"
DOUT = ROOT / "data/sv/SV13"
VIEW_CHAIN = "5TM1|A|1"
H12 = (538, 547)
TRACE = (528, 550)
CLASSES = ("agonist", "partial_agonist", "SERM", "SERD", "antagonist", "other")
MISSING = ("?", ".")
COV = {"C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "S": 1.05, "CL": 1.02, "BR": 1.20, "I": 1.39}
COV_OTHER = 0.80
BOND_TOL = 0.45


def core_mask(pos):
    cfg = json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())
    m = np.zeros(len(pos), bool)
    for lo, hi in cfg["core_ranges"]:
        m |= (pos >= lo) & (pos <= hi)
    return m, cfg["core_ranges"]


def make_view(keys, pos, ca):
    i = list(keys).index(VIEW_CHAIN)
    cm, ranges = core_mask(pos)
    hm = (pos >= H12[0]) & (pos <= H12[1])
    sel = (cm | hm) & np.isfinite(ca[i]).all(1)
    X = ca[i][sel]
    centre = X.mean(0)
    _, _, vt = np.linalg.svd(X - centre, full_matrices=False)
    x, y = vt[0].copy(), vt[1].copy()
    h12c = ca[i][hm & np.isfinite(ca[i]).all(1)].mean(0) - centre
    if h12c @ x < 0:
        x = -x
    if h12c @ y < 0:
        y = -y
    z = np.cross(x, y)
    R = np.vstack([x, y, z])
    info = dict(view_chain=VIEW_CHAIN, n_atoms=int(sel.sum()), n_core=int((cm & sel).sum()),
                n_h12=int((hm & sel).sum()), core_ranges=ranges, centre=centre.tolist(), rotation=R.tolist(),
                h12_centroid_screen=(R @ h12c).tolist())
    return R, centre, info


def project(xyz, R, centre):
    return (np.asarray(xyz, float) - centre) @ R.T


def representatives(cc, ct, co):
    M = np.load(ROOT / "results/rmsd_entry.npy")
    idx = pd.read_csv(ROOT / "results/rmsd_entry_index.csv")
    row = pd.Series(idx.row.values, index=idx.chain_key)
    res = ct.set_index("chain_key").resolution_A
    comp = co.set_index("chain_key").complete
    rows = []
    for cls in CLASSES:
        s = cc[(cc.action_class == cls) & cc.primary_set]
        if not len(s):
            rows.append(dict(action_class=cls, n_chains=0, n_entries=0, n_complete_h12=0, representative=""))
            continue
        keys = s.chain_key.tolist()
        r = row.loc[keys].values
        sub = M[np.ix_(r, r)]
        med = np.array([np.median(np.delete(sub[k], k)) if len(keys) > 1 else 0.0 for k in range(len(keys))])
        cand = pd.DataFrame(dict(chain_key=keys, median_rmsd=med, resolution_A=res.loc[keys].values))
        cand = cand.sort_values(["median_rmsd", "resolution_A", "chain_key"]).reset_index(drop=True)
        best = cand.iloc[0]
        n_tied = int((cand.median_rmsd == best.median_rmsd).sum())
        b = s[s.chain_key == best.chain_key].iloc[0]
        rows.append(dict(action_class=cls, n_chains=len(s), n_entries=s.entry_id.nunique(),
                         n_complete_h12=int(comp.loc[keys].sum()), representative=best.chain_key,
                         rep_entry=b.entry_id, rep_ligand=b.ccd, rep_ligand_codes=b.ligand_codes,
                         rep_cluster=b.cluster, rep_s08_cluster=int(b.s08_cluster),
                         rep_median_rmsd_A=float(best.median_rmsd), rep_resolution_A=float(best.resolution_A),
                         rep_h12_complete=int(comp.loc[best.chain_key]),
                         n_tied_at_min=n_tied, second_chain=cand.chain_key[1] if len(cand) > 1 else "",
                         second_median_rmsd_A=float(cand.median_rmsd[1]) if len(cand) > 1 else np.nan,
                         rule="lowest median RMSD (results/rmsd_entry.npy) to the other chains of the class; "
                              "ties by resolution, then chain_key"))
    return pd.DataFrame(rows)


def ligand_atoms(chain_key, entry_id, transforms):
    lig = pd.read_csv(ROOT / "results/chain_ligand_of_interest.csv")
    lig = lig[lig.chain_key == chain_key]
    _, atoms = read_cif(ROOT / "data/raw" / f"{entry_id}.cif")
    atoms = atoms[atoms.pdbx_PDB_model_num == "1"]
    atoms = atoms[~atoms.type_symbol.str.strip().str.upper().isin(["H", "D"])].copy()
    ins = atoms.pdbx_PDB_ins_code.where(~atoms.pdbx_PDB_ins_code.isin(MISSING), "")
    atoms["res_key"] = atoms.auth_seq_id.astype(str) + ins
    inst = []
    for lr in lig.itertuples():
        _, asym, rk = lr.ligand_key.split("|")
        la = atoms[(atoms.label_asym_id == asym) & (atoms.label_comp_id == lr.ligand_id) & (atoms.res_key == rk)]
        alt = sorted(set(la.label_alt_id) - set(MISSING))
        inst.append({"ligand_id": lr.ligand_id, "ligand_key": lr.ligand_key, "alt_ids": ";".join(alt),
                     "mean_occupancy": float(la.occupancy.mean()) if len(la) else np.nan, "df": la})
    kept, rule = choose_instances(inst) if inst else ([], "no ligand of interest")
    if not kept:
        return pd.DataFrame(), rule
    la = pd.concat([k["df"] for k in kept])
    xyz = to_core(la[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float), transforms[chain_key])
    return pd.DataFrame(dict(chain_key=chain_key, ligand_id=la.label_comp_id.values, atom=la.label_atom_id.values,
                             element=la.type_symbol.str.strip().str.upper().values, x=xyz[:, 0], y=xyz[:, 1],
                             z=xyz[:, 2])), rule


def bonds(df):
    xyz = df[["x", "y", "z"]].to_numpy()
    r = df.element.map(lambda e: COV.get(e, COV_OTHER)).to_numpy()
    out = []
    for i in range(len(df)):
        d = np.linalg.norm(xyz[i + 1:] - xyz[i], axis=1)
        for j in np.where(d < r[i] + r[i + 1:] + BOND_TOL)[0]:
            out.append((i, i + 1 + j, float(d[j])))
    return out


def main():
    DOUT.mkdir(parents=True, exist_ok=True)
    keys, pos, ca = load_ca_core()
    kidx = {k: i for i, k in enumerate(keys)}
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", low_memory=False)
    cc = pd.read_csv(ROOT / "results/sv/SV12/chain_class.csv")
    cc["primary_set"] = cc.primary_set.astype(str).eq("True")
    co = pd.read_csv(ROOT / "results/sv/SV01a/chain_outcomes.csv")
    bz = pd.read_csv(ROOT / "results/sv/SV01b/mobility_per_residue.csv")
    transforms = load_transforms()
    R, centre, vinfo = make_view(keys, pos, ca)
    reps = representatives(cc, ct, co)
    reps.to_csv(OUT / "representatives.csv", index=False)

    tw = (pos >= TRACE[0]) & (pos <= TRACE[1])
    trows = []
    for cls in CLASSES:
        s = cc[(cc.action_class == cls) & cc.primary_set]
        rep = reps.set_index("action_class").representative.get(cls, "")
        for r in s.itertuples():
            i = kidx[r.chain_key]
            for role, m in (("class_trace", tw), ("rep_full", np.ones(len(pos), bool))):
                if role == "rep_full" and r.chain_key != rep:
                    continue
                for j in np.where(m)[0]:
                    p = ca[i, j]
                    trows.append(dict(action_class=cls, chain_key=r.chain_key, entry_id=r.entry_id,
                                      cluster=r.cluster, role=role, is_representative=r.chain_key == rep,
                                      position=int(pos[j]), observed=bool(np.isfinite(p).all()),
                                      x=p[0], y=p[1], z=p[2]))
    tr = pd.DataFrame(trows)
    pj = project(tr[["x", "y", "z"]].to_numpy(), R, centre)
    tr["sx"], tr["sy"], tr["depth"] = pj[:, 0], pj[:, 1], pj[:, 2]
    bzi = bz.set_index(["chain_key", "esr1_position"]).b_z
    tr["b_z"] = [bzi.get((k, p), np.nan) for k, p in zip(tr.chain_key, tr.position)]
    tr.to_csv(DOUT / "panel_traces.csv", index=False)

    # ligands of the representatives, cross-checked with the SV14 atom archive
    a14 = np.load(ROOT / "data/sv/SV14/atoms_core.npz")
    off = np.concatenate([[0], np.cumsum(a14["lig_count"])])
    k14 = {k: i for i, k in enumerate(a14["chain_keys"].astype(str))}
    latoms, lbonds, lcheck = [], [], {}
    for r in reps[reps.representative != ""].itertuples():
        la, rule = ligand_atoms(r.representative, r.rep_entry, transforms)
        i = k14[r.representative]
        ref = a14["lig_xyz"][off[i]:off[i + 1]].astype(float)
        got = la[["x", "y", "z"]].to_numpy()
        dmax = float(cKDTree(ref).query(got)[0].max()) if len(ref) and len(got) else np.nan
        lcheck[r.representative] = dict(n_atoms=len(got), n_atoms_sv14=int(len(ref)), max_abs_diff_A=dmax,
                                        instance_rule=rule, ligand_ids=sorted(set(la.ligand_id)) if len(la) else [])
        la = la.reset_index(drop=True)
        la.insert(0, "action_class", r.action_class)
        la["atom_index"] = np.arange(len(la))
        p = project(la[["x", "y", "z"]].to_numpy(), R, centre)
        la["sx"], la["sy"], la["depth"] = p[:, 0], p[:, 1], p[:, 2]
        latoms.append(la)
        for a, b, dd in bonds(la):
            lbonds.append(dict(action_class=r.action_class, chain_key=r.representative, i=a, j=b, dist_A=dd))
    la = pd.concat(latoms, ignore_index=True)
    lb = pd.DataFrame(lbonds)
    la.to_csv(DOUT / "panel_ligand_atoms.csv", index=False)
    lb.to_csv(DOUT / "panel_ligand_bonds.csv", index=False)

    obs = tr[tr.observed]
    xs = np.r_[obs.sx.values, la.sx.values]
    ys = np.r_[obs.sy.values, la.sy.values]
    lim = dict(xlim=[float(xs.min() - 2), float(xs.max() + 2)], ylim=[float(ys.min() - 2), float(ys.max() + 2)])
    rep_bz = obs[obs.is_representative & obs.position.between(*H12)].b_z
    meta = dict(view=vinfo, limits=lim, ligand_check=lcheck,
                b_z_range_rep_h12=[float(rep_bz.min()), float(rep_bz.max())],
                colour_range_b_z=[float(np.floor(rep_bz.min())), float(np.ceil(rep_bz.max()))],
                n_trace_rows=int(len(tr)), n_observed=int(tr.observed.sum()),
                bond_rule=f"d < r_cov(i) + r_cov(j) + {BOND_TOL} A", n_bonds=int(len(lb)))
    (DOUT / "panel_meta.json").write_text(json.dumps(meta, indent=1))

    Rm = np.asarray(vinfo["rotation"])
    lines = [
        "# SV13 view: one view for every Part C panel (RUN_SV1028 SV13 step 8)",
        "frame: core v3 (data/sv/SV03/core_v3_transforms.csv; x_core = (x_deposited - m) @ R_fit + c, "
        "scripts/sv/common/core_frame.py)",
        f"selection: {VIEW_CHAIN} (cl7 medoid) C-alpha at the core v3 positions {vinfo['core_ranges']} "
        f"({vinfo['n_core']} atoms) + H12 {H12[0]}-{H12[1]} ({vinfo['n_h12']} atoms); {vinfo['n_atoms']} atoms",
        "method: principal axes of the selection (PyMOL 'orient' convention): screen x = PC1, screen y = PC2, "
        "screen z = x cross y (towards the viewer); signs: the selection's H12 centroid at +x and +y",
        "projection: orthographic, [sx, sy, depth] = ROT @ (x_core - CENTRE); drawn with matplotlib (PyMOL not "
        "available: step 12 fallback)",
        "CENTRE (A, core v3 frame): " + " ".join(f"{v:.4f}" for v in centre),
        "ROT row 1 (screen x): " + " ".join(f"{v:.6f}" for v in Rm[0]),
        "ROT row 2 (screen y): " + " ".join(f"{v:.6f}" for v in Rm[1]),
        "ROT row 3 (screen z): " + " ".join(f"{v:.6f}" for v in Rm[2]),
        f"limits (A, screen): x {lim['xlim'][0]:.2f} {lim['xlim'][1]:.2f}; y {lim['ylim'][0]:.2f} "
        f"{lim['ylim'][1]:.2f} (box of everything drawn in any panel + 2 A); equal aspect",
        f"H12 centroid of {VIEW_CHAIN} on screen (A): " + " ".join(f"{v:.2f}" for v in vinfo["h12_centroid_screen"]),
        "used unchanged for every panel of figures/sv/SV13/SV13_h12_class_panels.png",
    ]
    (OUT / "view.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(reps.to_string())
    print(json.dumps(lcheck, indent=1))


if __name__ == "__main__":
    main()
