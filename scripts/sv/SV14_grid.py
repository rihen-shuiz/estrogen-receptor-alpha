"""SV14 stage 2: grid pocket volume, ligand volume and filling per chain (RUN_SV1028 SV14 steps 1-4).

Input: data/sv/SV14/atoms_core.npz, chains.csv (SV14_extract.py; core v3 frame), results/sv/SV00/chain_table.csv,
results/sv/SV11/site_residues.csv (ligand-pocket site, in_site True).

[DEFAULT] as written in the run file:
  - grid spacing 0.4 A over the box of the ligand-pocket residues' heavy atoms + 6 A (grid lattice anchored at
    multiples of 0.4 A in the core v3 frame, so every chain uses the same lattice);
  - Bondi radii C 1.70, N 1.55, O 1.52, S 1.80, F 1.47, Cl 1.75, Br 1.85, I 1.98, other 1.80 A, no hydrogens;
  - probe centres = grid points >= r_vdW + 1.4 A from every receptor atom; cavity = grid points within 1.4 A of a
    probe centre (exact Euclidean distance transform, scipy.ndimage.distance_transform_edt);
  - buried = >= 10 of 14 rays (+-x, +-y, +-z, 8 cube diagonals) hit a receptor vdW sphere within 10 A; rays are
    marched on the grid (axes: 25 steps of 0.4 A = 10.0 A; diagonals: 14 steps of 0.693 A = 9.70 A, the last grid
    point within 10 A); the occupancy grid extends 26 points (10.4 A) beyond the box so rays see atoms outside it;
  - pocket = 6-connected component (scipy.ndimage.label) of buried cavity points overlapping the ligand of
    interest (vdW points of its heavy atoms); several overlapping components -> the one with the largest overlap
    (ties: larger component). Chains without a ligand: the component holding the buried cavity point nearest the
    E2 reference centroid (median, per coordinate, of the EST heavy-atom centroids over wild-type EST chains with a
    complete ligand);
  - volume = points x 0.064 A^3; pocket_open = the component is 6-adjacent to cavity points that are not buried
    (count recorded); sensitivity thresholds >= 8 and >= 12 rays;
  - ligand volume = vdW points of the ligand's modelled heavy atoms on the same lattice (counted on the extended
    grid, so ligand parts outside the box are included and count as outside the pocket);
  - filling % = 100 x ligand points inside the pocket component / pocket points; inside share = ligand points
    inside / all ligand points.
Output: data/sv/SV14/grid_per_chain.csv, data/sv/SV14/e2_centroid.json
Run: .venv/bin/python -B scripts/sv/SV14_grid.py
"""
from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/sv/SV14"
H = 0.4
VOXEL = H ** 3                     # 0.064 A^3
PROBE = 1.4
BOX_MARGIN = 6.0
RAY_LEN = 10.0
THRESHOLDS = (8, 10, 12)
PRIMARY_T = 10
PAD = 26                           # occupancy grid margin in points (10.4 A >= 10 A axis rays)
EXT = 4                            # probe-centre margin in points (1.6 A >= 1.4 A)
AXES = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
DIAGS = [(i, j, k) for i in (1, -1) for j in (1, -1) for k in (1, -1)]
DIRS = [(d, int(np.floor(RAY_LEN / H + 1e-9))) for d in AXES] + \
       [(d, int(np.floor(RAY_LEN / (H * np.sqrt(3)) + 1e-9))) for d in DIAGS]
CHUNK = 192


def paint(shape, origin, xyz, r, strict=False):
    """Boolean grid: True where some atom has |p - x| <= r (strict: < r). Exact distances, atoms in chunks."""
    g = np.zeros(shape, bool)
    if not len(r):
        return g
    shape_a = np.array(shape)
    far = origin + (shape_a - 1) * H
    keep = np.all((xyz + r[:, None] >= origin) & (xyz - r[:, None] <= far), axis=1)
    xyz, r = xyz[keep].astype(float), r[keep].astype(float)
    if not len(r):
        return g
    R = int(np.ceil(r.max() / H)) + 1
    off = np.arange(-R, R + 1)
    base = np.round((xyz - origin) / H).astype(int)
    flat = g.ravel()
    nx, ny, nz = shape
    for c in range(0, len(r), CHUNK):
        b, x, rr = base[c:c + CHUNK], xyz[c:c + CHUNK], r[c:c + CHUNK]
        ii = [b[:, a:a + 1] + off for a in range(3)]                      # (m, S) per axis
        d2a = []
        for a in range(3):
            d = (origin[a] + ii[a] * H - x[:, a:a + 1]) ** 2
            d[(ii[a] < 0) | (ii[a] >= shape_a[a])] = np.inf
            d2a.append(d)
        d2 = d2a[0][:, :, None, None] + d2a[1][:, None, :, None] + d2a[2][:, None, None, :]
        lim = (rr ** 2)[:, None, None, None]
        m = d2 < lim - 1e-9 if strict else d2 <= lim + 1e-9
        idx = (ii[0][:, :, None, None] * ny + ii[1][:, None, :, None]) * nz + ii[2][:, None, None, :]
        flat[idx[m]] = True
    return g


def chain_grids(rec_xyz, rec_rad, rec_pos, lig_xyz, lig_rad, site):
    """Box, receptor occupancy, cavity, ray hits and ligand points of one chain (steps 1-3)."""
    sm = np.isin(rec_pos, site)
    lo = np.floor((rec_xyz[sm].min(0) - BOX_MARGIN) / H) * H
    hi = np.ceil((rec_xyz[sm].max(0) + BOX_MARGIN) / H) * H
    n = np.round((hi - lo) / H).astype(int) + 1
    # receptor vdW occupancy on the padded grid (rays)
    o_pad = lo - PAD * H
    n_pad = tuple(n + 2 * PAD)
    occ = paint(n_pad, o_pad, rec_xyz, rec_rad)
    # probe centres on box + EXT, cavity = within 1.4 A of a probe centre
    o_ext = lo - EXT * H
    n_ext = tuple(n + 2 * EXT)
    excl = paint(n_ext, o_ext, rec_xyz, rec_rad + PROBE, strict=True)
    probe = ~excl
    if probe.any():
        dist = ndimage.distance_transform_edt(~probe) * H
        cav_ext = dist <= PROBE + 1e-9
    else:
        cav_ext = np.zeros(n_ext, bool)
    cav = cav_ext[EXT:-EXT, EXT:-EXT, EXT:-EXT]
    # buriedness: rays marched on the occupancy grid
    hits = np.zeros(tuple(n), np.int8)
    nx, ny, nz = n
    for (dx, dy, dz), K in DIRS:
        acc = np.zeros(tuple(n), bool)
        for k in range(1, K + 1):
            sx, sy, sz = PAD + k * dx, PAD + k * dy, PAD + k * dz
            acc |= occ[sx:sx + nx, sy:sy + ny, sz:sz + nz]
        hits += acc
    # ligand vdW points (extended grid: whole ligand; inner part for overlap)
    lig_pad = paint(n_pad, o_pad, lig_xyz, lig_rad) if len(lig_rad) else np.zeros(n_pad, bool)
    lig = lig_pad[PAD:-PAD, PAD:-PAD, PAD:-PAD]
    return {"n_site_atoms": int(sm.sum()), "lo": lo, "n": n, "occ": occ, "cav": cav, "hits": hits,
            "lig_pad": lig_pad, "lig": lig}


def select_component(lab, ncomp, lig, cand, has_ligand, lo, e2c):
    """Pocket label (0 = none), rule text and extra (overlapping components, or distance to the E2 centroid)."""
    L, how, extra = 0, "", np.nan
    if not ncomp:
        return L, how, extra
    sizes = np.bincount(lab.ravel(), minlength=ncomp + 1)
    if has_ligand:
        ov = np.bincount(lab[lig & cand], minlength=ncomp + 1)
        ov[0] = 0
        if ov.max() > 0:
            best = np.flatnonzero(ov == ov.max())
            L = int(best[np.argmax(sizes[best])])
            how = "ligand overlap"
            extra = int((ov > 0).sum())
        else:
            how = "no component overlaps the ligand"
    else:
        pts = np.argwhere(cand)
        d = np.linalg.norm(lo + pts * H - e2c, axis=1)
        j = int(np.argmin(d))
        L = int(lab[tuple(pts[j])])
        how = "nearest buried cavity point to E2 centroid"
        extra = float(d[j])
    return L, how, extra


def pocket_chain(args):
    key, rec_xyz, rec_rad, rec_pos, lig_xyz, lig_rad, site, e2c = args
    t0 = time.time()
    gr = chain_grids(rec_xyz, rec_rad, rec_pos, lig_xyz, lig_rad, site)
    lo, n, occ, cav, hits, lig_pad, lig = (gr[k] for k in ("lo", "n", "occ", "cav", "hits", "lig_pad", "lig"))
    out = {"chain_key": key, "n_site_heavy_atoms": gr["n_site_atoms"], "box_nx": int(n[0]), "box_ny": int(n[1]),
           "box_nz": int(n[2])}
    n_lig = int(lig_pad.sum())
    out.update(n_cavity_points=int(cav.sum()), ligand_points=n_lig, ligand_volume_A3=n_lig * VOXEL,
               ligand_points_outside_box=int(n_lig - lig.sum()),
               ligand_points_in_receptor_vdw=int((lig_pad & occ).sum()),
               ligand_points_in_cavity=int((lig & cav).sum()))
    if len(lig_rad):
        lc = lig_xyz.mean(0)
        out.update(lig_cx=lc[0], lig_cy=lc[1], lig_cz=lc[2], lig_centroid_to_e2_A=float(np.linalg.norm(lc - e2c)))
    six = ndimage.generate_binary_structure(3, 1)
    for t in THRESHOLDS:
        buried = hits >= t
        cand = cav & buried
        lab, ncomp = ndimage.label(cand, structure=six)
        sfx = f"_b{t}"
        out["n_components" + sfx] = int(ncomp)
        out["n_buried_cavity_points" + sfx] = int(cand.sum())
        L, how, extra = select_component(lab, ncomp, lig, cand, bool(len(lig_rad)), lo, e2c)
        out["pocket_rule" + sfx] = how
        if t == PRIMARY_T:
            out["n_overlapping_components"] = extra if len(lig_rad) else np.nan
            out["e2_centroid_to_nearest_point_A"] = extra if not len(lig_rad) else np.nan
        if L == 0:
            out.update({"pocket_points" + sfx: np.nan, "pocket_volume_A3" + sfx: np.nan,
                        "border_open_points" + sfx: np.nan, "ligand_points_inside" + sfx: np.nan,
                        "pocket_touches_box" + sfx: np.nan})
            continue
        comp = lab == L
        npts = int(comp.sum())
        border = ndimage.binary_dilation(comp, structure=six) & ~comp & cav & ~buried
        touches = bool(comp[0].any() or comp[-1].any() or comp[:, 0].any() or comp[:, -1].any()
                       or comp[:, :, 0].any() or comp[:, :, -1].any())
        inside = int((lig & comp).sum())
        out.update({"pocket_points" + sfx: npts, "pocket_volume_A3" + sfx: npts * VOXEL,
                    "border_open_points" + sfx: int(border.sum()), "ligand_points_inside" + sfx: inside,
                    "pocket_touches_box" + sfx: touches})
        if t == PRIMARY_T:
            pc = (lo + np.argwhere(comp) * H).mean(0)
            out.update(pocket_cx=pc[0], pocket_cy=pc[1], pocket_cz=pc[2],
                       pocket_centroid_to_e2_A=float(np.linalg.norm(pc - e2c)))
    out["seconds"] = round(time.time() - t0, 2)
    return out


def main():
    a = np.load(DATA / "atoms_core.npz")
    keys = a["chain_keys"].astype(str)
    rs = np.concatenate([[0], np.cumsum(a["rec_count"])])
    ls = np.concatenate([[0], np.cumsum(a["lig_count"])])
    ch = pd.read_csv(DATA / "chains.csv")
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv")
    s = pd.read_csv(ROOT / "results/sv/SV11/site_residues.csv")
    site = np.array(sorted(s[(s.site == "ligand pocket") & (s.in_site.astype(str) == "True")].uniprot_position
                           .astype(int)))
    # E2 reference centroid: wild-type EST chains with a complete ligand
    m = ch.merge(ct[["chain_key", "construct_mutations"]], on="chain_key")
    wt_est = m[(m.construct_mutations == "none") & (m.ligand_ids_used == "EST") & (m.ligand_incomplete == False)]  # noqa: E712
    idx = {k: i for i, k in enumerate(keys)}
    cents = np.array([a["lig_xyz"][ls[idx[k]]:ls[idx[k] + 1]].mean(0) for k in wt_est.chain_key])
    e2c = np.median(cents, axis=0)
    (DATA / "e2_centroid.json").write_text(json.dumps(
        {"definition": "per-coordinate median of EST heavy-atom centroids (core v3 frame), wild-type EST chains "
                       "with a complete ligand", "chains": list(wt_est.chain_key), "n_chains": int(len(wt_est)),
         "centroid_A": [float(v) for v in e2c],
         "max_distance_of_chain_centroid_to_median_A": float(np.linalg.norm(cents - e2c, axis=1).max())}, indent=1))
    only = sys.argv[1:]  # optional chain keys for a quick check
    jobs = []
    for i, k in enumerate(keys):
        if only and k not in only:
            continue
        jobs.append((k, a["rec_xyz"][rs[i]:rs[i + 1]].astype(float), a["rec_rad"][rs[i]:rs[i + 1]].astype(float),
                     a["rec_pos"][rs[i]:rs[i + 1]], a["lig_xyz"][ls[i]:ls[i + 1]].astype(float),
                     a["lig_rad"][ls[i]:ls[i + 1]].astype(float), site, e2c))
    t0 = time.time()
    if len(jobs) <= 4:
        res = [pocket_chain(j) for j in jobs]
    else:
        with Pool(10) as pool:
            res = list(pool.imap(pocket_chain, jobs, chunksize=4))
    df = pd.DataFrame(res)
    out = DATA / ("grid_per_chain.csv" if not only else "grid_check.csv")
    df.to_csv(out, index=False)
    print(f"{len(df)} chains in {time.time() - t0:.0f} s -> {out.relative_to(ROOT)}; E2 centroid {e2c.round(2)} "
          f"from {len(wt_est)} chains")


if __name__ == "__main__":
    main()
