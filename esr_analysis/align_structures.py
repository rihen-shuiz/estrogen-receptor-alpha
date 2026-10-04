"""Fit corresponding C-alpha atoms to an explicit reference chain.

Inputs are canonical-position coordinate mappings in angstroms. Kabsch fitting
uses shared core positions, then applies that same rigid transform to all atoms.
Outputs are aligned coordinates and fit details; no paths or reference policy
are inferred. Fitting describes geometry, not receptor identity or quality.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from pathlib import Path

import numpy as np


def load_coordinates(path: str | Path) -> dict[int, np.ndarray]:
    """Read a historical NPZ cache, rejecting invalid canonical positions/atoms."""
    with np.load(path, allow_pickle=False) as data:
        try:
            raw_positions = np.asarray(data["positions"])
            raw_coords = np.asarray(data["coords"])
        except KeyError as exc:
            raise ValueError(f"NPZ coordinate cache missing required array: {exc.args[0]}") from exc
    try:
        positions_as_float = raw_positions.astype(float)
        coords = raw_coords.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError("NPZ coordinate cache positions and coordinates must be numeric") from exc
    if positions_as_float.ndim != 1:
        raise ValueError("NPZ coordinate cache positions must be a 1D array")
    if not np.isfinite(positions_as_float).all():
        raise ValueError("NPZ coordinate cache positions must be finite")
    if not np.equal(positions_as_float, np.floor(positions_as_float)).all():
        raise ValueError("NPZ coordinate cache positions must be integer-like")
    positions = positions_as_float.astype(int)
    if np.any(positions < 1):
        raise ValueError("NPZ coordinate cache positions must be positive canonical indices")
    if len(np.unique(positions)) != len(positions):
        raise ValueError("NPZ coordinate cache positions must be unique")
    if coords.ndim != 2 or coords.shape[1:] != (3,):
        raise ValueError("NPZ coordinate cache coordinates must have shape N x 3")
    if len(positions) != len(coords):
        raise ValueError("NPZ coordinate cache positions and coordinate rows must have equal length")
    if not np.isfinite(coords).all():
        raise ValueError("NPZ coordinate cache coordinates must be finite")
    return {int(p): coords[i] for i, p in enumerate(positions)}


def kabsch_align(mobile: np.ndarray, reference: np.ndarray):
    """Return aligned atoms, proper rotation, mobile centroid, reference centroid.

    Atom rows must already correspond. Row-vector coordinates transform as
    ``(mobile - mobile_centroid) @ rotation + reference_centroid``.
    """
    if mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError("mobile/reference must both be N x 3 with matching shapes")
    cm = mobile.mean(axis=0)
    cr = reference.mean(axis=0)
    p = mobile - cm
    q = reference - cr
    h = p.T @ q
    u, s, vt = np.linalg.svd(h)
    r = u @ vt
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = u @ vt
    aligned = p @ r + cr
    return aligned, r, cm, cr


def apply_transform(coords: np.ndarray, r: np.ndarray, cm: np.ndarray, cr: np.ndarray):
    """Apply a fitted row-vector rigid transform to any atom selection."""
    return (coords - cm) @ r + cr


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Measure corresponding atoms in angstroms, without any additional fitting."""
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != 3:
        raise ValueError("RMSD requires corresponding atom arrays with equal N x 3 shapes")
    if len(a) == 0:
        return float("nan")
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def fit_to_reference(
    coords: Mapping[int, np.ndarray],
    reference: Mapping[int, np.ndarray],
    core_positions: Iterable[int],
    min_core: int = 60,
) -> tuple[dict[int, np.ndarray], dict]:
    """Fit the shared core; return every mobile atom and an auditable transform.

    Insufficient shared core returns an empty mapping and status
    ``insufficient_core`` with a reason. The caller chooses the reference and
    coverage policy. Historical alignment used a three-atom fitting guard after
    screening chains for at least 60 core atoms; pass ``min_core=3`` to replay it.
    """
    if min_core < 3:
        raise ValueError("min_core must be at least 3")
    common_core = sorted(set(core_positions) & set(reference) & set(coords))
    fit_info = {
        "status": "insufficient_core",
        "reason": None,
        "common_core_positions": common_core,
        "n_common_core_ca": len(common_core),
        "min_core": min_core,
        "core_alignment_rmsd_A": float("nan"),
        "rotation": None,
        "mobile_centroid": None,
        "reference_centroid": None,
    }
    if len(common_core) < min_core:
        fit_info["reason"] = f"{len(common_core)} shared core CA atoms; requires {min_core}"
        return {}, fit_info

    mobile = np.array([coords[p] for p in common_core])
    ref = np.array([reference[p] for p in common_core])
    aligned_core, rotation, mobile_centroid, reference_centroid = kabsch_align(mobile, ref)
    all_positions = sorted(coords)
    xyz = np.array([coords[p] for p in all_positions])
    aligned_xyz = apply_transform(xyz, rotation, mobile_centroid, reference_centroid)
    aligned_dict = {p: aligned_xyz[i] for i, p in enumerate(all_positions)}
    fit_info.update(
        status="ok",
        core_alignment_rmsd_A=rmsd(aligned_core, ref),
        rotation=rotation,
        mobile_centroid=mobile_centroid,
        reference_centroid=reference_centroid,
    )
    return aligned_dict, fit_info


def main():
    import json
    import pandas as pd
    from ._task_io import parser, prepare, coordinate_mapping, json_write
    from .run_io import read_table, positions
    p = parser('Block 09: conserved-core fit of explicit chains; saved transforms applied to every LBD CA')
    p.add_argument('--chain-keys', nargs='+', default=['3Q95|B|1', '7NFB|B|1'])
    args = p.parse_args()
    root, source, out, settings = prepare(args)
    mapping = read_table(source, 'residue_mapping.csv')
    reference = json.loads((source / 'reference.json').read_text())['chain_key']
    ref = coordinate_mapping(mapping, reference)
    canonical = np.arange(settings['esr1_lbd_start'], settings['esr1_lbd_end'] + 1)
    coords, transforms, summaries = [], {}, []
    for key in args.chain_keys:
        mobile = coordinate_mapping(mapping, key)
        aligned, transform = fit_to_reference(mobile, ref, positions(settings['esr1_core_ranges']), settings['min_core_ca_for_alignment'])
        if transform['status'] != 'ok':
            raise ValueError(f"Cannot fit {key}: {transform['reason']}")
        transforms[key] = transform
        coords.append([aligned.get(int(p), np.full(3, np.nan)) for p in canonical])
        summaries.append({'chain_key': key, 'n_common_core_ca': transform['n_common_core_ca'], 'core_alignment_rmsd_A': transform['core_alignment_rmsd_A']})
    np.savez_compressed(out / 'aligned_coordinates.npz', chain_keys=np.asarray(args.chain_keys), positions=canonical, coords=np.asarray(coords))
    json_write(out / 'alignment_transforms.json', transforms)
    pd.DataFrame(summaries).to_csv(out / 'global_alignment.csv', index=False)
    print(f'Block 09: {len(coords)} explicit chains; no pairwise refits')


if __name__ == "__main__":
    main()
