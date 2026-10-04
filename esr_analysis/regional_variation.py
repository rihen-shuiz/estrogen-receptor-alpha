"""Block 13: measure structural ensemble variation after the saved conserved-core fit.

Coordinates are experimental ESR1 C-alpha observations, not an MD simulation.
The wide profiles use each position's available entry observations. MDAnalysis
uses only positions shared by every frame, in canonical order, without fitting.
B-factors remain raw square-angstrom values and do not establish dynamics.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .run_io import read_table, save_table, write_json


# MemoryReader stores float32 coordinates. Allow this quantisation when checking
# MDAnalysis against the original float64 NumPy profile (angstroms).
RMSF_ATOL_A = 1e-5
RMSF_RTOL = 1e-5


def ensemble_rmsf(coords):
    """Population RMSF around the mean; whole missing atoms contribute nothing."""
    observed = np.isfinite(coords).all(axis=2)
    counts = observed.sum(axis=0)
    mean = np.divide(
        np.where(observed[..., None], coords, 0).sum(axis=0), counts[:, None],
        out=np.full(coords.shape[1:], np.nan), where=counts[:, None] > 0,
    )
    squared = np.where(observed, np.sum((coords - mean) ** 2, axis=2), 0)
    rmsf = np.sqrt(np.divide(
        squared.sum(axis=0), counts, out=np.full(len(counts), np.nan), where=counts > 0,
    ))
    return counts, rmsf


def _entry(key):
    return key.split("|", 1)[0]


def _frames(selected, reference):
    """Place reference first, replacing another selection from its entry."""
    if not selected or len(set(selected)) != len(selected):
        raise ValueError("Selected chain index must be nonempty and unique")
    if len({_entry(key) for key in selected}) != len(selected):
        raise ValueError("Selected chain index must contain one chain per entry")
    replaced = [key for key in selected if _entry(key) == _entry(reference) and key != reference]
    frames = [reference] + [key for key in selected if _entry(key) != _entry(reference)]
    return frames, replaced


def _mdanalysis_rmsf(coords, positions, numpy_rmsf):
    """Only unavailable imports permit fallback; calculation errors propagate."""
    summary = {"rmsf_atol_A": RMSF_ATOL_A, "rmsf_rtol": RMSF_RTOL,
               "rmsf_agreement": None, "rmsf_max_abs_difference_A": None}
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis.rms import RMSF
    except ImportError as exc:
        return None, dict(summary, status="unavailable", reason=str(exc))
    if not len(positions):
        return None, dict(summary, status="no_common_atoms", version=mda.__version__)
    universe = mda.Universe.empty(
        len(positions), n_residues=len(positions), atom_resindex=np.arange(len(positions)),
    )
    universe.add_TopologyAttr("names", ["CA"] * len(positions))
    universe.add_TopologyAttr("resids", positions)
    universe.load_new(np.asarray(coords, dtype=np.float32), order="fac")
    values = RMSF(universe.atoms).run().results.rmsf
    summary.update(
        status="ok", version=mda.__version__,
        rmsf_agreement=bool(np.allclose(values, numpy_rmsf, atol=RMSF_ATOL_A, rtol=RMSF_RTOL)),
        rmsf_max_abs_difference_A=float(np.max(np.abs(values - numpy_rmsf))),
    )
    if not summary["rmsf_agreement"]:
        summary["status"] = "rmsf_disagreement"
    return values, summary


def run_local(root, out, settings):
    """Read aligned run artifacts; save all-chain displacements and entry profiles."""
    out = Path(out)
    with np.load(out / "aligned_coordinates.npz", allow_pickle=False) as data:
        keys = data["chain_keys"].astype(str).tolist()
        positions = np.asarray(data["positions"])
        coords = np.asarray(data["coords"], dtype=float)
    canonical = np.arange(settings["esr1_lbd_start"], settings["esr1_lbd_end"] + 1)
    if not np.array_equal(positions, canonical):
        raise ValueError("Aligned positions must contain the complete canonical ESR1 LBD in order")
    if coords.shape != (len(keys), len(positions), 3) or len(set(keys)) != len(keys):
        raise ValueError("Aligned coordinates require unique chains and shape (nchains, npositions, 3)")
    if np.isinf(coords).any() or np.any(np.isnan(coords).any(axis=2) != np.isnan(coords).all(axis=2)):
        raise ValueError("Missing aligned atoms must have three NaNs; observed atoms must be finite")
    reference = json.loads((out / "reference.json").read_text())["chain_key"]
    selected = (out / "rmsd_entry_index.txt").read_text().split()
    frames, replaced = _frames(selected, reference)
    lookup = {key: i for i, key in enumerate(keys)}
    if (set(frames) | set(selected)) - set(keys):
        raise ValueError("Reference/selected chain absent from aligned_coordinates.npz")
    # The saved transform table is a prerequisite, not a request to fit again.
    transforms = read_table(out, "global_alignment.csv")
    if transforms.chain_key.duplicated().any() or set(keys) - set(transforms.chain_key):
        raise ValueError("Every aligned chain needs one saved global alignment summary")

    mapping = read_table(out, "residue_mapping.csv")
    mapping = mapping.loc[mapping.chain_key.isin(keys) & mapping.esr1_position.isin(positions)]
    if mapping.duplicated(["chain_key", "esr1_position"]).any():
        raise ValueError("Residue mapping must be unique per chain and ESR1 position")
    observed = np.isfinite(coords).all(axis=2)
    ref = coords[lookup[reference]]
    displacement = np.linalg.norm(coords - ref, axis=2)
    chain_indices, position_indices = np.where(observed)
    local = pd.DataFrame({
        "chain_key": np.asarray(keys)[chain_indices], "esr1_position": positions[position_indices],
        "displacement_A": displacement[chain_indices, position_indices],
    }).merge(
        mapping[["chain_key", "esr1_position", "B_factor_A2", "has_ca"]],
        on=["chain_key", "esr1_position"], how="left", validate="many_to_one", indicator=True,
    )
    if (local._merge != "both").any() or not local.has_ca.astype(str).str.lower().eq("true").all():
        raise ValueError("Observed aligned atoms need observed residue mapping rows")
    save_table(out, "local_displacements.csv", local.drop(columns=["has_ca", "_merge"]))

    ensemble_indices = [lookup[key] for key in frames]
    ensemble = coords[ensemble_indices]
    counts, rmsf = ensemble_rmsf(ensemble)
    distances = displacement[ensemble_indices]
    distance_counts = np.isfinite(distances).sum(axis=0)
    mean_displacement = np.divide(
        np.nansum(distances, axis=0), distance_counts,
        out=np.full(len(positions), np.nan), where=distance_counts > 0,
    )
    median_b = local.loc[local.chain_key.isin(frames)].groupby("esr1_position").B_factor_A2.median()
    reference_b = mapping.loc[mapping.chain_key == reference].set_index("esr1_position").B_factor_A2
    profiles = pd.DataFrame({
        "esr1_position": positions, "n_observed": counts, "ensemble_RMSF_A": rmsf,
        "mean_displacement_A": mean_displacement, "median_B_A2": median_b.reindex(positions).to_numpy(),
        "reference_B_A2": reference_b.reindex(positions).to_numpy(),
    })
    save_table(out, "ensemble_profiles.csv", profiles)
    regions = {f"core_{start}_{end}": (positions >= start) & (positions <= end)
               for start, end in settings["esr1_core_ranges"]}
    start, end = settings["h12_range"]
    regions["H12"] = (positions >= start) & (positions <= end)
    regions["other_LBD"] = ~np.logical_or.reduce(list(regions.values()))
    regional = [{"region": name, "n_positions": int(mask.sum()),
                 "median_displacement_A": profiles.loc[mask, "mean_displacement_A"].median()}
                for name, mask in regions.items()]
    save_table(out, "regional_summary.csv", regional)

    common = np.isfinite(ensemble).all(axis=(0, 2))
    trajectory = ensemble[:, common].astype(np.float32)
    mda_values, mda_summary = _mdanalysis_rmsf(trajectory, positions[common], rmsf[common])
    np.savez_compressed(
        out / "common_atom_trajectory.npz", chain_keys=np.asarray(frames), positions=positions[common],
        coords=trajectory, numpy_RMSF_A=rmsf[common],
        MDAnalysis_RMSF_A=np.full(common.sum(), np.nan) if mda_values is None else mda_values,
    )
    save_table(out, "frame_map.csv", [{
        "frame": i, "chain_key": key, "pdb_id": _entry(key), "is_reference": key == reference,
        "selection": "reference" if key == reference else "selected_entry_chain",
    } for i, key in enumerate(frames)])
    summary = {
        "status": mda_summary["status"], "reference_chain_key": reference, "frames": len(frames),
        "selected_entries": len(selected), "replaced_selected_chain_keys": replaced,
        "frame_chain_keys": frames, "all_aligned_chains": len(keys),
        "n_positions": len(positions), "common_atoms_count": int(common.sum()),
        "removed_positions": positions[~common].tolist(), "mdanalysis": mda_summary,
        "rmsf_agreement": mda_summary["rmsf_agreement"],
        "cohort": "one experimental ESR1 chain per entry; reference first; replaces same-entry selection",
        "alignment": "saved conserved-core alignment; no additional fit",
        "trajectory": "common_atom_trajectory.npz; float32 frames are experimental structures, not MD",
        "profile_missingness": "per-position available observations; population RMSF; absent reference gives NaN displacement",
        "regional_calculation": "median of per-position ensemble mean reference displacements; n_positions is configured extent",
        "units": {"coordinates_displacement_RMSF": "angstrom", "B_factor": "angstrom squared"},
    }
    write_json(out / "local_summary.json", summary)
    if mda_summary["rmsf_agreement"] is False:
        raise RuntimeError("MDAnalysis RMSF disagrees with NumPy; inspect local_summary.json")
    print(f"Local variation: {len(frames)} entry frames, {int(common.sum())}/{len(positions)} common CA atoms; "
          f"MDAnalysis {summary['status']}. Profiles and displacements: {out}")
    return summary


def regional_bfactors(source_folder, output_folder):
    """Produce the broader block-13 normalized B/helix outputs using completed code.

    Inputs are the saved run's mapping/chains/entry/local-displacement tables and
    existing source-linked helix annotations. Population SD uses ddof=0. This is
    separate from ensemble coordinate variation; neither is interpreted as MD.
    """
    from ._regional_bfactors import run
    run(source_folder=Path(source_folder), output_folder=Path(output_folder), annotation_file=Path(source_folder) / "helix_annotations.csv")
    return read_table(Path(output_folder), "helix_chain_summaries.csv")


def main():
    from ._task_io import parser, prepare
    from ._regional_bfactors import run
    args = parser('Block 13: existing population RMSF and coverage-aware within-chain B z; experimental ensemble, not MD').parse_args()
    root, source, out, settings = prepare(args)
    for name in ['aligned_coordinates.npz', 'reference.json', 'rmsd_entry_index.txt', 'global_alignment.csv', 'residue_mapping.csv']:
        path = source / name
        if not path.is_file():
            raise FileNotFoundError(f'Block 13 prerequisite missing: {path}')
    import shutil
    for name in ['aligned_coordinates.npz', 'reference.json', 'rmsd_entry_index.txt', 'global_alignment.csv', 'residue_mapping.csv']:
        shutil.copy2(source / name, out / name)
    run_local(root, out, settings)
    run(source_folder=source, output_folder=out, annotation_file=source / 'helix_annotations.csv')


if __name__ == "__main__":
    main()
