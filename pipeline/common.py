"""Shared paths, loaders and rules for the pipeline (pipeline/01_*.py ... 11_*.py).

Inputs live in data/ (frozen tables and the deposited backbone coordinates); every stage writes to results/<stage>/.
Seed 42 for every random step. Distances in Angstrom.

Groups: cl7 and cl4 are the two main average-linkage groups of the 788 x 788 C-alpha RMSD matrix at 1.5 A (lid
closed over the pocket / lid in the coactivator groove); the five small groups (14 chains) are pooled as "minor".

Ordinary variation (the rule for calling an effect real): a difference counts only if it is larger than BOTH the
within-group standard deviation of the same quantity and the 95th percentile of the absolute difference between
copies of one crystal (chains of the same entry in the same group).
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("ESR1_ROOT", Path(__file__).resolve().parents[1]))
DATA = ROOT / "data"
RESULTS = ROOT / "results"
SEED = 42

INPUTS = {
    # structures
    "chain_table": "data/structures/chain_table.csv",
    "entries": "data/structures/entries.csv",
    "chain_status": "data/structures/chain_status.csv",
    "attrition": "data/structures/attrition_funnel.csv",
    "september_clusters": "data/structures/september_cluster_stats.csv",
    "october_clusters": "data/structures/october_cluster_stats.csv",
    "frozen_groups": "data/structures/frozen_groups.csv",
    "additives": "data/structures/additives_per_entry.csv",
    # coordinates
    "backbone": "data/coordinates/backbone.npz",
    "core": "data/coordinates/core_v3.json",
    "reference": "data/coordinates/reference.json",
    "helix_bounds": "data/coordinates/helix_boundaries.csv",
    # lid tables from the 28 Oct work
    "geometry": "data/lid/geometry_per_chain.csv",
    "mobility": "data/lid/mobility_per_chain.csv",
    "lattice": "data/lid/lattice_contacts_per_chain.csv",
    "cl4_subgroups": "data/lid/cl4_subgroups.csv",
    # outliers (five detectors, October)
    "outlier_flags": "data/outliers/outlier_flags.csv",
    "outlier_verdicts": "data/outliers/outlier_verdicts.csv",
    # maps made in October
    "october_map_quality": "data/maps/october_map_quality.csv",
    "october_pca_scores": "data/maps/october_coord_pca_scores.csv",
    # ligands and pockets
    "ligands_ranked": "data/ligands/ligands_ranked.csv",
    "ligand_descriptors": "data/ligands/ligand_descriptors.csv",
    "ligand_ids": "data/ligands/ligand_ids.csv",
    "chain_class": "data/ligands/chain_class.csv",
    "action_class": "data/ligands/action_class.csv",
    "p2rank": "data/ligands/p2rank_pockets.csv",
    "pocket_sites": "data/ligands/pocket_sites.csv",
    "pocket_filling": "data/ligands/pocket_filling_per_chain.csv",
    "pocket_filling_tests": "data/ligands/pocket_filling_tests.csv",
    "affinity": "data/ligands/affinity_per_ligand.csv",
    # ESR2
    "rmsd_esr2": "data/esr2/rmsd_esr2.npy",
    "rmsd_esr2_index": "data/esr2/rmsd_esr2_index.txt",
    "esr2_chains": "data/esr2/esr2_chain_table.csv",
    "rmsd_joint": "data/esr2/rmsd_joint.npy",
    "joint_index": "data/esr2/joint_index.csv",
    # produced by stage 02
    "aligned": "results/02_lid_groups/aligned_ca.npz",
    "rmsd": "results/02_lid_groups/rmsd.npy",
    "rmsd_without_lid": "results/02_lid_groups/rmsd_without_528_550.npy",
    "rmsd_index": "results/02_lid_groups/chains.csv",
    "groups": "results/02_lid_groups/groups.csv",
}


def p(key: str) -> Path:
    return ROOT / INPUTS[key]


def out_dir(stage: str) -> Path:
    d = RESULTS / stage
    d.mkdir(parents=True, exist_ok=True)
    return d


def read(key: str, **kw) -> pd.DataFrame:
    """CSV reader that keeps CCD code 'NA' (sodium) as text."""
    return pd.read_csv(p(key), keep_default_na=False, na_values=[""], **kw)


def group_of(cluster) -> str:
    return {7: "cl7", 4: "cl4"}.get(int(cluster), "minor")


def rmsd(key: str = "rmsd"):
    m = np.load(p(key))
    idx = read("rmsd_index")
    assert m.shape == (len(idx), len(idx))
    return m, idx


def additives_by_entry() -> pd.DataFrame:
    """Distinct additive CCD codes per entry (ions, solvents, buffers, cryoprotectants; ligands_ranked is_additive)."""
    lr = read("ligands_ranked")
    rows = [(e, r.ligand_id) for _, r in lr[lr.is_additive.astype(str) == "True"].iterrows()
            for e in str(r.entry_ids).split(";") if e and e != "nan"]
    a = pd.DataFrame(rows, columns=["entry_id", "ccd"])
    return a.groupby("entry_id").ccd.apply(lambda s: ";".join(sorted(set(s)))).rename("additive_codes").reset_index()


def chains() -> pd.DataFrame:
    """788 chains in matrix order: metadata, groups, lid geometry, action class, lid order / mobility, contacts."""
    _, idx = rmsd()
    t = idx[["chain_key", "entry_id"]].copy()
    ct = read("chain_table")
    t = t.merge(ct.drop(columns=["entry_id"]), on="chain_key", how="left", validate="one_to_one")
    gr = read("groups")[["chain_key", "cut_1.0", "cut_1.5"]]
    t = t.merge(gr, on="chain_key", how="left", validate="one_to_one")
    assert (t["cut_1.5"] == t.cluster_1p5).all(), "recomputed groups differ from the chain table"
    t["group"] = t["cut_1.5"].map(group_of)
    g = read("geometry")
    keep = ["chain_key", "dist_H12c_H11axis_A", "dist_H12c_H3axis_A", "angle_H12_H3_deg", "angle_H12_H11_deg",
            "hinge_H11C_H12_deg", "kink_H11_deg", "len_H12", "dist_H12c_H11axis_A_w544"]
    t = t.merge(g[keep], on="chain_key", how="left", validate="one_to_one")
    cc = read("chain_class")[["chain_key", "action_class"]]
    t = t.merge(cc, on="chain_key", how="left", validate="one_to_one")
    t["action_class"] = t.action_class.fillna("unlabelled")
    mob = read("mobility")[["chain_key", "h12_status", "mobility_h12_z", "mobility_z", "n_obs_538_547",
                            "tls_refined"]]
    t = t.merge(mob, on="chain_key", how="left", validate="one_to_one")
    lat = read("lattice")[["chain_key", "sym_atoms", "asu_atoms", "sym_receptor_atoms", "n_residues_sym_contact"]]
    t = t.merge(lat, on="chain_key", how="left", validate="one_to_one")
    of = read("outlier_flags")[["chain_key", "n_methods", "consensus_outlier"]]
    t = t.merge(of, on="chain_key", how="left", validate="one_to_one")
    ov = read("outlier_verdicts")[["chain_key", "verdict", "reasons"]]
    t = t.merge(ov.rename(columns={"verdict": "outlier_verdict", "reasons": "outlier_reasons"}), on="chain_key", how="left",
                validate="one_to_one")
    t = t.merge(additives_by_entry(), on="entry_id", how="left", validate="many_to_one")
    t["additive_codes"] = t.additive_codes.fillna("")
    lig = read("ligand_ids")[["ccd", "parent_pref_name"]]
    lr = read("ligands_ranked")[["ligand_id", "ligand_name"]]
    nm = lr.rename(columns={"ligand_id": "ccd", "ligand_name": "ccd_name"}).merge(lig.drop_duplicates("ccd"),
                                                                                   on="ccd", how="outer")
    nm["ligand_name"] = nm.parent_pref_name.fillna(nm.ccd_name)
    t = t.merge(nm[["ccd", "ligand_name"]].drop_duplicates("ccd").rename(columns={"ccd": "ligand_id_primary"}),
                on="ligand_id_primary", how="left", validate="many_to_one")
    assert len(t) == len(idx) and (t.chain_key.values == idx.chain_key.values).all()
    return t


def copy_diffs(df: pd.DataFrame, col: str) -> np.ndarray:
    """|difference| of `col` between every pair of chains of one entry that sit in the same group."""
    out = []
    for _, d in df.dropna(subset=[col]).groupby(["entry_id", "group"]):
        v = d[col].to_numpy(float)
        if len(v) > 1:
            i, j = np.triu_indices(len(v), 1)
            out.append(np.abs(v[i] - v[j]))
    return np.concatenate(out) if out else np.array([])


def ordinary_variation(df: pd.DataFrame, col: str, group: str | None = None) -> dict:
    d = df if group is None else df[df.group == group]
    cd = copy_diffs(d, col)
    sd = float(d[col].std()) if d[col].notna().sum() > 1 else float("nan")
    p95 = float(np.percentile(cd, 95)) if len(cd) else float("nan")
    return {"within_sd": sd, "copy_p95": p95, "n_copy_pairs": int(len(cd)),
            "threshold": float(np.nanmax([sd, p95])) if np.isfinite([sd, p95]).any() else float("nan")}


def bh(pvals) -> np.ndarray:
    pv = np.asarray(pvals, float)
    q = np.full_like(pv, np.nan)
    ok = np.isfinite(pv)
    x = pv[ok]
    if len(x):
        o = np.argsort(x)
        r = x[o] * len(x) / (np.arange(len(x)) + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        qq = np.empty_like(x)
        qq[o] = np.minimum(r, 1)
        q[ok] = qq
    return q


def masked_rmsd(x: np.ndarray, pos_mask=None) -> np.ndarray:
    """C-alpha RMSD over the positions observed in both chains of each pair (x: chains x positions x 3, NaN gaps)."""
    if pos_mask is not None:
        x = x[:, pos_mask]
    obs = np.isfinite(x[..., 0])
    n = len(x)
    s, c = np.zeros((n, n)), np.zeros((n, n))
    for k in range(x.shape[1]):
        o = obs[:, k]
        xp = np.where(o[:, None], x[:, k], 0.0)
        sq = (xp ** 2).sum(1)
        d2 = sq[:, None] + sq[None, :] - 2 * xp @ xp.T
        m = np.outer(o, o)
        s += np.where(m, np.maximum(d2, 0), 0)
        c += m
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.sqrt(s / c)
    np.fill_diagonal(r, 0)
    return r


def write_numbers(out: Path, values: dict, name: str = "numbers.json"):
    def conv(v):
        if isinstance(v, (np.integer,)):
            return int(v)
        if isinstance(v, (np.floating, float)):
            return None if not np.isfinite(v) else float(v)
        if isinstance(v, (np.bool_,)):
            return bool(v)
        if isinstance(v, dict):
            return {str(k): conv(x) for k, x in v.items()}
        if isinstance(v, (list, tuple, np.ndarray)):
            return [conv(x) for x in v]
        return v
    (out / name).write_text(json.dumps(conv(values), indent=1))
