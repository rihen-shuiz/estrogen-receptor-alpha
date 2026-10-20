"""Tracker S07: all-heavy-atom RMSD on a 50-entry subset, compared with the C-alpha matrix (validation).

Subset: one chain per entry, the historical author chain (esr_analysis/resources/frozen/rmsd_entry_index.txt).
50 entries stratified by the S08 cut_1.5 cluster (results/entry_clusters.csv): every entry from a cluster with
fewer than 10 entries, the rest drawn in proportion from the large clusters with a fixed seed.
Coordinates: deposited mmCIFs in data/raw (downloaded from RCSB when missing, compared with data/raw/sha256sums.txt).
Atoms are mapped to ESR1 positions through results/lbd_extract/residue_mapping.csv (label_asym_id, label_seq_id)
and moved with the saved S03 core-fit transform (results/alignment/alignment_transforms.json); no refit.
Check: the transformed C-alpha atoms must equal results/alignment/aligned_coordinates.npz.
Pair distance: polymer heavy atoms (no H/D; modified residues kept) in the LBD range that both chains share by (ESR1 position, residue name,
atom name). Residues with different residue names in the two chains are skipped, so the comparison is atom for atom.
ca_rmsd is read from results/rmsd_entry.npy (S04, not recomputed); ca_rmsd_same_residues uses only the residues
that entered the all-atom value.
Outputs: results/rmsd_allatom_subset.csv (one row per pair) and results/rmsd_allatom_subset_summary.json
(Spearman between C-alpha and all-atom over pairs, per-chain neighbour ranking, cluster agreement).
Run from the repo root: PYTHONPATH=. .venv/Scripts/python scripts/allatom_rmsd_subset.py
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet, fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score

from esr_analysis.align_structures import apply_transform
from esr_analysis.download_structures import download_one
from esr_analysis.read_structure import read_cif

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
RAW = ROOT / "data/raw"
N_ENTRIES = 50
SMALL_CLUSTER = 10  # clusters with fewer entries than this are taken whole
SEED = 20261020
CUT = "cut_1.5"
MATCH_TOL_A = 1e-3  # transformed CA vs saved aligned CA


def choose_subset():
    historical = (ROOT / "esr_analysis/resources/frozen/rmsd_entry_index.txt").read_text().split()
    clusters = pd.read_csv(ROOT / "results/entry_clusters.csv")
    pool = clusters[clusters.chain_key.isin(historical)][["entry_id", "chain_key", CUT]]
    sizes = pool[CUT].value_counts()
    small = pool[pool[CUT].isin(sizes[sizes < SMALL_CLUSTER].index)]
    large = pool[~pool.index.isin(small.index)]
    remaining = N_ENTRIES - len(small)
    shares = large[CUT].value_counts()
    quota = (shares / shares.sum() * remaining).round().astype(int)
    quota.iloc[0] += remaining - quota.sum()
    rng = np.random.default_rng(SEED)
    picked = [small]
    for cluster, n in quota.items():
        members = large[large[CUT] == cluster]
        picked.append(members.iloc[np.sort(rng.choice(len(members), size=n, replace=False))])
    subset = pd.concat(picked).sort_values("chain_key").reset_index(drop=True)
    assert len(subset) == N_ENTRIES and subset.entry_id.is_unique
    return subset


def fetch(entry_ids):
    sums = dict(line.split()[::-1] for line in (RAW / "sha256sums.txt").read_text().splitlines() if line.strip())
    status = {}
    for pdb_id in entry_ids:
        path = RAW / f"{pdb_id}.cif"
        ok, how = download_one(pdb_id, path)
        if not ok:
            raise RuntimeError(f"{pdb_id}: download failed ({how})")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        status[pdb_id] = {"download": how, "sha256_matches_manifest": sums.get(path.name) == digest}
    return status


def chain_atoms(chain_key, mapping, transform, aligned_ca):
    pdb_id, asym, model = chain_key.split("|")
    _, atoms = read_cif(RAW / f"{pdb_id}.cif")
    atoms = atoms[(atoms.pdbx_PDB_model_num == model) & (atoms.label_asym_id == asym)
                  & ~atoms.type_symbol.isin(["H", "D"])]
    rows = mapping[mapping.chain_key == chain_key].dropna(subset=["esr1_position"])
    position = dict(zip(rows.label_seq_id.astype(str), rows.esr1_position.astype(int)))
    atoms = atoms.assign(position=atoms.label_seq_id.map(position)).dropna(subset=["position"])
    atoms = atoms[atoms.position.between(SETTINGS["esr1_lbd_start"], SETTINGS["esr1_lbd_end"])]
    xyz = apply_transform(atoms[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float), np.array(transform["rotation"]),
                          np.array(transform["mobile_centroid"]), np.array(transform["reference_centroid"]))
    keys = list(zip(atoms.position.astype(int), atoms.label_comp_id, atoms.label_atom_id))
    table = dict(zip(keys, xyz))
    ca = {k[0]: v for k, v in table.items() if k[2] == "CA"}
    worst = max(np.abs(ca[p] - aligned_ca[p]).max() for p in ca if p in aligned_ca)
    missing = set(aligned_ca) - set(ca)
    if worst > MATCH_TOL_A or missing:
        raise ValueError(f"{chain_key}: transformed CA disagree with aligned_coordinates.npz "
                         f"(max {worst:.4f} A, {len(missing)} saved CA missing)")
    return table


def pair_rmsd(a, b):
    residue_a = {(p, r) for p, r, _ in a}
    residue_b = {(p, r) for p, r, _ in b}
    same = {p for p, r in residue_a & residue_b}
    shared = [k for k in a if k in b and k[0] in same]
    ca = [k for k in shared if k[2] == "CA"]
    d = np.array([a[k] - b[k] for k in shared])
    dca = np.array([a[k] - b[k] for k in ca])
    return (float(np.sqrt((d ** 2).sum(1).mean())), len(shared),
            float(np.sqrt((dca ** 2).sum(1).mean())), len(ca), len(residue_a) - len(residue_a & residue_b))


def main():
    subset = choose_subset()
    status = fetch(subset.entry_id)
    keys = subset.chain_key.tolist()

    matrix_keys = (ROOT / "results/alignment/rmsd_entry_index.txt").read_text().split()
    matrix = np.load(ROOT / "results/rmsd_entry.npy")
    rows = [matrix_keys.index(k) for k in keys]
    ca_matrix = matrix[np.ix_(rows, rows)]

    transforms = json.loads((ROOT / "results/alignment/alignment_transforms.json").read_text())
    npz = np.load(ROOT / "results/alignment/aligned_coordinates.npz", allow_pickle=False)
    npz_keys = npz["chain_keys"].astype(str).tolist()
    mapping = pd.read_csv(ROOT / "results/lbd_extract/residue_mapping.csv", low_memory=False,
                          usecols=["chain_key", "label_seq_id", "esr1_position"], dtype={"label_seq_id": str})
    mapping = mapping[mapping.chain_key.isin(keys)]
    atoms = {}
    for key in keys:
        coords = npz["coords"][npz_keys.index(key)]
        aligned_ca = {int(p): c for p, c in zip(npz["positions"], coords) if np.isfinite(c).all()}
        atoms[key] = chain_atoms(key, mapping, transforms[key], aligned_ca)

    records = []
    allatom = np.zeros_like(ca_matrix)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            value, n_atoms, ca_same, n_ca, n_skipped = pair_rmsd(atoms[keys[i]], atoms[keys[j]])
            allatom[i, j] = allatom[j, i] = value
            records.append({"chain_a": keys[i], "chain_b": keys[j], "entry_a": keys[i].split("|")[0],
                            "entry_b": keys[j].split("|")[0], "ca_rmsd": ca_matrix[i, j], "allatom_rmsd": value,
                            "ca_rmsd_same_residues": ca_same, "n_heavy_atoms": n_atoms, "n_residues": n_ca,
                            "n_residues_skipped_name_differs": n_skipped,
                            f"{CUT}_a": subset[CUT][i], f"{CUT}_b": subset[CUT][j]})
    pairs = pd.DataFrame(records)
    pairs.to_csv(ROOT / "results/rmsd_allatom_subset.csv", index=False)

    rho, p = spearmanr(pairs.ca_rmsd, pairs.allatom_rmsd)
    rho_same, _ = spearmanr(pairs.ca_rmsd_same_residues, pairs.allatom_rmsd)
    per_chain = [spearmanr(np.delete(ca_matrix[i], i), np.delete(allatom[i], i))[0] for i in range(len(keys))]
    z_ca = linkage(squareform(ca_matrix, checks=False), method="average")
    z_aa = linkage(squareform(allatom, checks=False), method="average")
    agreement = {}
    for k in [2, 3, 4, 5, 7]:
        agreement[str(k)] = float(adjusted_rand_score(fcluster(z_ca, k, "maxclust"), fcluster(z_aa, k, "maxclust")))
    cophenetic_rho = spearmanr(cophenet(z_ca), cophenet(z_aa))[0]
    summary = {
        "task": "S07", "n_entries": len(keys), "n_pairs": len(pairs), "seed": SEED, "cluster_column": CUT,
        "subset_by_cluster": {str(k): int(v) for k, v in subset[CUT].value_counts().sort_index().items()},
        "chains": keys,
        "spearman_ca_vs_allatom_pairs": {"rho": float(rho), "p": float(p)},
        "spearman_ca_same_residues_vs_allatom_pairs": float(rho_same),
        "per_chain_neighbour_spearman": {"median": float(np.median(per_chain)), "min": float(np.min(per_chain)),
                                         "worst_chain": keys[int(np.argmin(per_chain))]},
        "cophenetic_spearman_ca_vs_allatom": float(cophenetic_rho),
        "adjusted_rand_same_n_clusters": agreement,
        "allatom_rmsd_A": {"median": float(pairs.allatom_rmsd.median()), "max": float(pairs.allatom_rmsd.max())},
        "ca_rmsd_A": {"median": float(pairs.ca_rmsd.median()), "max": float(pairs.ca_rmsd.max())},
        "n_heavy_atoms_min": int(pairs.n_heavy_atoms.min()),
        "downloads": status,
        "note": "same S03 core-fit superposition for both measures; all-atom adds side chains, no refit",
    }
    (ROOT / "results/rmsd_allatom_subset_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"S07: {len(keys)} chains, {len(pairs)} pairs; Spearman CA vs all-atom rho {rho:.3f} "
          f"(same residues {rho_same:.3f}); per-chain median {np.median(per_chain):.3f}; "
          f"ARI {agreement}; sha256 mismatches {sum(not s['sha256_matches_manifest'] for s in status.values())}")


if __name__ == "__main__":
    main()
