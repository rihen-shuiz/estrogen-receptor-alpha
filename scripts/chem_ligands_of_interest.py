"""Chemistry on ligands of interest only (tracker task C11).

Input: results/chain_ligand_of_interest.csv (selected chains, from task Q18) and
results/ligands.csv (SMILES). The chemistry set is the distinct ligand_id values with a
ligand in that file. For this set the script recomputes ECFP4 (Morgan radius 2, 2048 bits),
the Tanimoto matrix and K-means for K=3..7 with the same settings as C07 and C10
(n_init=20, random_state=42, silhouette on binary features, Euclidean distance).
Saved full-set matrices are loaded only to check the new values; they are not replaced.

Outputs (new files; earlier results are kept):
  results/ligands_of_interest.csv   one row per ligand of interest; matrix_row indexes the .npy files
  results/loi_fp_ecfp4.npy          ECFP4 bits for the clustered ligands
  results/loi_tanimoto_ecfp4.npy    Tanimoto matrix for the clustered ligands
  results/kmeans_scores_loi.csv     silhouette and inertia for K=3..7
  results/ligand_clusters_loi.csv   cluster labels for each K; chem_cluster = best-silhouette K

Run from anywhere: python scripts/chem_ligands_of_interest.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

ROOT = Path(__file__).resolve().parents[1]
K_RANGE = range(3, 8)
RADIUS, N_BITS = 2, 2048
SEED, N_INIT = 42, 20
GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS, fpSize=N_BITS)


def ecfp4(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return np.array(GENERATOR.GetFingerprintAsNumPy(mol), dtype=np.uint8)


def tanimoto(bits):
    bits = bits.astype(float)
    inter = bits @ bits.T
    size = bits.sum(axis=1)
    union = size[:, None] + size[None, :] - inter
    return np.divide(inter, union, out=np.zeros_like(inter), where=union > 0)


def main():
    RDLogger.DisableLog("rdApp.*")
    chains = pd.read_csv(ROOT / "results/chain_ligand_of_interest.csv", dtype=str, keep_default_na=False)
    selected = chains[chains.ligand_id != ""]
    counts = (selected.groupby("ligand_id")
              .agg(n_chains=("chain_key", "nunique"), n_entries=("entry_id", "nunique"))
              .reset_index())
    ligands = pd.read_csv(ROOT / "results/ligands.csv", dtype=str, keep_default_na=False)
    table = counts.merge(ligands[["ligand_id", "ligand_name", "smiles", "canonical_isomeric_smiles",
                                  "chemistry_status"]], on="ligand_id", how="left", validate="one_to_one")
    if table.smiles.isna().any():
        raise ValueError("ligand of interest missing from results/ligands.csv")
    table = table.sort_values("ligand_id", kind="stable").reset_index(drop=True)

    fingerprints = [ecfp4(smiles) for smiles in table.canonical_isomeric_smiles]
    table["in_clustering"] = [fp is not None for fp in fingerprints]
    matrix_row, n = [], 0
    for fp in fingerprints:
        matrix_row.append(n if fp is not None else "")
        n += fp is not None
    table["matrix_row"] = matrix_row
    clustered = table[table.in_clustering]
    bits = np.array([fp for fp in fingerprints if fp is not None])
    tanimoto_matrix = tanimoto(bits)

    saved_fp = np.load(ROOT / "results/fp_ecfp4.npy")
    saved_tanimoto = np.load(ROOT / "results/tanimoto_ecfp4.npy")
    index = pd.read_csv(ROOT / "results/fp_index.csv", dtype=str, keep_default_na=False)
    saved_row = dict(zip(index.ligand_id, index.row.astype(int)))
    rows_saved = [saved_row[ligand_id] for ligand_id in clustered.ligand_id]
    fp_matches_saved = np.array_equal(bits, saved_fp[rows_saved])
    tanimoto_max_diff = float(np.abs(tanimoto_matrix - saved_tanimoto[np.ix_(rows_saved, rows_saved)]).max())

    scores, labels = [], {}
    for k in K_RANGE:
        model = KMeans(n_clusters=k, n_init=N_INIT, random_state=SEED).fit(bits.astype(float))
        silhouette = silhouette_score(bits.astype(float), model.labels_, metric="euclidean")
        scores.append({"k": k, "silhouette": silhouette, "inertia": model.inertia_, "status": "tested", "reason": ""})
        labels[f"cluster_k{k}"] = model.labels_ + 1
    best_k = max(scores, key=lambda row: row["silhouette"])["k"]

    clusters = pd.DataFrame({"ligand_id": clustered.ligand_id.to_numpy(),
                             "chem_cluster": labels[f"cluster_k{best_k}"]})
    for k in K_RANGE:
        clusters[f"cluster_k{k}"] = labels[f"cluster_k{k}"]

    table.to_csv(ROOT / "results/ligands_of_interest.csv", index=False)
    np.save(ROOT / "results/loi_fp_ecfp4.npy", bits)
    np.save(ROOT / "results/loi_tanimoto_ecfp4.npy", tanimoto_matrix)
    pd.DataFrame(scores).to_csv(ROOT / "results/kmeans_scores_loi.csv", index=False)
    clusters.to_csv(ROOT / "results/ligand_clusters_loi.csv", index=False)

    print(f"ligands of interest: {len(table)}  clustered: {len(clustered)}  chain-ligand rows: {len(selected)}")
    print(f"ECFP4 bits equal saved rows: {fp_matches_saved}; Tanimoto max abs diff vs saved: {tanimoto_max_diff:.3g}")
    print(pd.DataFrame(scores).to_string(index=False))
    print(f"best silhouette K: {best_k}")


if __name__ == "__main__":
    main()
