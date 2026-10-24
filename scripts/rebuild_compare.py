"""Compare every rebuilt upstream file (results/rebuild/) with its frozen counterpart.

Targets: Q04 (lbd_chains + lbd_extract/*), Q06 entries, Q08 ligands, C01 parse report, C02-C04 arrays and
fp_index, C05 UMAP, C06 t-SNE, C07/C10 K-means, E03 numbering map, E04 matched pairs, X01 ESR2 entries,
X13 ChEMBL download, X14 ChEMBL-PDB matches.
Tables: row counts, key sets, per-column exact or numeric (|d| <= 1e-6) equality. Arrays: shape, max |d|.
Clustering: ARI per K plus the score table. Embeddings (stochastic): finite-row mask, Procrustes disparity
and mean 10-nearest-neighbour overlap between frozen and rebuilt coordinates.
Output: results/rebuild/rebuild_report.csv (target, frozen_path, rebuilt_path, script, status, detail)
Run from the repo root: PYTHONPATH=. python scripts/rebuild_compare.py
"""
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy.spatial import procrustes
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors

from rebuild_common import CCD_LOCAL, OUT, REPORT, ROOT, compare_arrays, compare_frames, read_csv, report


def table(target, frozen, key, script, rebuilt=None, ignore=()):
    rebuilt = rebuilt or "results/rebuild/" + frozen.split("/", 1)[1]
    if not (ROOT / rebuilt).exists():
        report(target, frozen, rebuilt, script, "not_rebuilt", "rebuilt file missing")
        return
    status, detail = compare_frames(read_csv(ROOT / frozen), read_csv(ROOT / rebuilt), key, ignore)
    report(target, frozen, rebuilt, script, status, detail)


def array(target, name, script):
    frozen, rebuilt = f"results/{name}", f"results/rebuild/{name}"
    status, detail = compare_arrays(np.load(ROOT / frozen), np.load(ROOT / rebuilt))
    report(target, frozen, rebuilt, script, status, detail)


def ligands():
    frozen, rebuilt, script = "results/ligands.csv", "results/rebuild/ligands.csv", "scripts/rebuild_ligands.py"
    fz, rb = read_csv(ROOT / frozen), read_csv(ROOT / rebuilt)
    status, detail = compare_frames(fz, rb, "ligand_id", ignore=["ccd_sha256"])
    # Local CCD files are checked out by git with core.autocrlf=true (CRLF); the frozen hashes are of LF bytes.
    merged = fz.merge(rb, on="ligand_id", suffixes=("_f", "_r"))
    differ = merged[merged.ccd_sha256_f != merged.ccd_sha256_r]
    lf_ok = 0
    for row in differ.itertuples():
        path = CCD_LOCAL / f"{row.ligand_id}.cif"
        if path.exists() and hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest() == row.ccd_sha256_f:
            lf_ok += 1
    detail += (f"; ccd_sha256: {len(differ)} differ, all local CCD files checked out with CRLF (git core.autocrlf);"
               f" LF-normalised hash equals frozen for {lf_ok}/{len(differ)}")
    if status in ("identical", "equivalent") and lf_ok != len(differ):
        status = "differs"
    report("Q08 ligands", frozen, rebuilt, script, status, detail)


def kmeans(target, suffix, script, frozen_suffix=None):
    frozen_suffix = suffix if frozen_suffix is None else frozen_suffix
    fs, rs = f"results/kmeans_scores{frozen_suffix}.csv", f"results/rebuild/kmeans_scores{suffix}.csv"
    fl, rl = f"results/ligand_clusters{frozen_suffix}.csv", f"results/rebuild/ligand_clusters{suffix}.csv"
    s_status, s_detail = compare_frames(read_csv(ROOT / fs), read_csv(ROOT / rs), "k")
    a, b = read_csv(ROOT / fl), read_csv(ROOT / rl)
    a["ligand_id"], b["ligand_id"] = a.ligand_id.astype(str), b.ligand_id.astype(str)
    m = a.merge(b, on="ligand_id", suffixes=("_f", "_r"))
    aris = {c: adjusted_rand_score(m[c + "_f"], m[c + "_r"]) for c in a.columns if c != "ligand_id"}
    l_status, l_detail = compare_frames(a, b, "ligand_id")
    detail = (f"scores: {s_status} ({s_detail}); labels: {l_status} ({l_detail}); common ligands={len(m)}; ARI "
              + ", ".join(f"{c}={v:.3f}" for c, v in aris.items()))
    statuses = {s_status, l_status}
    status = "identical" if statuses == {"identical"} else "equivalent" if statuses <= {"identical", "equivalent"} else "differs"
    report(target, f"{fs} + {fl}", f"{rs} + {rl}", script, status, detail)


def embedding(target, name, cols, script):
    frozen, rebuilt = f"results/{name}", f"results/rebuild/{name}"
    a, b = read_csv(ROOT / frozen), read_csv(ROOT / rebuilt)
    a["ligand_id"], b["ligand_id"] = a.ligand_id.astype(str), b.ligand_id.astype(str)
    m = a.merge(b, on="ligand_id", suffixes=("_f", "_r"))
    fa = m[[c + "_f" for c in cols]].to_numpy(float)
    fb = m[[c + "_r" for c in cols]].to_numpy(float)
    ok_a, ok_b = np.isfinite(fa).all(1), np.isfinite(fb).all(1)
    both = ok_a & ok_b
    _, _, disparity = procrustes(fa[both], fb[both])
    k = 10
    na = NearestNeighbors(n_neighbors=k + 1).fit(fa[both]).kneighbors(fa[both], return_distance=False)[:, 1:]
    nb = NearestNeighbors(n_neighbors=k + 1).fit(fb[both]).kneighbors(fb[both], return_distance=False)[:, 1:]
    overlap = np.mean([len(set(x) & set(y)) / k for x, y in zip(na, nb)])
    nonfinite_f = sorted(m.ligand_id[~ok_a])
    nonfinite_r = sorted(m.ligand_id[~ok_b])
    detail = (f"stochastic embedding, coordinates not expected to match; rows frozen={len(a)} rebuilt={len(b)}; "
              f"finite frozen={int(ok_a.sum())} rebuilt={int(ok_b.sum())}; non-finite rows same set={nonfinite_f == nonfinite_r} "
              f"({','.join(nonfinite_f)}); max|d| coords={np.nanmax(np.abs(fa - fb)):.3g}; Procrustes disparity={disparity:.3f}; "
              f"mean {k}-NN overlap={overlap:.3f}")
    report(target, frozen, rebuilt, script, "differs", detail)


def main():
    if REPORT.exists():
        REPORT.unlink()
    s = "scripts/rebuild_lbd_extract.py"
    table("Q04 lbd_chains", "results/lbd_chains.csv", "chain_key", s)
    table("Q04 lbd_extract/entries", "results/lbd_extract/entries.csv", "pdb_id", s, ignore=["source_path"])
    table("Q04 lbd_extract/receptor_chains", "results/lbd_extract/receptor_chains.csv", "chain_key", s)
    table("Q04 lbd_extract/residue_mapping", "results/lbd_extract/residue_mapping.csv", ["chain_key", "label_seq_id"], s)
    table("Q06 entries", "results/entries.csv", "entry_id", "scripts/rebuild_entries_qc.py")
    ligands()
    table("C01 smiles_parse_report", "results/smiles_parse_report.csv", "ligand_id", "scripts/rebuild_ligands.py")
    s = "scripts/rebuild_fingerprints.py"
    table("C02/C03 fp_index", "results/fp_index.csv", "row", s)
    array("C02 fp_ecfp4", "fp_ecfp4.npy", s)
    array("C03 fp_maccs", "fp_maccs.npy", s)
    array("C04 tanimoto_ecfp4", "tanimoto_ecfp4.npy", s)
    s = "scripts/rebuild_embeddings.py"
    embedding("C05 umap_coords", "umap_coords.csv", ["umap1", "umap2"], s)
    embedding("C06 tsne_coords", "tsne_coords.csv", ["tsne1", "tsne2"], s)
    s = "scripts/rebuild_kmeans.py"
    kmeans("C07 kmeans all ligands", "", s)
    kmeans("C10 kmeans noadd (Q16 is_additive)", "_noadd", s)
    kmeans("C10 variant: settings solvent_additive_codes only", "_noadd_settings", s, frozen_suffix="_noadd")
    s = "scripts/rebuild_esr1_esr2_maps.py"
    table("E03 numbering_map_esr1_esr2", "results/numbering_map_esr1_esr2.csv", ["esr1_pos", "esr2_pos"], s)
    table("E04 matched_pairs", "results/matched_pairs.csv", ["esr1_entry", "esr2_entry"], s)
    table("X01 esr2_entries", "results/esr2_entries.csv", "entry_id", "scripts/rebuild_esr2_entries.py")
    s = "scripts/rebuild_chembl_activities.py"
    if (OUT / "chembl_esr1_activities.csv").exists():
        table("X13 chembl_esr1_activities", "data/chembl_esr1_activities.csv", "activity_id", s,
              rebuilt="results/rebuild/chembl_esr1_activities.csv")
    else:
        note = (OUT / "chembl_download_incomplete.txt")
        report("X13 chembl_esr1_activities", "data/chembl_esr1_activities.csv", "results/rebuild/chembl_esr1_activities.csv",
               s, "not_rebuilt", note.read_text().strip() if note.exists() else "download not run")
    table("X14 chembl_pdb_matches", "results/chembl_pdb_matches.csv", ["ligand_id", "molecule_chembl_id"],
          "scripts/rebuild_chembl_pdb_matches.py")
    print(pd.read_csv(REPORT).status.value_counts().to_dict())


if __name__ == "__main__":
    main()
