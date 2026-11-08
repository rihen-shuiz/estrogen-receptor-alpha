"""SV00 step 3 + 5: counts against the original run and the reproduction gate.

Counts (step 3): from the regenerated Q04 extraction (data/sv/lbd_chains.csv, SV00_extract.py) versus the frozen
results/lbd_chains.csv and the numbers in the run file (ESR1 487 entries; 901 chains in 444 entries; 788 kept chains
in 351 entries; 113 excluded with reasons; ESR2 39 entries, 86 chains).
Gate (step 5), the existing code path, imported unchanged:
  S03 scripts/align_selected_chains.fit_all (esr_analysis.align_structures.fit_to_reference, reference
      esr_analysis/resources/frozen/reference.json, core v3 from alignment_core_v3.json, min core from settings),
  S04 esr_analysis.calculate_rmsd.compute_distances on esr_analysis._task_io.aligned_caches (min_pairwise_common_ca),
  S08 esr_analysis.cluster_structures.structural_groups (average linkage, cuts 1.0/1.5/2.0/2.5, primary 1.5 A),
fed with the regenerated residue mapping instead of results/lbd_extract/residue_mapping.csv. Chains = the S03 rule
(identity esr1 and eligible_historical) applied to the regenerated chain table.
Compared with results/rmsd_entry.npy (+ results/rmsd_entry_index.csv) and results/entry_clusters.csv.
Differing chain (fixed before looking at results): a kept chain whose regenerated LBD C-alpha set or coordinates differ
from the frozen residue mapping (any position with a different has_ca, or |dxyz| > 0.0005 A, half the mmCIF rounding
step), or whose 1.5 A cluster differs from the frozen one after matching each recomputed cluster to the frozen cluster
it shares most chains with, or which is missing on one side.
Verdict: pass = ARI(1.5 A) == 1.0; [DEFAULT] tolerated = <= 8 differing chains, every one in an entry whose latest
revision (data/sv/manifest.csv) is after the original download date; otherwise BLOCKED for coordinates.
Outputs: results/sv/SV00/{counts.csv, excluded_chains.csv, coordinate_diff_by_chain.csv, rmsd_recomputed.npy,
rmsd_recomputed_index.csv, clusters_recomputed.csv, gate.json, gate_chains.csv},
data/sv/gate/{aligned_coordinates.npz, alignment_transforms.json, global_alignment.csv, overlap.npy}.
Run: .venv/bin/python -I scripts/sv/SV00_gate.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import adjusted_rand_score  # noqa: E402

import align_selected_chains as S03  # noqa: E402
from esr_analysis._task_io import aligned_caches, coordinate_mapping, json_write  # noqa: E402
from esr_analysis.calculate_rmsd import compute_distances  # noqa: E402
from esr_analysis.cluster_structures import structural_groups  # noqa: E402

OUT = ROOT / "results/sv/SV00"
GATE = ROOT / "data/sv/gate"
SETTINGS = S03.SETTINGS
CUTS = [1.0, 1.5, 2.0, 2.5]
COORD_TOL_A = 0.0005
MAX_DIFFER = 8
EXPECTED = {"esr1_entries": 487, "esr1_chains": 901, "esr1_chain_entries": 444, "kept_chains": 788,
            "kept_entries": 351, "excluded_chains": 113, "esr2_entries": 39, "esr2_chains": 86}


def read(path, **kw):
    return pd.read_csv(path, keep_default_na=False, na_values=[""], low_memory=False, **kw)


def chain_counts(lbd, esr1_ids, esr2_ids):
    e1 = lbd[lbd.identity == "esr1"]
    kept = e1[e1.eligible_historical.astype(str).str.lower() == "true"]
    e2 = lbd[lbd.identity == "esr2"]
    return {"esr1_entries": len(esr1_ids), "esr1_chains": len(e1), "esr1_chain_entries": e1.entry_id.nunique(),
            "kept_chains": len(kept), "kept_entries": kept.entry_id.nunique(), "excluded_chains": len(e1) - len(kept),
            "esr2_entries": len(esr2_ids), "esr2_chains": len(e2[e2.entry_id.isin(esr2_ids)])}


def cluster_mismatch(new, old):
    """Chains whose recomputed cluster is not the frozen cluster it overlaps most."""
    tab = pd.crosstab(new, old)
    best = tab.idxmax(axis=1)
    return new.map(best) != old


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    GATE.mkdir(parents=True, exist_ok=True)
    esr1_ids = [x.strip().upper() for x in (ROOT / "data/pdb_ids_esr1.txt").read_text().split()]
    esr2_ids = [x.strip().upper() for x in (ROOT / "data/pdb_ids_esr2.txt").read_text().split()]
    lbd_new = read(ROOT / "data/sv/lbd_chains.csv")
    lbd_old = read(ROOT / "results/lbd_chains.csv")
    c_new, c_old = chain_counts(lbd_new, esr1_ids, esr2_ids), chain_counts(lbd_old, esr1_ids, esr2_ids)
    counts = pd.DataFrame([{"quantity": k, "run_file": EXPECTED[k], "frozen_tables": c_old[k], "regenerated": c_new[k],
                            "match": EXPECTED[k] == c_old[k] == c_new[k]} for k in EXPECTED])
    counts.to_csv(OUT / "counts.csv", index=False)

    # excluded ESR1 chains, each with its reason (regenerated and frozen)
    reasons_old = read(ROOT / "results/chain_status.csv")[["chain_key", "completeness_status", "completeness_reason"]]
    e1n = lbd_new[lbd_new.identity == "esr1"]
    exc = e1n[e1n.eligible_historical.astype(str).str.lower() != "true"][
        ["chain_key", "entry_id", "auth_asym_id", "mapped_lbd_ca", "mapped_core_ca", "eligibility_reason"]]
    exc = exc.merge(reasons_old, on="chain_key", how="outer", indicator=True)
    exc["_merge"] = exc["_merge"].map({"both": "both", "left_only": "regenerated_only", "right_only": "frozen_only"})
    exc = exc[(exc._merge != "frozen_only") | (exc.completeness_status != "complete")]
    exc.rename(columns={"_merge": "present_in"}).to_csv(OUT / "excluded_chains.csv", index=False)

    # coordinates: regenerated vs frozen residue mapping, kept chains
    frozen_keys = read(ROOT / "results/rmsd_entry_index.csv").chain_key.tolist()
    kept_new = sorted(e1n[e1n.eligible_historical.astype(str).str.lower() == "true"].chain_key)
    cols = ["chain_key", "esr1_position", "has_ca", "x_A", "y_A", "z_A"]
    m_new = read(ROOT / "data/sv/lbd_extract/residue_mapping.csv", usecols=cols)
    m_old = read(ROOT / "results/lbd_extract/residue_mapping.csv", usecols=cols)
    allkeys = sorted(set(frozen_keys) | set(kept_new))
    m_new, m_old = m_new[m_new.chain_key.isin(allkeys)], m_old[m_old.chain_key.isin(allkeys)]
    lbd_rng = (SETTINGS["esr1_lbd_start"], SETTINGS["esr1_lbd_end"])
    rows = []
    gn, go = dict(tuple(m_new.groupby("chain_key"))), dict(tuple(m_old.groupby("chain_key")))
    for key in allkeys:
        a, b = gn.get(key), go.get(key)
        if a is None or b is None:
            rows.append({"chain_key": key, "status": "missing_regenerated" if a is None else "missing_frozen"})
            continue
        a = a[a.has_ca.astype(bool) & a.esr1_position.between(*lbd_rng)].set_index("esr1_position")
        b = b[b.has_ca.astype(bool) & b.esr1_position.between(*lbd_rng)].set_index("esr1_position")
        common = a.index.intersection(b.index)
        d = np.abs(a.loc[common, ["x_A", "y_A", "z_A"]].to_numpy(float) - b.loc[common, ["x_A", "y_A", "z_A"]].to_numpy(float))
        dmax = float(d.max()) if d.size else np.nan
        only_new, only_old = len(a.index.difference(b.index)), len(b.index.difference(a.index))
        same = only_new == 0 and only_old == 0 and (np.isnan(dmax) or dmax <= COORD_TOL_A)
        rows.append({"chain_key": key, "status": "identical" if same else "different", "n_ca_regenerated": len(a),
                     "n_ca_frozen": len(b), "ca_only_regenerated": only_new, "ca_only_frozen": only_old,
                     "max_abs_dxyz_A": dmax})
    cdiff = pd.DataFrame(rows)
    cdiff.to_csv(OUT / "coordinate_diff_by_chain.csv", index=False)

    # S03 fit with the regenerated mapping (existing functions)
    mapping = read(ROOT / "data/sv/lbd_extract/residue_mapping.csv", dtype={"chain_key": str})
    reference_key = json.loads(S03.REFERENCE.read_text())["chain_key"]
    reference = coordinate_mapping(mapping, reference_key)
    mobiles = {key: coordinate_mapping(mapping, key) for key in kept_new}
    coords, transforms, summary = S03.fit_all(mobiles, reference, reference_key, kept_new, S03.CORE_CONFIG["core_ranges"])
    canonical = np.arange(S03.LBD[0], S03.LBD[1] + 1)
    fits = GATE / "aligned_coordinates.npz"
    np.savez_compressed(fits, chain_keys=np.asarray(kept_new), positions=canonical, coords=coords)
    json_write(GATE / "alignment_transforms.json", transforms)
    summary.to_csv(GATE / "global_alignment.csv", index=False)
    old_fit = np.load(ROOT / "results/alignment/aligned_coordinates.npz", allow_pickle=False)
    old_keys = old_fit["chain_keys"].astype(str).tolist()
    shared = [k for k in kept_new if k in set(old_keys)]
    ia, ib = [kept_new.index(k) for k in shared], [old_keys.index(k) for k in shared]
    aligned_diff = float(np.nanmax(np.abs(coords[ia] - old_fit["coords"][ib])))
    aligned_same_gaps = bool(np.array_equal(np.isnan(coords[ia]), np.isnan(old_fit["coords"][ib])))

    # S04 distances
    caches = aligned_caches(fits, kept_new)
    kept, matrix, overlap, excluded = compute_distances(kept_new, caches, SETTINGS["min_pairwise_common_ca"])
    if np.isnan(matrix).any():
        raise ValueError("NaN pairs remain after pruning")
    np.save(OUT / "rmsd_recomputed.npy", matrix)
    np.save(GATE / "overlap.npy", overlap)
    pd.DataFrame({"row": range(len(kept)), "chain_key": kept, "entry_id": [k.split("|")[0] for k in kept]}).to_csv(
        OUT / "rmsd_recomputed_index.csv", index=False)

    old = np.load(ROOT / "results/rmsd_entry.npy")
    pos_old = {k: i for i, k in enumerate(frozen_keys)}
    both = [k for k in kept if k in pos_old]
    rn = [kept.index(k) for k in both]
    ro = [pos_old[k] for k in both]
    dmat = np.abs(matrix[np.ix_(rn, rn)] - old[np.ix_(ro, ro)])
    max_d_rmsd = float(dmat.max())
    n_pairs_diff = int((dmat[np.triu_indices(len(both), 1)] > 1e-6).sum())

    # S08 clusters
    primary = SETTINGS["primary_structural_threshold_A"]
    assign, z = structural_groups(kept, matrix, CUTS, primary)
    cl = pd.DataFrame({"chain_key": kept, "entry_id": [k.split("|")[0] for k in kept]})
    for cut in CUTS:
        cl[f"cut_{cut:.1f}"] = assign[f"cluster_{cut:.1f}A"].to_numpy()
    cl.to_csv(OUT / "clusters_recomputed.csv", index=False)
    fz = read(ROOT / "results/entry_clusters.csv")
    jj = cl.merge(fz, on="chain_key", suffixes=("_new", "_frozen"))
    ari = {f"{cut:.1f}": float(adjusted_rand_score(jj[f"cut_{cut:.1f}_frozen"], jj[f"cut_{cut:.1f}_new"])) for cut in CUTS}
    labels_identical = {f"{cut:.1f}": bool((jj[f"cut_{cut:.1f}_frozen"] == jj[f"cut_{cut:.1f}_new"]).all()) for cut in CUTS}
    jj["cluster_differs_1.5"] = cluster_mismatch(jj["cut_1.5_new"], jj["cut_1.5_frozen"]).to_numpy()

    # differing chains
    man = read(ROOT / "data/sv/manifest.csv").set_index("entry")
    diff_keys = set(cdiff.loc[cdiff.status != "identical", "chain_key"])
    diff_keys |= set(jj.loc[jj["cluster_differs_1.5"], "chain_key"])
    diff_keys |= set(frozen_keys) ^ set(kept)
    gc = pd.DataFrame({"chain_key": pd.Series(sorted(diff_keys), dtype=object)})
    gc["entry_id"] = gc.chain_key.str.split("|").str[0]
    gc = gc.merge(cdiff[["chain_key", "status", "max_abs_dxyz_A", "ca_only_regenerated", "ca_only_frozen"]],
                  on="chain_key", how="left")
    gc = gc.merge(jj[["chain_key", "cut_1.5_frozen", "cut_1.5_new", "cluster_differs_1.5"]], on="chain_key", how="left")
    gc["in_frozen_matrix"] = gc.chain_key.isin(frozen_keys)
    gc["in_recomputed_matrix"] = gc.chain_key.isin(kept)
    gc["latest_revision_date"] = gc.entry_id.map(man.latest_revision_date)
    gc["sha_matches_frozen"] = gc.entry_id.map(man.sha_matches_frozen)
    gc["revised_after_original"] = gc.entry_id.map(man.revised_after_original).astype(str).str.lower().eq("true")
    gc.to_csv(OUT / "gate_chains.csv", index=False)

    if ari["1.5"] == 1.0:
        verdict = "pass"
    elif len(gc) <= MAX_DIFFER and gc.revised_after_original.all():
        verdict = "tolerated"
    else:
        verdict = "blocked_for_coordinates"
    remediated = sorted(gc.chain_key) if verdict == "tolerated" else []
    res = {"verdict": verdict, "ari_by_cut": ari, "labels_identical_by_cut": labels_identical,
           "max_abs_delta_rmsd_A": max_d_rmsd, "n_pairs_delta_gt_1e-6": n_pairs_diff,
           "n_chains_recomputed": len(kept), "n_chains_frozen": len(frozen_keys), "n_chains_compared": len(both),
           "n_distance_exclusions": len(excluded), "reference_chain": reference_key,
           "core": S03.CORE_CONFIG["name"], "aligned_coords_max_abs_diff_A": aligned_diff,
           "aligned_coords_same_gaps": aligned_same_gaps,
           "n_chains_coordinates_differ": int((cdiff.status != "identical").sum()),
           "n_chains_cluster_differs_1.5": int(jj["cluster_differs_1.5"].sum()),
           "n_differing_chains": len(gc), "differing_chains": sorted(gc.chain_key),
           "remediated_chains": remediated, "max_differing_tolerated": MAX_DIFFER,
           "original_download_date": str(man.original_download_date.iloc[0]),
           "clusters_per_cut_recomputed": {f"{c:.1f}": int(cl[f"cut_{c:.1f}"].nunique()) for c in CUTS},
           "clusters_per_cut_frozen": {f"{c:.1f}": int(fz[f"cut_{c:.1f}"].nunique()) for c in CUTS}}
    (OUT / "gate.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k not in ("differing_chains", "remediated_chains")}, indent=1))
    print(counts.to_string(index=False))


if __name__ == "__main__":
    main()
