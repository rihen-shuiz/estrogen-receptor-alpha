"""SV05 (RUN_SV1028): stable core v4 from complete helices, compared with core v3.

Step text (docs/sv/briefs/SV05.md), applied as written:
1. Candidates = positions inside the SV02 global consensus helices (results/sv/SV02/consensus_global.csv) observed in
   >= 98 % of the 788 kept chains [DEFAULT] (SV00 results/sv/SV00/observed_mask.npy + _rows/_cols), minus 528-550.
2. Iterative superposition of all chains on the candidates, starting from the core v3 frame (SV00 gate fit, S03 code
   path, data/sv/gate/aligned_coordinates.npz via scripts/sv/common/core_frame.load_ca_core). One iteration:
     mean structure mu = nan-mean over chains of the candidate C-alpha;
     every chain is Kabsch-fitted (esr_analysis.align_structures.kabsch_align) on its observed candidates onto mu and the
     transform is applied to all its C-alpha;
     the ensemble is then re-anchored: one Kabsch fit of the new candidate mean onto the starting (v3-frame) mean,
     applied to every chain, so the frame cannot drift (a common rigid motion leaves every spread unchanged).
   Mean shift = mean over all observed (chain, candidate) C-alpha of |x_new - x_old| (A); stop when < 0.01 A
   (cap MAX_ITER iterations). Spread per position = sqrt(mean over observing chains of |x - mu|^2), the Y10 'sd'
   formula, computed for all 246 positions in the converged frame (selection uses candidates only).
3. [DEFAULT] core v4 = candidates with spread <= 0.5 A; < 80 positions -> BLOCKED. Also the 123 lowest-spread
   candidates (the size of v3), ties broken by position.
4. Compare with v3 (esr_analysis/resources/config/alignment_core_v3.json core_ranges): Jaccard and per-helix counts.
5. Recompute with v4 through the code path of results/core_sensitivity.csv (scripts/core_sensitivity.py, tracker K09):
   S03 scripts/align_selected_chains.fit_all (fit_to_reference, reference 7NFB|B|1, min core from settings) on the frozen
   residue mapping, S04 esr_analysis.calculate_rmsd.compute_distances on esr_analysis._task_io.aligned_caches
   (min_pairwise_common_ca), S08 esr_analysis.cluster_structures.structural_groups (average linkage, 1.0/1.5/2.0/2.5 A);
   ARI against results/entry_clusters.csv (S08), Spearman of all pairwise RMSDs against results/rmsd_entry.npy, cluster
   sizes. The same fit_all call with core v3 is checked against results/alignment/aligned_coordinates.npz first.
   Y10-style displacement profile (scripts/ext_displacement_profile.py formulas, re-implemented here so nothing is
   written to results/ext or figures/ext): S08 1.5 A clusters cl7/cl4, complete-H12 medoids (results/chain_status.csv
   h12_label == complete, smallest summed RMSD within the cluster) on the v4 matrix, v4-aligned C-alpha. The same
   function on the v3 inputs is checked against results/ext/cluster_displacement_profile.csv.
Nothing existing is overwritten; the superposition core of SV03/SV04/SV07 stays v3 (run file section 2).
Outputs: results/sv/SV05/{core_v4.json, spread.csv, ari.csv, rmsd_v4.npy, rmsd_v4_index.csv, clusters_v4.csv,
cluster_sizes.csv, core_comparison.csv, core_sensitivity_v4.csv, displacement_profile_v4.csv, iterations.csv,
numbers_core_v4.json}, data/sv/SV05/{aligned_coordinates_v4.npz, alignment_transforms_v4.json, global_alignment_v4.csv,
iterative_fit_coords.npz}.
Run: .venv/bin/python -B -I scripts/sv/SV05_core_v4.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts/sv"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402
from sklearn.metrics import adjusted_rand_score  # noqa: E402

import align_selected_chains as S03  # noqa: E402
from common.core_frame import load_ca_core  # noqa: E402
from esr_analysis._task_io import aligned_caches, coordinate_mapping, json_write  # noqa: E402
from esr_analysis.align_structures import apply_transform, kabsch_align  # noqa: E402
from esr_analysis.calculate_rmsd import compute_distances  # noqa: E402
from esr_analysis.cluster_structures import structural_groups  # noqa: E402

OUT = ROOT / "results/sv/SV05"
DATA = ROOT / "data/sv/SV05"
SETTINGS = S03.SETTINGS
CORE_V3 = json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())
CUTS = [1.0, 1.5, 2.0, 2.5]
OBS_MIN = 0.98            # [DEFAULT] step 1
EXCLUDE = (528, 550)      # step 1
SHIFT_TOL_A = 0.01        # step 2
MAX_ITER = 200            # cap (not reached = converged)
SPREAD_MAX_A = 0.5        # [DEFAULT] step 3
MIN_CORE_V4 = 80          # [DEFAULT] step 3: fewer -> BLOCKED
N_LIKE_V3 = 123           # step 3: size of v3


def ranges_of(pos):
    pos = sorted(int(p) for p in pos)
    out = []
    for p in pos:
        if out and p == out[-1][1] + 1:
            out[-1][1] = p
        else:
            out.append([p, p])
    return out


def helix_labels(positions):
    """Per position: SV02 consensus helix (global) if inside it, and the helix-table range name (SV02 table_start/end)."""
    hb = pd.read_csv(ROOT / "results/sv/SV02/helix_boundaries.csv")
    cons, table = [], []
    for p in positions:
        c = hb[(hb.global_start <= p) & (hb.global_end >= p)].helix_name.tolist()
        t = hb[(hb.table_start <= p) & (hb.table_end >= p)].helix_name.tolist()
        cons.append(c[0] if c else "")
        table.append(t[0] if t else "")
    return cons, table


def iterative_fit(X0, cand_cols):
    """Generalised superposition on the candidate columns, start = v3 frame; returns coords, log, mean-structure."""
    X = X0.copy()
    C = np.asarray(cand_cols)
    mu0 = np.nanmean(X[:, C], axis=0)
    log = []
    for it in range(1, MAX_ITER + 1):
        old = X.copy()
        mu = np.nanmean(X[:, C], axis=0)
        rms_to_mean = []
        for i in range(len(X)):
            sub = X[i, C]
            ok = np.isfinite(sub[:, 0])
            _, r, cm, cr = kabsch_align(sub[ok], mu[ok])
            fin = np.isfinite(X[i, :, 0])
            X[i, fin] = apply_transform(X[i, fin], r, cm, cr)
            rms_to_mean.append(float(np.sqrt(np.mean(np.sum((X[i, C][ok] - mu[ok]) ** 2, axis=1)))))
        # re-anchor the ensemble on the starting v3-frame mean (one common rigid motion)
        mu_new = np.nanmean(X[:, C], axis=0)
        _, r, cm, cr = kabsch_align(mu_new, mu0)
        fin = np.isfinite(X[..., 0])
        X[fin] = apply_transform(X[fin], r, cm, cr)
        d = np.linalg.norm(X[:, C] - old[:, C], axis=2)
        shift = float(np.nanmean(d))
        log.append({"iteration": it, "mean_shift_A": shift, "max_shift_A": float(np.nanmax(d)),
                    "mean_rms_to_mean_A": float(np.mean(rms_to_mean))})
        if shift < SHIFT_TOL_A:
            break
    return X, pd.DataFrame(log)


def spread(X):
    mu = np.nanmean(X, axis=0)
    dev2 = np.nansum((X - mu) ** 2, axis=2)
    dev2[~np.isfinite(X[..., 0])] = np.nan
    return np.sqrt(np.nanmean(dev2, axis=0))


def complete_medoid(d, cluster, h12_label, c):
    m = np.where((cluster == c) & (h12_label == "complete"))[0]
    return m[np.argmin(d[np.ix_(m, m)].sum(1))]


def y10_profile(x, pos, d, cluster, h12_label, keys):
    """scripts/ext_displacement_profile.py (Y10) formulas, returned as a table."""
    out = pd.DataFrame({"position": pos})
    mu, sd, ob = {}, {}, {}
    for c in (7, 4):
        xc = x[cluster == c]
        mu[c] = np.nanmean(xc, 0)
        dev = np.sqrt(np.nansum((xc - mu[c]) ** 2, 2))
        dev[~np.isfinite(xc[..., 0])] = np.nan
        sd[c] = np.sqrt(np.nanmean(dev ** 2, 0))
        ob[c] = np.isfinite(xc[..., 0]).mean(0)
    m7, m4 = complete_medoid(d, cluster, h12_label, 7), complete_medoid(d, cluster, h12_label, 4)
    out["d_mean_A"] = np.linalg.norm(mu[7] - mu[4], axis=1)
    out["d_medoid_A"] = np.linalg.norm(x[m7] - x[m4], axis=1)
    out["sd_cl7_A"], out["sd_cl4_A"] = sd[7], sd[4]
    out["obs_cl7"], out["obs_cl4"] = ob[7], ob[4]
    sq = out.d_mean_A ** 2
    out["share_sq"] = sq / sq.sum()
    out["medoid_cl7"], out["medoid_cl4"] = keys[m7], keys[m4]
    return out


def profile_numbers(p):
    h12 = (p.position >= SETTINGS["h12_range"][0]) & (p.position <= SETTINGS["h12_range"][1])
    tail = p.position >= 528
    top = p.sort_values("d_mean_A", ascending=False).head(5)
    return {"share_528_550": float(p.share_sq[tail].sum()), "share_h12_538_547": float(p.share_sq[h12].sum()),
            "max_d_mean_outside_528_A": float(p.d_mean_A[~tail].max()),
            "position_max_d_mean_outside_528": int(p.position[~tail][p.d_mean_A[~tail].idxmax()]),
            "top_positions": [int(v) for v in top.position], "top_d_mean_A": [round(float(v), 3) for v in top.d_mean_A],
            "medoid_cl7": str(p.medoid_cl7.iloc[0]), "medoid_cl4": str(p.medoid_cl4.iloc[0])}


def best_match_changed(new, old):
    tab = pd.crosstab(new, old)
    best = tab.idxmax(axis=1)
    return int((pd.Series(new).map(best).values != pd.Series(old).values).sum())


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    numbers = {}

    # ---------------- step 1: candidates
    mask = np.load(ROOT / "results/sv/SV00/observed_mask.npy")
    rows = pd.read_csv(ROOT / "results/sv/SV00/observed_mask_rows.csv")
    pos = pd.read_csv(ROOT / "results/sv/SV00/observed_mask_cols.csv").esr1_position.to_numpy()
    keys_g, pos_g, X0 = load_ca_core()
    if list(keys_g) != rows.chain_key.tolist() or not np.array_equal(pos_g, pos):
        raise ValueError("SV00 mask and gate coordinates are not in the same order")
    if not np.array_equal(np.isfinite(X0[..., 0]), mask.astype(bool)):
        raise ValueError("SV00 mask disagrees with the gate coordinates")
    n_chains = len(rows)
    obs_frac = mask.mean(axis=0)
    cons = pd.read_csv(ROOT / "results/sv/SV02/consensus_global.csv")
    in_cons = np.zeros(len(pos), bool)
    for r in cons.itertuples():
        in_cons |= (pos >= r.start) & (pos <= r.end)
    excluded = (pos >= EXCLUDE[0]) & (pos <= EXCLUDE[1])
    cand = in_cons & (obs_frac >= OBS_MIN) & ~excluded
    cand_cols = np.where(cand)[0]
    numbers["n_consensus_positions"] = int(in_cons.sum())
    numbers["n_consensus_below_98pct"] = int((in_cons & (obs_frac < OBS_MIN) & ~excluded).sum())
    numbers["consensus_below_98pct_positions"] = [int(p) for p in pos[in_cons & (obs_frac < OBS_MIN) & ~excluded]]
    numbers["n_consensus_in_528_550"] = int((in_cons & excluded).sum())
    numbers["n_candidates"] = int(cand.sum())

    # ---------------- step 2: iterative superposition
    X, log = iterative_fit(X0, cand_cols)
    log.to_csv(OUT / "iterations.csv", index=False)
    converged = bool(log.mean_shift_A.iloc[-1] < SHIFT_TOL_A)
    np.savez_compressed(DATA / "iterative_fit_coords.npz", chain_keys=np.asarray(keys_g), positions=pos, coords=X)
    sp = spread(X)
    sp0 = spread(X0)
    numbers.update(iterations=int(len(log)), converged=converged, final_mean_shift_A=float(log.mean_shift_A.iloc[-1]))

    # ---------------- step 3: core v4
    v3_pos = set(S03.positions(CORE_V3["core_ranges"]))
    cons_lab, table_lab = helix_labels(pos)
    tab = pd.DataFrame({"position": pos, "consensus_helix": cons_lab, "table_helix": table_lab,
                        "in_consensus_helix": in_cons, "n_observed": mask.sum(0), "obs_frac": obs_frac,
                        "excluded_528_550": excluded, "candidate": cand, "spread_A": sp, "spread_v3_frame_A": sp0})
    c = tab[tab.candidate].sort_values(["spread_A", "position"])
    tab["rank_among_candidates"] = np.nan
    tab.loc[c.index, "rank_among_candidates"] = np.arange(1, len(c) + 1)
    tab["in_core_v3"] = tab.position.isin(v3_pos)
    tab["in_core_v4"] = tab.candidate & (tab.spread_A <= SPREAD_MAX_A)
    tab["in_lowest_123"] = tab.candidate & (tab.rank_among_candidates <= N_LIKE_V3)
    tab.to_csv(OUT / "spread.csv", index=False)
    v4 = tab.position[tab.in_core_v4].astype(int).tolist()
    low123 = tab.position[tab.in_lowest_123].astype(int).tolist()
    blocked = len(v4) < MIN_CORE_V4
    numbers.update(n_core_v4=len(v4), blocked=blocked, n_lowest=len(low123),
                   spread_max_in_lowest_123_A=float(tab.spread_A[tab.in_lowest_123].max()),
                   candidate_spread_min_A=float(c.spread_A.min()), candidate_spread_median_A=float(c.spread_A.median()),
                   candidate_spread_max_A=float(c.spread_A.max()))

    # ---------------- step 4: compare with v3
    def jac(a, b):
        a, b = set(a), set(b)
        return len(a & b) / len(a | b)

    numbers["jaccard_v4_v3"] = jac(v4, v3_pos)
    numbers["jaccard_lowest123_v3"] = jac(low123, v3_pos)
    comp = []
    groups = [g for g in pd.read_csv(ROOT / "results/sv/SV02/helix_boundaries.csv").helix_name] + [""]
    for g in groups:
        sel = tab[tab.table_helix == g]
        s3, s4, s1 = set(sel.position[sel.in_core_v3]), set(sel.position[sel.in_core_v4]), set(sel.position[sel.in_lowest_123])
        if not (s3 or s4 or s1):
            continue
        comp.append({"helix": g or "outside helix-table ranges", "v3_n": len(s3), "v4_n": len(s4), "both_v3_v4": len(s3 & s4),
                     "v3_only": len(s3 - s4), "v4_only": len(s4 - s3),
                     "v3_only_positions": " ".join(f"{a}-{b}" if a != b else str(a) for a, b in ranges_of(s3 - s4)),
                     "v4_only_positions": " ".join(f"{a}-{b}" if a != b else str(a) for a, b in ranges_of(s4 - s3)),
                     "lowest123_n": len(s1), "lowest123_and_v3": len(s1 & s3),
                     "lowest123_only_positions": " ".join(f"{a}-{b}" if a != b else str(a) for a, b in ranges_of(s1 - s3)),
                     "v3_not_lowest123_positions": " ".join(f"{a}-{b}" if a != b else str(a) for a, b in ranges_of(s3 - s1)),
                     "v4_ranges": " ".join(f"{a}-{b}" if a != b else str(a) for a, b in ranges_of(s4))})
    comp = pd.DataFrame(comp)
    comp.to_csv(OUT / "core_comparison.csv", index=False)

    core = {"name": "alignment_core_v4", "status": "proposal from SV05; superposition core of SV03/SV04/SV07 stays v3",
            "core_ranges": ranges_of(v4), "positions": v4, "residues_total": len(v4),
            "helices": sorted({h for h, p in zip(tab.table_helix, tab.position) if p in v4}),
            "consensus_helices": sorted({h for h, p in zip(tab.consensus_helix, tab.position) if p in v4}),
            "rule": (f"positions in the SV02 global consensus helices observed in >= {OBS_MIN:.0%} of {n_chains} kept chains, "
                     f"excluding {EXCLUDE[0]}-{EXCLUDE[1]}; iterative superposition on them from the core v3 frame until the "
                     f"mean C-alpha shift < {SHIFT_TOL_A} A; core v4 = spread <= {SPREAD_MAX_A} A; BLOCKED if < {MIN_CORE_V4}"),
            "blocked": blocked, "n_candidates": int(cand.sum()), "iterations": int(len(log)), "converged": converged,
            "final_mean_shift_A": float(log.mean_shift_A.iloc[-1]),
            "lowest_spread_123": {"core_ranges": ranges_of(low123), "positions": low123, "residues_total": len(low123),
                                  "max_spread_A": numbers["spread_max_in_lowest_123_A"]},
            "jaccard_vs_v3": numbers["jaccard_v4_v3"], "jaccard_lowest123_vs_v3": numbers["jaccard_lowest123_v3"],
            "v3_core_ranges": CORE_V3["core_ranges"],
            "source": ("results/sv/SV02/consensus_global.csv, results/sv/SV00/observed_mask.npy, data/sv/gate/"
                       "aligned_coordinates.npz; spread table results/sv/SV05/spread.csv; script scripts/sv/SV05_core_v4.py"),
            "h12_settings_range": SETTINGS["h12_range"]}
    (OUT / "core_v4.json").write_text(json.dumps(core, indent=2) + "\n")
    if blocked:
        (OUT / "numbers_core_v4.json").write_text(json.dumps(numbers, indent=1) + "\n")
        print(f"BLOCKED: core v4 has {len(v4)} positions (< {MIN_CORE_V4})")
        return

    # ---------------- step 5: S03 -> S04 -> S08 with core v4
    keys = (ROOT / "results/alignment/rmsd_entry_index.txt").read_text().split()
    if S03.eligible_chain_keys() != keys or keys != list(keys_g):
        raise ValueError("chain order differs between S03, the matrix index and the SV00 gate")
    mapping = pd.read_csv(S03.MAPPING, dtype={"chain_key": str}, keep_default_na=False, na_values=[""], low_memory=False)
    reference_key = json.loads(S03.REFERENCE.read_text())["chain_key"]
    reference = coordinate_mapping(mapping, reference_key)
    mobiles = {k: coordinate_mapping(mapping, k) for k in keys}
    # path check: the same call with core v3 must reproduce results/alignment/aligned_coordinates.npz
    v3_coords, _, _ = S03.fit_all(mobiles, reference, reference_key, keys, CORE_V3["core_ranges"])
    frozen = np.load(ROOT / "results/alignment/aligned_coordinates.npz")
    path_check = {"max_abs_diff_A": float(np.nanmax(np.abs(v3_coords - frozen["coords"]))),
                  "same_gaps": bool(np.array_equal(np.isnan(v3_coords), np.isnan(frozen["coords"])))}
    numbers["s03_path_check_v3"] = path_check
    v4_coords, v4_transforms, v4_summary = S03.fit_all(mobiles, reference, reference_key, keys, core["core_ranges"])
    canonical = np.arange(SETTINGS["esr1_lbd_start"], SETTINGS["esr1_lbd_end"] + 1)
    npz = DATA / "aligned_coordinates_v4.npz"
    np.savez_compressed(npz, chain_keys=np.asarray(keys), positions=canonical, coords=v4_coords)
    json_write(DATA / "alignment_transforms_v4.json", v4_transforms)
    v4_summary.to_csv(DATA / "global_alignment_v4.csv", index=False)
    numbers["v4_fit_status"] = v4_summary.status.value_counts().to_dict()
    numbers["v4_core_rmsd_median_A"] = float(v4_summary.core_alignment_rmsd_A.median())
    numbers["v4_core_rmsd_max_A"] = float(v4_summary.core_alignment_rmsd_A.max())
    numbers["v4_n_common_core_min"] = int(v4_summary.n_common_core_ca.min())

    caches = aligned_caches(npz, keys)
    kept, mat, _, excl = compute_distances(keys, caches, SETTINGS["min_pairwise_common_ca"])
    if kept != keys or np.isnan(mat).any():
        raise ValueError(f"v4 matrix does not cover the 788 chains ({len(kept)} kept, {len(excl)} excluded)")
    np.save(OUT / "rmsd_v4.npy", mat)
    pd.DataFrame({"row": range(len(keys)), "chain_key": keys, "entry_id": [k.split("|")[0] for k in keys]}).to_csv(
        OUT / "rmsd_v4_index.csv", index=False)
    primary = np.load(ROOT / "results/rmsd_entry.npy")
    s08 = pd.read_csv(ROOT / "results/entry_clusters.csv")
    if s08.chain_key.tolist() != keys:
        raise ValueError("entry_clusters.csv row order does not match the matrix")
    up = np.triu_indices(len(keys), k=1)
    diff = np.abs(mat[up] - primary[up])
    rho = float(spearmanr(primary[up], mat[up])[0])
    numbers.update(n_distance_exclusions=len(excl), spearman_pairwise_v4_vs_v3=rho,
                   median_abs_diff_A=float(np.median(diff)), max_abs_diff_A=float(diff.max()),
                   p95_abs_diff_A=float(np.percentile(diff, 95)),
                   n_pairs_abs_diff_gt_0_5A=int((diff > 0.5).sum()), n_pairs=int(len(diff)))

    assign, _ = structural_groups(keys, mat, CUTS, SETTINGS["primary_structural_threshold_A"])
    cl = assign.rename(columns={f"cluster_{c:.1f}A": f"cut_{c:.1f}" for c in CUTS})
    cl.insert(1, "entry_id", [k.split("|")[0] for k in keys])
    for c in CUTS:
        cl[f"s08_cut_{c:.1f}"] = s08[f"cut_{c:.1f}"].values
    cl.to_csv(OUT / "clusters_v4.csv", index=False)

    ari_rows, size_rows = [], []
    for c in CUTS:
        a, b = s08[f"cut_{c:.1f}"].to_numpy(), cl[f"cut_{c:.1f}"].to_numpy()
        sa, sb = pd.Series(a).value_counts(), pd.Series(b).value_counts()
        ari_rows.append({"cut_A": c, "n_chains": len(keys), "n_clusters_v3_s08": int(sa.size), "n_clusters_v4": int(sb.size),
                         "ari_v4_vs_s08": float(adjusted_rand_score(a, b)),
                         "chains_changed_best_match": best_match_changed(b, a),
                         "top3_fraction_v3": float(sa.head(3).sum() / len(keys)), "top3_fraction_v4": float(sb.head(3).sum() / len(keys)),
                         "sizes_v3_s08": ";".join(str(int(v)) for v in sa.values),
                         "sizes_v4": ";".join(str(int(v)) for v in sb.values)})
        ct = pd.crosstab(b, a)
        for lab, n in sb.items():
            best = ct.loc[lab].idxmax()
            size_rows.append({"cut_A": c, "core": "v4", "cluster": int(lab), "n_chains": int(n),
                              "best_match_other_core_cluster": int(best), "n_shared": int(ct.loc[lab, best]),
                              "size_of_best_match": int(sa[best])})
        for lab, n in sa.items():
            best = ct[lab].idxmax()
            size_rows.append({"cut_A": c, "core": "v3 (S08)", "cluster": int(lab), "n_chains": int(n),
                              "best_match_other_core_cluster": int(best), "n_shared": int(ct.loc[best, lab]),
                              "size_of_best_match": int(sb[best])})
    ari = pd.DataFrame(ari_rows)
    ari["spearman_pairwise_v4_vs_v3"] = rho
    ari["median_abs_diff_A"] = numbers["median_abs_diff_A"]
    ari["max_abs_diff_A"] = numbers["max_abs_diff_A"]
    ari.to_csv(OUT / "ari.csv", index=False)
    pd.DataFrame(size_rows).to_csv(OUT / "cluster_sizes.csv", index=False)

    # K09-format table (results/core_sensitivity.csv columns), v3 primary + v4; the K09 file itself is not touched
    sens = []
    for name, size, labels, sp_ in [(CORE_V3["name"] + " (primary)", len(v3_pos), {c: s08[f"cut_{c:.1f}"] for c in CUTS}, 1.0),
                                    ("alignment_core_v4", len(v4), {c: cl[f"cut_{c:.1f}"] for c in CUTS}, rho)]:
        row = {"core": name, "core_residues": size, "n_chains": len(keys), "spearman_pairwise": sp_,
               "median_abs_diff_A": numbers["median_abs_diff_A"] if "v4" in name else 0.0,
               "max_abs_diff_A": numbers["max_abs_diff_A"] if "v4" in name else 0.0}
        for c in CUTS:
            tag = f"{c:.1f}".replace(".", "_")
            row[f"n_clusters_{tag}"] = int(pd.Series(labels[c]).nunique())
            row[f"ari_{tag}"] = float(adjusted_rand_score(s08[f"cut_{c:.1f}"], labels[c]))
        row["ari_vs_primary"] = row["ari_1_5"]
        sens.append(row)
    pd.DataFrame(sens).to_csv(OUT / "core_sensitivity_v4.csv", index=False)

    # ---------------- step 5: Y10-style displacement profile on v4
    status = pd.read_csv(ROOT / "results/chain_status.csv").set_index("chain_key").loc[keys]
    h12_label = status.h12_label.to_numpy()
    cluster = s08["cut_1.5"].to_numpy().astype(int)
    karr = np.asarray(keys)
    y10 = pd.read_csv(ROOT / "results/ext/cluster_displacement_profile.csv")
    # check: the function on the v3 inputs reproduces the Y10 file
    p3 = y10_profile(frozen["coords"].astype(float), canonical, primary, cluster, h12_label, karr)
    numcols = ["d_mean_A", "d_medoid_A", "sd_cl7_A", "sd_cl4_A", "obs_cl7", "obs_cl4", "share_sq"]
    numbers["y10_reproduction_max_abs_diff"] = float(np.nanmax(np.abs(p3[numcols].to_numpy() - y10[numcols].to_numpy())))
    numbers["y10_reproduction_same_medoids"] = bool(p3.medoid_cl7.iloc[0] == y10.medoid_cl7.iloc[0]
                                                    and p3.medoid_cl4.iloc[0] == y10.medoid_cl4.iloc[0])
    p4 = y10_profile(v4_coords, canonical, mat, cluster, h12_label, karr)
    prof = p4.copy()
    for col in ["d_mean_A", "d_medoid_A", "sd_cl7_A", "sd_cl4_A", "share_sq"]:
        prof[f"{col}_v3_y10"] = y10[col].to_numpy()
    prof["delta_d_mean_A"] = prof.d_mean_A - prof.d_mean_A_v3_y10
    prof["medoid_cl7_v3_y10"], prof["medoid_cl4_v3_y10"] = y10.medoid_cl7.iloc[0], y10.medoid_cl4.iloc[0]
    prof["in_core_v3"] = prof.position.isin(v3_pos)
    prof["in_core_v4"] = prof.position.isin(v4)
    prof.to_csv(OUT / "displacement_profile_v4.csv", index=False)
    numbers["profile_v4"] = profile_numbers(p4)
    numbers["profile_v3_y10"] = profile_numbers(y10)
    numbers["profile_max_abs_delta_d_mean_A"] = float(np.nanmax(np.abs(prof.delta_d_mean_A)))
    numbers["profile_position_max_abs_delta"] = int(prof.position[prof.delta_d_mean_A.abs().idxmax()])
    numbers["profile_max_abs_delta_outside_528_A"] = float(np.nanmax(np.abs(prof.delta_d_mean_A[prof.position < 528])))
    (OUT / "numbers_core_v4.json").write_text(json.dumps(numbers, indent=1, default=lambda o: o.item()) + "\n")
    print(json.dumps(numbers, indent=1, default=lambda o: o.item()))


if __name__ == "__main__":
    main()
