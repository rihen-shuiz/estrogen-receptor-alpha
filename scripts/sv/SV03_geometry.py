"""SV03 steps 1-3 (RUN_SV1028): helix lengths with censoring, helix axes and per-chain geometry descriptors.

Inputs (tracker cards SV00 / SV02 `out`, rule 1.10):
- results/sv/SV02/helix_boundaries.csv: global consensus helices (H1 312-321, H3 342-361, H4/5 372-394, H7 422-437,
  H8 442-455, H9 473-492, H10/11 497-526, H12 538-544); helix names from SV02 (ext_common.HELIX_NAMES).
- data/sv/SV02/dssp_residues.csv.gz: pydssp state per residue of the full chain with its pydssp segment (chain break
  = new segment), the residue table behind results/sv/SV02/sse_codes_per_chain.csv and sse_per_chain.csv.
- C-alpha in the core v3 frame: data/sv/gate/aligned_coordinates.npz with the SV00 gate transforms
  data/sv/gate/alignment_transforms.json (S03 code path, core v3, reference 7NFB|B|1); deposited backbone
  data/sv/backbone/backbone.npz for the frame check.
- results/sv/SV00/chain_table.csv: s08_cluster, ligand_chemotype, ligand heavy atoms, basic amine, carboxylic acid,
  peptide_bound, construct_mutations.

Rules fixed before any output was looked at:
- Helix ranges = SV02 global consensus helices, except H12: the run file (section 2) keeps the H12 window 538-547
  for all primary analyses, so H12 = 538-547; the SV02 data-derived H12 (538-544, a proposal) is kept as a
  sensitivity copy of the H12 columns (suffix _w544), not tested.
- "H11" in the descriptors = the SV02 consensus helix H10/11 (497-526), which SV02 keeps as one helix. Its two
  halves (step 2) split the consensus range at its midpoint: H11N 497-511, H11C 512-526.
- Step 1, length per chain x helix: chain helix runs = maximal runs of consecutive UniProt positions with pydssp 'H'
  inside one pydssp segment (as SV02). Runs overlapping the helix range: start = first position of the first run,
  end = last position of the last run, length = end - start + 1 (n_helical = H residues in those runs).
  Terminus state (start: residues start-1, start-2; end: end+1, end+2):
    censored = one of the two residues is not observed (no pydssp state) or lies across a chain break;
    real     = both observed in the same segment and not helical (the step's 'real end');
    gap      = both observed, the next one not helical but the second one helical (a 1-residue break before
               another H run): neither censored nor a real end.
  Uncensored length = both termini 'real' (only these lengths are compared); censoring fraction per helix = chains
  with >= 1 censored terminus / chains with a run in the range.
- Step 2 [DEFAULT] axis: scripts/sv/common/helix_geometry.fit_axis on the C-alpha (core v3 frame) of the positions
  in the helix range that are observed and pydssp 'H' in the chain; >= 6 such C-alpha, >= 2 consecutive 4-residue
  windows; helix centroid = mean of the window centroids (on the axis), reported in the core v3 frame.
- Step 3 descriptors: angle H12-H3, H12-H11, H12-H4/5, H11-H3 (oriented N -> C axes, 0-180 deg); distance from the
  H12 centroid to the H3 axis and to the H11 axis (A); H11 kink = angle H11N-H11C; H11-H12 hinge angle =
  angle H11C-H12 (the H11 half next to the hinge loop; the whole-H11 angle is angle H12-H11).
Outputs: results/sv/SV03/helix_per_chain.csv (long: chain x helix incl. H11 halves and H12 538-544),
  geometry_per_chain.csv (wide, one row per kept chain), helix_length_summary.csv, h12_window_sensitivity.csv
  (H12 descriptors 538-547 vs 538-544, descriptive), axis_method_check.csv (axis rule on ideal helices of 6-12
  C-alpha), frame_check.json, numbers_geometry.json; data/sv/SV03/core_v3_transforms.csv (frozen copy of the frame for SV04 / SV13).
Run: .venv/bin/python -B -I scripts/sv/SV03_geometry.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from common import core_frame as cf  # noqa: E402
from common.helix_geometry import axis_angle, fit_axis, point_axis_distance  # noqa: E402

OUT = ROOT / "results/sv/SV03"
DATA = ROOT / "data/sv/SV03"
H12_PRIMARY = (538, 547)
H12_SENS_NAME = "H12_w544"


def sha16(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16]


def helix_ranges():
    hb = pd.read_csv(ROOT / "results/sv/SV02/helix_boundaries.csv")
    rng = {r.helix_name: (int(r.global_start), int(r.global_end)) for r in hb.itertuples()}
    h12_sv02 = rng["H12"]
    rng["H12"] = H12_PRIMARY
    s, e = rng["H10/11"]
    half = (s + e) // 2   # 497-526 -> 511
    halves = {"H11N": (s, half), "H11C": (half + 1, e)}
    return rng, halves, h12_sv02


def frame(ct_keys):
    """Freeze the gate transforms, check the convention against the gate's aligned C-alpha."""
    DATA.mkdir(parents=True, exist_ok=True)
    tr = cf.transforms_from_gate()
    tr.to_csv(cf.FROZEN, index=False)
    T = cf.load_transforms()
    bb = np.load(ROOT / "data/sv/backbone/backbone.npz")
    keys, pos, ca_core = cf.load_ca_core()
    assert (bb["chain_keys"].astype(str) == keys).all() and (bb["positions"] == pos).all()
    ca_dep = bb["coords"][:, :, 1, :]
    maxd, mask_same = 0.0, True
    for i, k in enumerate(keys):
        x = cf.to_core(ca_dep[i], T[k])
        ok = ~np.isnan(ca_core[i, :, 0])
        mask_same &= bool((np.isnan(ca_dep[i, :, 0]) == ~ok).all())
        maxd = max(maxd, float(np.abs(x[ok] - ca_core[i][ok]).max()))
    chk = {"source": {"transforms": "data/sv/gate/alignment_transforms.json",
                      "aligned_ca": "data/sv/gate/aligned_coordinates.npz",
                      "backbone": "data/sv/backbone/backbone.npz"},
           "sha256_16": {p: sha16(ROOT / p) for p in ("data/sv/gate/alignment_transforms.json",
                                                      "data/sv/gate/aligned_coordinates.npz",
                                                      "data/sv/backbone/backbone.npz",
                                                      "data/sv/SV03/core_v3_transforms.csv")},
           "core": "alignment_core_v3 (338-362, 371-396, 420-439, 441-456, 496-531), reference 7NFB|B|1 (SV00 gate)",
           "convention": "x_core = (x_deposited - mobile_centroid) @ R + reference_centroid (row vectors)",
           "n_chains": int(len(keys)), "status_counts": tr.status.value_counts().to_dict(),
           "max_abs_diff_backbone_CA_vs_aligned_A": maxd, "observed_mask_identical": mask_same,
           "chain_order_equals_chain_table": bool(list(keys) == list(ct_keys)),
           "frozen_copy": "data/sv/SV03/core_v3_transforms.csv",
           "reuse": "scripts/sv/common/core_frame.py: load_transforms(), to_core(), load_ca_core()"}
    (OUT / "frame_check.json").write_text(json.dumps(chk, indent=1))
    return keys, pos, ca_core, chk


def chain_residues(res):
    """pos -> (segment, sse) per chain, and H runs inside one segment (SV02 rule)."""
    obs, runs = {}, {}
    for key, g in res.groupby("chain_key", sort=False):
        g = g[g.esr1_position.notna()].sort_values("esr1_position")
        p = g.esr1_position.astype(int).to_numpy()
        seg = g.segment.to_numpy()
        h = (g.sse == "H").to_numpy()
        obs[key] = {int(a): (int(b), c) for a, b, c in zip(p, seg, g.sse.to_numpy())}
        rr, start = [], None
        for i in range(len(p)):
            cont = i > 0 and p[i] == p[i - 1] + 1 and seg[i] == seg[i - 1]
            if start is not None and (not h[i] or not cont):
                rr.append((start, int(p[i - 1])))
                start = None
            if h[i] and start is None:
                start = int(p[i])
        if start is not None:
            rr.append((start, int(p[-1])))
        runs[key] = rr
    return obs, runs


def terminus(o, p, nbs):
    seg = o[p][0]
    if any(nb not in o or o[nb][0] != seg for nb in nbs):
        return "censored"
    if all(o[nb][1] != "H" for nb in nbs):
        return "real"
    return "gap"


def helix_row(key, o, rr, name, s, e, ca_i, pos_index, with_length=True):
    ov = [r for r in rr if r[1] >= s and r[0] <= e]
    n_obs = sum(1 for p in range(s, e + 1) if p in o)
    rec = {"helix": name, "range_start": s, "range_end": e, "n_observed_in_range": n_obs, "n_runs": len(ov)}
    if with_length:
        if ov:
            a, b = ov[0][0], ov[-1][1]
            st, et = terminus(o, a, (a - 1, a - 2)), terminus(o, b, (b + 1, b + 2))
            rec.update(status="helix", start=a, end=b, length=b - a + 1,
                       n_helical=int(sum(y - x + 1 for x, y in ov)), start_state=st, end_state=et,
                       uncensored=(st == "real" and et == "real"),
                       censored=(st == "censored" or et == "censored"))
        else:
            rec.update(status="not observed" if n_obs == 0 else "not helical", start=np.nan, end=np.nan,
                       length=np.nan, n_helical=0, start_state="", end_state="", uncensored=False, censored=False)
    # axis: observed helical C-alpha inside the range
    hp = [p for p in range(s, e + 1) if p in o and o[p][1] == "H" and p in pos_index
          and not np.isnan(ca_i[pos_index[p], 0])]
    xyz = np.array([ca_i[pos_index[p]] for p in hp]).reshape(-1, 3)
    ax = fit_axis(hp, xyz)
    rec.update(n_ca_axis=ax["n_ca"], n_windows=ax["n_windows"], axis_status=ax["status"],
               centroid_x=ax["centre"][0], centroid_y=ax["centre"][1], centroid_z=ax["centre"][2],
               axis_ux=ax["direction"][0], axis_uy=ax["direction"][1], axis_uz=ax["direction"][2],
               rms_window_offset_A=ax["rms_window_offset"], axis_span_A=ax["span_A"])
    return rec, ax


def descriptors(ax, h12="H12"):
    a = {k: v for k, v in ax.items()}

    def ok(n):
        return n in a and a[n]["status"] == "ok"

    def ang(m, n):
        return axis_angle(a[m]["direction"], a[n]["direction"]) if ok(m) and ok(n) else np.nan

    def dist(m, n):
        return point_axis_distance(a[m]["centre"], a[n]["centre"], a[n]["direction"]) if ok(m) and ok(n) else np.nan

    return {"angle_H12_H3_deg": ang(h12, "H3"), "angle_H12_H11_deg": ang(h12, "H10/11"),
            "angle_H12_H45_deg": ang(h12, "H4/5"), "angle_H11_H3_deg": ang("H10/11", "H3"),
            "dist_H12c_H3axis_A": dist(h12, "H3"), "dist_H12c_H11axis_A": dist(h12, "H10/11"),
            "kink_H11_deg": ang("H11N", "H11C"), "hinge_H11C_H12_deg": ang("H11C", h12)}


def tag(name):
    return name.replace("/", "")


def axis_method_check():
    """Error of the [DEFAULT] axis rule on ideal alpha-helices (radius 2.3 A, rise 1.5 A, 100 deg per residue) with
    n = 6-12 C-alpha: angle to the true axis and distance of the fitted centre from it. The error does not depend on
    the helix phase (a rotation about the true axis), so one helix per n is enough."""
    rows = []
    for n in range(6, 13):
        k = np.arange(n)
        th = np.radians(100.0 * k)
        xyz = np.c_[2.3 * np.cos(th), 2.3 * np.sin(th), 1.5 * k]
        a = fit_axis(np.arange(n), xyz)
        rows.append({"n_ca": n, "n_windows": a["n_windows"], "angle_error_deg": axis_angle(a["direction"], [0, 0, 1]),
                     "centre_offset_A": point_axis_distance(a["centre"], [0, 0, 0], [0, 0, 1]),
                     "rms_window_offset_A": a["rms_window_offset"]})
    pd.DataFrame(rows).to_csv(OUT / "axis_method_check.csv", index=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    axis_method_check()
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", keep_default_na=False, na_values=[""])
    rng, halves, h12_sv02 = helix_ranges()
    keys, pos, ca_core, chk = frame(ct.chain_key)
    pos_index = {int(p): i for i, p in enumerate(pos)}
    kidx = {k: i for i, k in enumerate(keys)}
    res = pd.read_csv(ROOT / "data/sv/SV02/dssp_residues.csv.gz")
    obs, runs = chain_residues(res)

    long_rows, wide_rows = [], []
    group = ct.s08_cluster.map({7: "cl7", 4: "cl4"}).fillna("minor")
    muts = ct.construct_mutations.fillna("none").astype(str).str.split(";")
    for i, r in ct.iterrows():
        key = r.chain_key
        o, rr, ca_i = obs[key], runs[key], ca_core[kidx[key]]
        base = {"chain_key": key, "entry_id": r.entry_id, "s08_cluster": int(r.s08_cluster), "group": group[i]}
        ax = {}
        w = dict(base)
        for name, (s, e) in list(rng.items()) + [(H12_SENS_NAME, h12_sv02)]:
            rec, a = helix_row(key, o, rr, name, s, e, ca_i, pos_index)
            long_rows.append({**base, **rec})
            ax[name] = a
            t = tag(name)
            w.update({f"len_{t}": rec["length"], f"ulen_{t}": rec["length"] if rec["uncensored"] else np.nan,
                      f"start_{t}": rec["start"], f"end_{t}": rec["end"],
                      f"state_{t}": f"{rec['start_state']}/{rec['end_state']}" if rec["status"] == "helix" else rec["status"],
                      f"axis_{t}": rec["axis_status"] == "ok",
                      f"centroid_{t}_x": rec["centroid_x"], f"centroid_{t}_y": rec["centroid_y"],
                      f"centroid_{t}_z": rec["centroid_z"]})
        for name, (s, e) in halves.items():
            rec, a = helix_row(key, o, rr, name, s, e, ca_i, pos_index, with_length=False)
            long_rows.append({**base, **rec, "status": "half of H10/11 (axis only)"})
            ax[name] = a
            w[f"axis_{name}"] = rec["axis_status"] == "ok"
        d = descriptors(ax)
        ds = descriptors({**ax, "H12": ax[H12_SENS_NAME]})
        w.update(d)
        w.update({f"{k}_w544": v for k, v in ds.items() if "H12" in k})
        m = muts[i]
        w.update({"ligand_chemotype": r.ligand_chemotype, "ligand_id_primary": r.ligand_id_primary,
                  "ligand_heavy_atoms": r.ligand_heavy_atoms, "ligand_basic_amine": r.ligand_basic_amine,
                  "ligand_carboxylic_acid": r.ligand_carboxylic_acid, "peptide_bound": r.peptide_bound,
                  "construct_mutations": r.construct_mutations,
                  "mut_Y537S": "Y537S" in m, "mut_L536S": "L536S" in m, "mut_D538G": "D538G" in m,
                  "mut_other": any(x not in ("Y537S", "L536S", "D538G", "none", "") for x in m)})
        wide_rows.append(w)

    lg = pd.DataFrame(long_rows)
    lg.to_csv(OUT / "helix_per_chain.csv", index=False)
    wd = pd.DataFrame(wide_rows)
    wd.to_csv(OUT / "geometry_per_chain.csv", index=False)

    # censoring summary per helix (step 1)
    summ = []
    for name in list(rng) + [H12_SENS_NAME]:
        x = lg[lg.helix == name]
        hx = x[x.status == "helix"]
        rec = {"helix": name, "range": f"{int(x.range_start.iloc[0])}-{int(x.range_end.iloc[0])}",
               "n_chains": len(x), "n_with_run": len(hx), "n_not_observed": int((x.status == "not observed").sum()),
               "n_not_helical": int((x.status == "not helical").sum()),
               "start_censored": int((hx.start_state == "censored").sum()),
               "end_censored": int((hx.end_state == "censored").sum()),
               "start_gap": int((hx.start_state == "gap").sum()), "end_gap": int((hx.end_state == "gap").sum()),
               "n_censored": int(hx.censored.sum()), "n_uncensored": int(hx.uncensored.sum()),
               "n_gap_only": int(((~hx.censored.astype(bool)) & (~hx.uncensored.astype(bool))).sum()),
               "censoring_fraction": hx.censored.mean() if len(hx) else np.nan,
               "uncensored_fraction": hx.uncensored.mean() if len(hx) else np.nan,
               "n_axis": int((x.axis_status == "ok").sum())}
        for g in ("cl7", "cl4", "minor"):
            u = hx[(hx.group == g) & hx.uncensored].length
            rec[f"n_uncensored_{g}"] = len(u)
            rec[f"median_uncensored_length_{g}"] = u.median() if len(u) else np.nan
            rec[f"q1_{g}"] = u.quantile(0.25) if len(u) else np.nan
            rec[f"q3_{g}"] = u.quantile(0.75) if len(u) else np.nan
            rec[f"censoring_fraction_{g}"] = hx[hx.group == g].censored.mean() if (hx.group == g).any() else np.nan
        summ.append(rec)
    sm_ = pd.DataFrame(summ)
    sm_.to_csv(OUT / "helix_length_summary.csv", index=False)

    # cross-check with SV02 sse_per_chain.csv (same ranges: every helix except H12 538-547; H12_w544 = SV02 H12)
    sv02 = pd.read_csv(ROOT / "results/sv/SV02/sse_per_chain.csv")
    cmp = {}
    for name in list(rng) + [H12_SENS_NAME]:
        ref = sv02[sv02.helix_name == ("H12" if name == H12_SENS_NAME else name)].set_index("chain_key")
        if name == "H12":
            continue
        mine = lg[lg.helix == name].set_index("chain_key").loc[ref.index]
        hx = ref.status == "helix"
        cmp[name] = {"status_equal": int((mine.status == ref.status).sum()),
                     "start_end_equal": int(((mine.start[hx] == ref.start[hx]) & (mine.end[hx] == ref.end[hx])).sum()),
                     "start_censored_equal": int((mine.start_state[hx].eq("censored") ==
                                                  ref.start_censored[hx].astype(str).eq("True")).sum()),
                     "end_censored_equal": int((mine.end_state[hx].eq("censored") ==
                                                ref.end_censored[hx].astype(str).eq("True")).sum()),
                     "n_helix": int(hx.sum()), "n": int(len(ref))}
    # H12 window sensitivity (descriptive only): 538-547 (primary) vs SV02 data-derived 538-544
    sens = []
    for c in [c for c in wd.columns if c.endswith("_w544") and not c.startswith(("len_", "ulen_", "start_", "end_",
                                                                                     "state_", "axis_", "centroid_"))]:
        p, s_ = wd[c[:-5]], wd[c]
        ok = p.notna() & s_.notna()
        sens.append({"descriptor": c[:-5], "n_both": int(ok.sum()), "n_primary_only": int((p.notna() & ~ok).sum()),
                     "n_w544_only": int((s_.notna() & ~ok).sum()),
                     "spearman_rho": p[ok].corr(s_[ok], method="spearman"),
                     "median_abs_diff": (p[ok] - s_[ok]).abs().median(), "max_abs_diff": (p[ok] - s_[ok]).abs().max()})
    for c in ("len_H12", "ulen_H12"):
        p, s_ = wd[c], wd[f"{c}_w544"]
        ok = p.notna() & s_.notna()
        sens.append({"descriptor": c, "n_both": int(ok.sum()), "n_primary_only": int((p.notna() & ~ok).sum()),
                     "n_w544_only": int((s_.notna() & ~ok).sum()), "spearman_rho": p[ok].corr(s_[ok], method="spearman"),
                     "median_abs_diff": (p[ok] - s_[ok]).abs().median(), "max_abs_diff": (p[ok] - s_[ok]).abs().max()})
    pd.DataFrame(sens).to_csv(OUT / "h12_window_sensitivity.csv", index=False)
    ax_ok = lg.groupby("helix").axis_status.value_counts().unstack(fill_value=0)
    num = {"ranges": {k: list(v) for k, v in rng.items()}, "halves": {k: list(v) for k, v in halves.items()},
           "h12_sensitivity_range": list(h12_sv02), "n_chains": int(len(wd)),
           "groups": wd.group.value_counts().to_dict(), "frame_max_abs_diff_A": chk["max_abs_diff_backbone_CA_vs_aligned_A"],
           "sv02_crosscheck": cmp, "axis_status": {h: ax_ok.loc[h].to_dict() for h in ax_ok.index},
           "n_h12_axis": int(wd.axis_H12.sum()),
           "h12_axis_n_ca_by_group": {g: {int(k): int(v) for k, v in x.n_ca_axis.value_counts().sort_index().items()}
                                      for g, x in lg[(lg.helix == "H12") & (lg.axis_status == "ok")].groupby("group")},
           "n_h12_map_points": int(wd[["dist_H12c_H11axis_A", "angle_H12_H3_deg"]].notna().all(axis=1).sum()),
           "descriptor_n": {c: int(wd[c].notna().sum()) for c in wd.columns
                            if c.startswith(("angle_", "dist_", "kink_", "hinge_", "ulen_"))}}
    (OUT / "numbers_geometry.json").write_text(json.dumps(num, indent=1, default=lambda v: v.item()
                                                          if hasattr(v, "item") else str(v)))
    print(json.dumps({k: num[k] for k in ("frame_max_abs_diff_A", "sv02_crosscheck", "n_h12_axis",
                                         "n_h12_map_points", "descriptor_n")}, indent=1, default=str))
    print(ax_ok.to_string())
    print(sm_[["helix", "range", "n_with_run", "start_censored", "end_censored", "start_gap", "end_gap",
               "n_uncensored", "censoring_fraction", "median_uncensored_length_cl7",
               "median_uncensored_length_cl4", "n_axis"]].to_string())


if __name__ == "__main__":
    main()
