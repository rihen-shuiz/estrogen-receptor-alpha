"""SV04 steps 1-3 (RUN_SV1028): helix midpoints for incomplete helices (Ye), truncation validation, eligibility.

Inputs (tracker cards SV00 / SV02 / SV03 `out`, rule 1.10; SV03 is DONE, so its code and frame are reused):
- results/sv/SV02/helix_boundaries.csv: consensus helices (global_start-global_end); H12 = 538-547, the primary H12
  window of run file section 2 (the SV02 data-derived 538-544 is a proposal that SV01c tests), as in SV03.
- C-alpha in the core v3 frame: scripts/sv/common/core_frame.load_ca_core() (SV00 gate fit, frozen by SV03 in
  data/sv/SV03/core_v3_transforms.csv; SV03 checked deposited backbone x frozen transforms = these coordinates).
- data/sv/SV02/dssp_residues.csv.gz: pydssp state per residue (read with SV03_geometry.chain_residues).
- results/sv/SV00/chain_table.csv (s08_cluster = S08 cut 1.5 A), results/sv/SV01a/chain_outcomes.csv (h12_status).

Rules fixed before any output was looked at:
- Step 1: C-alpha of a helix in a chain = positions of the consensus range that are observed and pydssp 'H' (the
  SV03 axis rule "observed helical C-alpha"; step 1 says "fit the axis as in SV03"), so every SV04 axis equals the
  SV03 axis of the same chain x helix. A residue that is observed but not helical counts as missing. A helix needs
  >= 6 such C-alpha (and >= 2 consecutive 4-residue windows, SV03).
- Step 2: scripts/sv/common/helix_midpoint.estimate_midpoint (projection on the axis regressed on residue number;
  axis points of the consensus start and end; midpoint = axis point of the middle residue (start + end) / 2).
- Step 3 validation: a helix is complete in a chain when every position of its consensus range is observed and
  helical (k = 0). On every complete helix, delete k = 1..6 residues from the N end, and separately from the C end,
  re-estimate with the same procedure (needs >= 6 C-alpha left) and measure the midpoint error (A) and the axis angle
  error (deg) against the estimate from the full helix. Error curve per helix x end (N, C, either = both pooled) x k.
- [DEFAULT] eligibility: the observed length of a chain's helix is n = number of its C-alpha used (k = L - n
  residues missing, L = consensus length). The helix is used for the chain only if n >= 6 and the median midpoint
  error at that observed length (pooled over N- and C-end deletions of k residues) is <= 1.0 A; k = 0 (complete) has
  error 0 by definition. Observed lengths with k > 6 were not validated and are not used.
Outputs: results/sv/SV04/validation_error.csv (error curve), midpoints.csv (chain x helix: midpoint, unit axis,
  rise, eligibility), numbers_midpoints.json; data/sv/SV04/validation_truncations.csv.gz (every truncation).
Run: .venv/bin/python -B -I scripts/sv/SV04_midpoints.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from common import core_frame as cf  # noqa: E402
from common.helix_geometry import MIN_CA, axis_angle  # noqa: E402
from common.helix_midpoint import estimate_midpoint  # noqa: E402
from SV03_geometry import chain_residues  # noqa: E402  (SV03 residue / run table, reused)

OUT = ROOT / "results/sv/SV04"
DATA = ROOT / "data/sv/SV04"
H12_PRIMARY = (538, 547)
K_MAX = 6
MAX_MEDIAN_ERR = 1.0   # [DEFAULT] A
HELIX_TAG = {"H1": "H1", "H3": "H3", "H4/5": "H45", "H7": "H7", "H8": "H8", "H9": "H9", "H10/11": "H1011",
             "H12": "H12"}


def helix_ranges():
    hb = pd.read_csv(ROOT / "results/sv/SV02/helix_boundaries.csv")
    rng = {r.helix_name: (int(r.global_start), int(r.global_end)) for r in hb.itertuples()}
    rng["H12"] = H12_PRIMARY
    return rng


def q(x, p):
    return float(np.quantile(x, p)) if len(x) else np.nan


def curve_rows(tr, rng):
    rows = []
    for name, (s, e) in rng.items():
        L = e - s + 1
        t = tr[tr.helix == name]
        n_complete = int(t.chain_key.nunique())
        for end in ("N", "C", "either"):
            for k in range(0, K_MAX + 1):
                x = t[(t.k_removed == k) & ((t.end == end) | (end == "either"))] if k else t[t.k_removed == 0]
                ok = x[x.status == "ok"]
                me, ae = ok.midpoint_error_A.to_numpy(), ok.axis_angle_error_deg.to_numpy()
                rows.append({"helix": name, "range_start": s, "range_end": e, "consensus_length": L, "end": end,
                             "k_removed": k, "observed_length": L - k, "n_complete_chains": n_complete,
                             "n_truncations": int(len(x)), "n_estimated": int(len(ok)),
                             "midpoint_err_median_A": q(me, 0.5), "midpoint_err_q25_A": q(me, 0.25),
                             "midpoint_err_q75_A": q(me, 0.75), "midpoint_err_q90_A": q(me, 0.9),
                             "midpoint_err_max_A": float(me.max()) if len(me) else np.nan,
                             "axis_err_median_deg": q(ae, 0.5), "axis_err_q25_deg": q(ae, 0.25),
                             "axis_err_q75_deg": q(ae, 0.75), "axis_err_q90_deg": q(ae, 0.9),
                             "axis_err_max_deg": float(ae.max()) if len(ae) else np.nan})
    c = pd.DataFrame(rows)
    # [DEFAULT] eligibility per helix x observed length (pooled N/C rows only)
    piv = c.end == "either"
    c["eligible_observed_length"] = np.where(
        piv, (c.observed_length >= MIN_CA) & (c.midpoint_err_median_A <= MAX_MEDIAN_ERR) & (c.n_estimated > 0), pd.NA)
    c.loc[piv & (c.k_removed == 0), "eligible_observed_length"] = c.loc[piv & (c.k_removed == 0), "n_complete_chains"] > 0
    return c


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    rng = helix_ranges()
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", keep_default_na=False, na_values=[""])
    oc = pd.read_csv(ROOT / "results/sv/SV01a/chain_outcomes.csv").set_index("chain_key")
    keys, pos, ca_core = cf.load_ca_core()
    assert list(keys) == list(ct.chain_key)
    pos_index = {int(p): i for i, p in enumerate(pos)}
    res = pd.read_csv(ROOT / "data/sv/SV02/dssp_residues.csv.gz")
    obs, _ = chain_residues(res)
    group = ct.s08_cluster.map({7: "cl7", 4: "cl4"}).fillna("minor")

    est_rows, trunc_rows = [], []
    for i, r in ct.iterrows():
        key, ca_i = r.chain_key, ca_core[i]
        o = obs[key]
        for name, (s, e) in rng.items():
            L = e - s + 1
            hp = [p for p in range(s, e + 1) if p in o and o[p][1] == "H" and not np.isnan(ca_i[pos_index[p], 0])]
            n_obs = sum(1 for p in range(s, e + 1) if not np.isnan(ca_i[pos_index[p], 0]))
            xyz = np.array([ca_i[pos_index[p]] for p in hp]).reshape(-1, 3)
            m = estimate_midpoint(hp, xyz, s, e)
            k = L - len(hp)
            est_rows.append({"chain_key": key, "entry_id": r.entry_id, "s08_cluster": int(r.s08_cluster),
                             "group": group[i], "h12_status": oc.loc[key, "h12_status"], "helix": name,
                             "helix_tag": HELIX_TAG[name], "range_start": s, "range_end": e, "consensus_length": L,
                             "mid_residue": m["mid_residue"], "n_observed_ca": n_obs, "n_fit_ca": len(hp),
                             "k_missing": k, "fit_first": hp[0] if hp else np.nan,
                             "fit_last": hp[-1] if hp else np.nan, "n_windows": m["n_windows"],
                             "axis_status": m["status"], "rise_per_residue_A": m["rise_per_residue_A"],
                             "fit_resid_sd_A": m["fit_resid_sd_A"],
                             **{f"mid_{a}": v for a, v in zip("xyz", m["midpoint"])},
                             **{f"u_{a}": v for a, v in zip("xyz", m["direction"])},
                             **{f"start_{a}": v for a, v in zip("xyz", m["start_point"])},
                             **{f"end_{a}": v for a, v in zip("xyz", m["end_point"])}})
            if k != 0:
                continue
            # step 3: truncations of a complete helix (all positions observed and helical)
            trunc_rows.append({"chain_key": key, "entry_id": r.entry_id, "group": group[i], "helix": name,
                               "end": "none", "k_removed": 0, "n_fit_ca": len(hp), "status": "ok",
                               "midpoint_error_A": 0.0, "axis_angle_error_deg": 0.0, "rise_per_residue_A":
                               m["rise_per_residue_A"]})
            for end in ("N", "C"):
                for kk in range(1, K_MAX + 1):
                    sel = slice(kk, None) if end == "N" else slice(None, len(hp) - kk)
                    hp_t, xyz_t = hp[sel], xyz[sel]
                    mt = estimate_midpoint(hp_t, xyz_t, s, e)
                    ok = mt["status"] == "ok"
                    trunc_rows.append({"chain_key": key, "entry_id": r.entry_id, "group": group[i], "helix": name,
                                       "end": end, "k_removed": kk, "n_fit_ca": len(hp_t), "status": mt["status"],
                                       "midpoint_error_A": float(np.linalg.norm(mt["midpoint"] - m["midpoint"]))
                                       if ok else np.nan,
                                       "axis_angle_error_deg": axis_angle(mt["direction"], m["direction"])
                                       if ok else np.nan, "rise_per_residue_A": mt["rise_per_residue_A"]})
    est = pd.DataFrame(est_rows)
    tr = pd.DataFrame(trunc_rows)
    tr.to_csv(DATA / "validation_truncations.csv.gz", index=False)
    curve = curve_rows(tr, rng)
    curve.to_csv(OUT / "validation_error.csv", index=False)

    # [DEFAULT] eligibility of each chain x helix by its observed length
    ev = curve[curve.end == "either"].set_index(["helix", "k_removed"])
    med = {(h, k): v for (h, k), v in ev.midpoint_err_median_A.items()}
    elig_k = {(h, k): bool(v) for (h, k), v in ev.eligible_observed_length.items()}

    def rule(rw):
        if rw.axis_status != "ok":
            return False, np.nan, "no axis (" + rw.axis_status.replace("no axis: ", "") + ")"
        if rw.k_missing > K_MAX:
            return False, np.nan, f"not validated (k = {rw.k_missing} > {K_MAX})"
        mm = med[(rw.helix, rw.k_missing)]
        if not elig_k[(rw.helix, rw.k_missing)]:
            return False, mm, f"median midpoint error {mm:.2f} A > {MAX_MEDIAN_ERR} A at k = {rw.k_missing}"
        return True, mm, "complete" if rw.k_missing == 0 else f"eligible (k = {rw.k_missing})"

    rr = est.apply(rule, axis=1, result_type="expand")
    est["validated_median_err_A"], est["eligible"], est["eligibility_reason"] = rr[1], rr[0], rr[2]
    est = est[[c for c in est.columns if c not in ("eligible", "validated_median_err_A", "eligibility_reason")] +
              ["validated_median_err_A", "eligible", "eligibility_reason"]]
    est.to_csv(OUT / "midpoints.csv", index=False)

    # cross-check: SV04 axes = SV03 axes (same C-alpha, same rule)
    h3 = pd.read_csv(ROOT / "results/sv/SV03/helix_per_chain.csv")
    h3 = h3[h3.helix.isin(rng)].set_index(["chain_key", "helix"])
    m4 = est.set_index(["chain_key", "helix"]).loc[h3.index]
    same_n = int((m4.n_fit_ca.to_numpy() == h3.n_ca_axis.to_numpy()).sum())
    okb = (h3.axis_status == "ok").to_numpy() & (m4.axis_status == "ok").to_numpy()
    u3 = h3[["axis_ux", "axis_uy", "axis_uz"]].to_numpy()[okb]
    u4 = m4[["u_x", "u_y", "u_z"]].to_numpy()[okb]
    max_ang = float(np.degrees(np.arccos(np.clip((u3 * u4).sum(1), -1, 1))).max())
    # distance of the SV04 midpoint from the SV03 axis line through the SV03 centroid (should be ~0)
    c3 = h3[["centroid_x", "centroid_y", "centroid_z"]].to_numpy()[okb]
    m4p = m4[["mid_x", "mid_y", "mid_z"]].to_numpy()[okb]
    d = m4p - c3
    perp = np.linalg.norm(d - (d * u3).sum(1)[:, None] * u3, axis=1)
    along = (d * u3).sum(1)

    num = {"ranges": {k: list(v) for k, v in rng.items()}, "n_chains": int(ct.shape[0]), "k_max": K_MAX,
           "max_median_err_A": MAX_MEDIAN_ERR,
           "sv03_crosscheck": {"n_rows": int(len(h3)), "n_fit_ca_equal": same_n, "n_axis_both": int(okb.sum()),
                               "n_axis_status_equal": int(((h3.axis_status == "ok").to_numpy() ==
                                                           (m4.axis_status == "ok").to_numpy()).sum()),
                               "max_axis_angle_diff_deg": max_ang,
                               "max_midpoint_offset_from_sv03_axis_A": float(perp.max()),
                               "midpoint_minus_sv03_centroid_along_axis_A": {
                                   h: {"median": float(np.median(along[(m4.reset_index().helix.to_numpy()[okb] == h)])),
                                       "max_abs": float(np.abs(along[(m4.reset_index().helix.to_numpy()[okb] == h)]).max())}
                                   for h in rng}},
           "complete_helices": {h: {"n": int(((est.helix == h) & (est.k_missing == 0)).sum()),
                                    **{g: int(((est.helix == h) & (est.k_missing == 0) & (est.group == g)).sum())
                                       for g in ("cl7", "cl4", "minor")}} for h in rng},
           "k_distribution": {h: {int(k): int(v) for k, v in est[est.helix == h].k_missing.value_counts().sort_index().items()}
                              for h in rng},
           "eligible_k": {h: [int(k) for (hh, k), v in elig_k.items() if hh == h and v] for h in rng},
           "eligible_by_helix": {h: {"eligible": int(est[(est.helix == h)].eligible.sum()),
                                     "axis_ok": int((est[(est.helix == h)].axis_status == "ok").sum()),
                                     "reasons": est[est.helix == h].eligibility_reason.str.replace(
                                         r" \d+\.\d+ A", " x A", regex=True).value_counts().to_dict()}
                                 for h in rng},
           "eligible_helices_per_chain": {int(k): int(v) for k, v in
                                          est.groupby("chain_key").eligible.sum().value_counts().sort_index().items()},
           "rise_per_residue_A_complete": {h: float(est[(est.helix == h) & (est.k_missing == 0)].rise_per_residue_A.median())
                                           for h in rng}}
    hs = est[est.helix == "H12"]
    num["h12_by_status"] = {st: {"n": int(len(x)), "axis_ok": int((x.axis_status == "ok").sum()),
                                 "eligible": int(x.eligible.sum()),
                                 "eligible_by_group": {g: int(x[x.group == g].eligible.sum()) for g in ("cl7", "cl4", "minor")},
                                 "n_by_group": {g: int((x.group == g).sum()) for g in ("cl7", "cl4", "minor")},
                                 "k_distribution": {int(k): int(v) for k, v in x.k_missing.value_counts().sort_index().items()}}
                            for st, x in hs.groupby("h12_status")}
    (OUT / "numbers_midpoints.json").write_text(json.dumps(num, indent=1, default=lambda v: v.item()
                                                           if hasattr(v, "item") else str(v)))
    pd.set_option("display.width", 220)
    print(json.dumps({k: num[k] for k in ("sv03_crosscheck", "complete_helices", "eligible_k", "h12_by_status")},
                     indent=1))
    print(curve[curve.end == "either"][["helix", "k_removed", "observed_length", "n_estimated",
                                         "midpoint_err_median_A", "midpoint_err_q90_A", "axis_err_median_deg",
                                         "eligible_observed_length"]].to_string())


if __name__ == "__main__":
    main()
