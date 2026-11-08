"""Helix midpoint estimate for incomplete helices (SV04, Ye's request), built on the SV03 axis rule.

estimate_midpoint(positions, coords, range_start, range_end):
  1. axis = scripts/sv/common/helix_geometry.fit_axis on the given C-alpha (SV03 [DEFAULT]: >= 6 C-alpha, centroids of
     consecutive 4-residue windows, principal line, oriented N -> C, passing through `centre`);
  2. projection of every given C-alpha on the axis, t_i = (x_i - centre) . u, regressed on residue number by ordinary
     least squares, t = a + b * residue (b = rise per residue, A);
  3. axis point of residue r: P(r) = centre + (a + b * r) * u; start / end point = P(range_start) / P(range_end);
     midpoint = P((range_start + range_end) / 2), the axis point of the middle residue of the consensus helix
     (= the middle of the start and end points, P is linear in r).
The chain is then represented per helix by the midpoint and the unit axis vector u.
"""
from __future__ import annotations

import numpy as np

from .helix_geometry import MIN_CA, fit_axis


def estimate_midpoint(positions, coords, range_start: int, range_end: int, min_ca: int = MIN_CA) -> dict:
    positions = np.asarray(positions, int)
    coords = np.asarray(coords, float).reshape(-1, 3)
    ax = fit_axis(positions, coords, min_ca)
    nan3 = np.full(3, np.nan)
    out = {"status": ax["status"], "n_ca": ax["n_ca"], "n_windows": ax["n_windows"], "direction": ax["direction"],
           "centre": ax["centre"], "rise_per_residue_A": np.nan, "intercept_A": np.nan, "fit_resid_sd_A": np.nan,
           "mid_residue": (range_start + range_end) / 2.0, "start_point": nan3, "end_point": nan3, "midpoint": nan3,
           "rms_window_offset": ax["rms_window_offset"]}
    if ax["status"] != "ok":
        return out
    c, u = ax["centre"], ax["direction"]
    t = (coords - c) @ u
    b, a = np.polyfit(positions.astype(float), t, 1)
    resid = t - (a + b * positions)

    def point(r):
        return c + (a + b * r) * u

    out.update(rise_per_residue_A=float(b), intercept_A=float(a), fit_resid_sd_A=float(resid.std(ddof=0)),
               start_point=point(range_start), end_point=point(range_end), midpoint=point(out["mid_residue"]))
    return out
