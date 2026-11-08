"""Helix axis geometry shared by SV03 and SV04 (RUN_SV1028).

SV03 step 2 [DEFAULT] axis rule:
  - a helix gets an axis only with >= 6 observed helical C-alpha (MIN_CA);
  - points = centroids of consecutive 4-residue windows (i, i+1, i+2, i+3 all present in the residue set);
  - axis = principal line through those points (first right-singular vector of the centred points), oriented N -> C
    (positive projection of the last window centroid minus the first); at least 2 window centroids are needed to
    define a line.
  - the axis passes through `centre` = mean of the window centroids (reported as the helix centroid, a point on the
    axis).
Angles between oriented axes are in 0-180 degrees; point-to-axis distances are perpendicular distances to the
infinite axis line. Angles and distances do not depend on the frame.
"""
from __future__ import annotations

import numpy as np

MIN_CA = 6
WINDOW = 4


def window_centroids(positions, coords, window: int = WINDOW):
    """positions: sorted int array; coords: [n, 3]. Centroids of every run of `window` consecutive positions."""
    positions = np.asarray(positions, int)
    coords = np.asarray(coords, float)
    pts, mids = [], []
    for i in range(len(positions) - window + 1):
        if positions[i + window - 1] - positions[i] == window - 1:
            pts.append(coords[i:i + window].mean(axis=0))
            mids.append(positions[i:i + window].mean())
    return np.array(pts).reshape(-1, 3), np.array(mids, float)


def fit_axis(positions, coords, min_ca: int = MIN_CA, window: int = WINDOW) -> dict:
    """Axis of one helix from its observed helical C-alpha (positions sorted, coords [n, 3] without NaN).

    Returns a dict with status ('ok' or the reason), n_ca, n_windows, centre [3], direction [3] (unit, N -> C),
    rms_window_offset (RMS distance of the window centroids from the axis line, A) and span_A (projection range of
    the window centroids along the axis).
    """
    positions = np.asarray(positions, int)
    coords = np.asarray(coords, float)
    out = {"n_ca": int(len(positions)), "n_windows": 0, "centre": np.full(3, np.nan), "direction": np.full(3, np.nan),
           "rms_window_offset": np.nan, "span_A": np.nan}
    if len(positions) < min_ca:
        out["status"] = f"no axis: < {min_ca} observed helical CA"
        return out
    pts, _ = window_centroids(positions, coords, window)
    out["n_windows"] = int(len(pts))
    if len(pts) < 2:
        out["status"] = "no axis: < 2 consecutive 4-residue windows"
        return out
    centre = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - centre)
    u = vt[0] / np.linalg.norm(vt[0])
    if np.dot(pts[-1] - pts[0], u) < 0:
        u = -u
    proj = (pts - centre) @ u
    resid = (pts - centre) - np.outer(proj, u)
    out.update(status="ok", centre=centre, direction=u,
               rms_window_offset=float(np.sqrt((resid ** 2).sum(axis=1).mean())), span_A=float(proj.max() - proj.min()))
    return out


def axis_angle(u, v) -> float:
    """Angle between two oriented unit axes, degrees (0-180); NaN if either is missing."""
    u, v = np.asarray(u, float), np.asarray(v, float)
    if np.isnan(u).any() or np.isnan(v).any():
        return np.nan
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v), -1.0, 1.0))))


def point_axis_distance(p, centre, direction) -> float:
    """Perpendicular distance (A) from point p to the infinite line through centre along direction."""
    p, c, u = (np.asarray(x, float) for x in (p, centre, direction))
    if np.isnan(p).any() or np.isnan(c).any() or np.isnan(u).any():
        return np.nan
    d = p - c
    return float(np.linalg.norm(d - np.dot(d, u) * u))
