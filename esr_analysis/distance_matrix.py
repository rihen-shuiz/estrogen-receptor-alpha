"""Validate the saved square RMSD contract used by clustering and PCA. No I/O."""
from __future__ import annotations
import numpy as np

DISTANCE_TOLERANCE_A = 1e-8


def validate_distance_matrix(
    rmsd_matrix: np.ndarray, expected_size: int | None = None,
) -> np.ndarray:
    """Validate the historical distance contract before linkage or PCA."""
    try:
        matrix = np.asarray(rmsd_matrix, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("RMSD distance matrix must contain numeric values") from exc
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("RMSD distance matrix must be 2D and square")
    if expected_size is not None and matrix.shape[0] != expected_size:
        raise ValueError("RMSD distance matrix size must match the number of chain keys")
    if not np.isfinite(matrix).all():
        raise ValueError("RMSD distance matrix must contain only finite values")
    if not np.allclose(matrix, matrix.T, rtol=0.0, atol=DISTANCE_TOLERANCE_A):
        max_asymmetry = float(np.max(np.abs(matrix - matrix.T)))
        raise ValueError(
            "RMSD distance matrix must be symmetric within "
            f"{DISTANCE_TOLERANCE_A:g} A (maximum asymmetry {max_asymmetry:g} A)"
        )
    if not np.allclose(np.diag(matrix), 0.0, rtol=0.0, atol=DISTANCE_TOLERANCE_A):
        max_diagonal = float(np.max(np.abs(np.diag(matrix))))
        raise ValueError(
            "RMSD distance matrix diagonal must be zero within "
            f"{DISTANCE_TOLERANCE_A:g} A (maximum diagonal {max_diagonal:g} A)"
        )
    if np.any(matrix < -DISTANCE_TOLERANCE_A):
        minimum = float(matrix.min())
        raise ValueError(
            "RMSD distance matrix must not contain negative distances beyond "
            f"{DISTANCE_TOLERANCE_A:g} A (minimum {minimum:g} A)"
        )
    return matrix

