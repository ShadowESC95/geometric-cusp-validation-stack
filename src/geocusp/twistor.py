from __future__ import annotations

import numpy as np


def standard_external_Z4(dtype=float) -> np.ndarray:
    return np.eye(4, dtype=dtype)


def line_from_xyzw(x, y, z0, w, dtype=None) -> np.ndarray:
    if dtype is None:
        dtype = np.result_type(x, y, z0, w, float)
    return np.array([[1, x, 0, -y], [0, z0, 1, w]], dtype=dtype)


def four_bracket(A: np.ndarray, B: np.ndarray, i: int, j: int, Z: np.ndarray | None = None):
    """4-bracket <A B i j>; supports real or complex arithmetic."""
    if Z is None:
        Z = standard_external_Z4(dtype=np.result_type(A, B))
    val = np.linalg.det(np.column_stack([A, B, Z[:, i] if Z.shape[1] == 4 else Z[i], Z[:, j] if Z.shape[1] == 4 else Z[j]]))
    return val


def all_four_brackets(line: np.ndarray, Z: np.ndarray | None = None) -> dict[str, complex]:
    if Z is None:
        Z = standard_external_Z4(dtype=np.result_type(line))
    A, B = line
    pairs = [(1,2),(2,3),(3,4),(1,4),(1,3),(2,4)]
    return {f"{i}{j}": four_bracket(A, B, i-1, j-1, Z) for i,j in pairs}


def cross_ratio(line: np.ndarray, Z: np.ndarray | None = None):
    b = all_four_brackets(line, Z)
    return b["12"] * b["34"] / (b["23"] * b["14"])


def mutual_bracket(line1: np.ndarray, line2: np.ndarray):
    A, B = line1; C, D = line2
    return np.linalg.det(np.column_stack([A, B, C, D]))


def positive_chart_bounds(line: np.ndarray) -> dict[str, complex]:
    return all_four_brackets(line)


def strict_cusp_line(i: int = 2, alpha: complex = 0.7) -> np.ndarray:
    Z = standard_external_Z4(dtype=complex)
    return np.vstack([Z[i-1], Z[i-2] + alpha * Z[i]])
