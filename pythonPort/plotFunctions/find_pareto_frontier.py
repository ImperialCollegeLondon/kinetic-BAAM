from __future__ import annotations
import numpy as np
from typing import Tuple


def find_pareto_frontier(input_arr: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Pareto frontier for minimization problems.

    Mirrors the MATLAB function by:
    - considering only the first two columns of input
    - removing duplicate rows
    - returning a flag mask over the original input and the Pareto points

    Returns
    -------
    flag : np.ndarray of bool, shape (n,)
        True where the corresponding row in the original input is Pareto-efficient
    pareto_points : np.ndarray, shape (k, 2)
        Unique 2D Pareto points (minimization)
    """
    if input_arr.ndim != 2 or input_arr.shape[1] < 2:
        raise ValueError("input_arr must be (n, >=2)")

    data = input_arr[:, :2]
    # unique rows
    data_unique = np.unique(data, axis=0)

    out = []
    m = data_unique.shape[0]
    for i in range(m):
        c = np.broadcast_to(data_unique[i], data_unique.shape)
        t = data_unique.copy()
        t[i, :] = np.inf
        smaller_idx = c >= t
        idx = np.sum(smaller_idx, axis=1) == data_unique.shape[1]
        if not np.any(idx):
            out.append(data_unique[i])

    pareto_points = np.array(out).reshape(-1, 2) if out else np.empty((0, 2))
    # map back to original rows
    def _row_in_pareto(row):
        return np.any(np.all(pareto_points == row[:2], axis=1)) if pareto_points.size else False

    flag = np.array([_row_in_pareto(row) for row in data])
    return flag, pareto_points
