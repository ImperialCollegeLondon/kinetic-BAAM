"""DSL.py

Imperial College London, Multiphase Systems Laboratory
Year: 2025
Author: Hassan Azzan (HA)

Purpose:
    Port of MATLAB `DSL.m` — implements the competitive DSL adsorption
    isotherm used by the k-BAAM model.

Function:
    DSL(P, y1, T, qsb_1, qsd_1, qsb_2, qsd_2, bo_1, do_1, bo_2, do_2,
        delUb_1, delUd_1, delUb_2, delUd_2) -> (q1_star, q2_star)

Notes:
    - P expected in Pa, T in K. Function is vectorized via numpy where
      possible.
"""
from __future__ import annotations
import numpy as np


def DSL(P, y1, T, qsb_1, qsd_1, qsb_2, qsd_2, bo_1, do_1, bo_2, do_2,
        delUb_1, delUd_1, delUb_2, delUd_2):
    """Compute competitive DSL equilibrium adsorbed amounts.

    Inputs follow the MATLAB function signature. P should be in Pa and T in K.
    The function is vectorized via numpy where appropriate.
    """
    R = 8.314

    P = np.asarray(P, dtype=float)
    y1 = np.asarray(y1, dtype=float)
    T = np.asarray(T, dtype=float)

    b1 = bo_1 * np.exp(-delUb_1 / (R * T))
    d1 = do_1 * np.exp(-delUd_1 / (R * T))
    b2 = bo_2 * np.exp(-delUb_2 / (R * T))
    d2 = do_2 * np.exp(-delUd_2 / (R * T))

    # avoid division by zero in denominators by using numpy operations
    c1 = P * y1 / (R * T)
    c2 = P * (1 - y1) / (R * T)

    term1_num = qsb_1 * b1 * c1
    term1_den = 1 + b1 * c1 + b2 * c2
    term1 = term1_num / term1_den

    term2_num = qsd_1 * d1 * c1
    term2_den = 1 + d1 * c1 + d2 * c2
    term2 = term2_num / term2_den

    q1_star = term1 + term2

    term1_num_2 = qsb_2 * b2 * c2
    term1_den_2 = 1 + b2 * c2 + b1 * c1
    term1_2 = term1_num_2 / term1_den_2

    term2_num_2 = qsd_2 * d2 * c2
    term2_den_2 = 1 + d1 * c1 + d2 * c2
    term2_2 = term2_num_2 / term2_den_2

    q2_star = term1_2 + term2_2

    return q1_star, q2_star
