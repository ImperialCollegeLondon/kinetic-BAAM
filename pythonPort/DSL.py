"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Imperial College London, United Kingdom
Multiphase Systems Laboratory
Year:     2025
MATLAB:   R2024a
Authors:  Hassan Azzan (HA), Ayca Yilmaz (AY)

Purpose:
Function that takes pressure, mole fraction of component 1, temperature,
and isotherm parameters. Outputs corresponding equilibrium adsorbed
amounts for both components.

Last modified:
- 2025-09-17, HA: Initial creation

Input arguments:
    - P: pressure [Pa]
    - y1: mole fraction of component 1 [-]
    - T: temperature [K]
    - DSL isotherm parameters

Output arguments:
    - q1_star: equilibrium adsorbed amount for component 1 [mol/kg]
    - q2_star: equilibrium adsorbed amount for component 2 [mol/kg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"""
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
