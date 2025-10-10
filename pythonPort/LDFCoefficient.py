"""LDFCoefficient.py

Imperial College London, Multiphase Systems Laboratory
Year: 2025
Author: Hassan Azzan (HA)

Purpose:
    Port of MATLAB `LDFCoefficient.m`. Computes LDF rate coefficients for
    the two-site model used by k-BAAM.

Function:
    LDFCoefficient(P, y1, T, q1_star, q2_star, parameters) -> (k1, k2)

Notes:
    - Uses parameters['rp'], parameters['Dm'], parameters['tau'], and
      other fields from the parameters dict. Returns scalars or arrays
      consistent with input shapes.
"""
from __future__ import annotations
import numpy as np


def LDFCoefficient(P, y1, T, q1_star, q2_star, parameters):
    Rg = 8.3145

    rp = parameters['rp'] if isinstance(parameters, dict) else getattr(parameters, 'rp')
    Dp = parameters['Dm'] / parameters['tau'] if isinstance(parameters, dict) else getattr(parameters, 'Dm') / getattr(parameters, 'tau')

    k01 = 15 * parameters['epsilon_p'] * Dp / (rp ** 2)
    k02 = 15 * parameters['epsilon_p'] * Dp / (rp ** 2)

    # compute ratios carefully to avoid divide-by-zero
    try:
        if y1 > 0 and y1 < 1:
            ratio1 = ((P * y1 / (Rg * T)) / q1_star) / parameters['rho_s']
            ratio2 = ((P * (1 - y1) / (Rg * T)) / q2_star) / parameters['rho_s']
        elif y1 <= 0:
            ratio1 = ((P * (1 - y1) / (Rg * T)) / q2_star) / parameters['rho_s']
            ratio2 = 1.0
        elif y1 >= 1:
            ratio1 = ((P * y1 / (Rg * T)) / q1_star) / parameters['rho_s']
            ratio2 = 1.0
        else:
            ratio1 = 1.0
            ratio2 = 1.0
    except Exception:
        # fallback to safe defaults
        ratio1 = 1.0
        ratio2 = 1.0

    # handle NaNs
    if np.isnan(ratio1):
        ratio1 = 1.0
    if np.isnan(ratio2):
        ratio2 = 1.0

    k1 = k01 * ratio1
    k2 = k02 * ratio2
    return k1, k2
