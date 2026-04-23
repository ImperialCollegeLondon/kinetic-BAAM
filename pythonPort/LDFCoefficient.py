"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Imperial College London, United Kingdom
Multiphase Systems Laboratory
Year:     2025
MATLAB:   R2024a
Authors:  Hassan Azzan (HA), Ayca Yilmaz (AY)

Purpose:
Function that takes pressure, mole fraction of component 1, temperature,
corresponding equilibrium adsorbed amounts, and parameters. Outputs
concentration dependent mass transfer coefficients.

Last modified:
- 2025-09-17, HA: Initial creation

Input arguments:
    - P: pressure [Pa]
    - y1: mole fraction of component 1 [-]
    - T: temperature [K]
    - q1_star: equilibrium adsorbed amount of component 1 at P, y1, T
    - q2_star: equilibrium adsorbed amount of component 2 at P, y1, T
    - parameters: contains adsorbent properties and process parameters

Output arguments:
    - k1: concentration dependent mass transfer coefficient for component 1 [1/s]
    - k2: concentration dependent mass transfer coefficient for component 2 [1/s]

Dependencies:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"""
from __future__ import annotations
import numpy as np


def LDFCoefficient(P, y1, T, q1_star, q2_star, parameters):
    Rg = 8.3145

    rp = parameters['rp'] if isinstance(parameters, dict) else getattr(parameters, 'rp')
    Dp = parameters['Dm'] / parameters['tau'] if isinstance(parameters, dict) else getattr(parameters, 'Dm') / getattr(parameters, 'tau')

    if parameters.get('LDFtest', False):
        Dp = Dp * parameters.get('LDFFactor', 1.0)

    if parameters.get('equilibrium', False):
        Dp = Dp * 100.0

    k01 = 15 * parameters['epsilon_p'] * Dp / (rp ** 2)
    k02 = 15 * parameters['epsilon_p'] * Dp / (rp ** 2)

    eps = 1e-20

    # compute ratios carefully to avoid divide-by-zero
    try:
        if y1 > 0 and y1 < 1:
            ratio1 = ((P * y1 / (Rg * T)) / (q1_star + eps)) / parameters['rho_s']
            ratio2 = ((P * (1 - y1) / (Rg * T)) / (q2_star + eps)) / parameters['rho_s']
        elif y1 <= 0:
            ratio2 = ((P * 1.0 / (Rg * T)) / (q2_star + eps)) / parameters['rho_s']
            ratio1 = 1.0 / (
                parameters['qsb_1'] * parameters['bo_1'] * np.exp(-parameters['delUb_1'] / (Rg * T))
                + parameters['qsd_1'] * parameters['do_1'] * np.exp(-parameters['delUd_1'] / (Rg * T))
            ) / parameters['rho_s']
        elif y1 >= 1:
            ratio1 = ((P * 1.0 / (Rg * T)) / (q1_star + eps)) / parameters['rho_s']
            ratio2 = 1.0 / (
                parameters['qsb_2'] * parameters['bo_2'] * np.exp(-parameters['delUb_2'] / (Rg * T))
                + parameters['qsd_2'] * parameters['do_2'] * np.exp(-parameters['delUd_2'] / (Rg * T))
            ) / parameters['rho_s']
        else:
            ratio1 = 1.0
            ratio2 = 1.0
    except Exception:
        # fallback to safe defaults
        ratio1 = 1.0
        ratio2 = 1.0

    # handle NaNs
    if np.isnan(ratio1):
        ratio1 = 1.0 / (
            parameters['qsb_1'] * parameters['bo_1'] * np.exp(-parameters['delUb_1'] / (Rg * T))
            + parameters['qsd_1'] * parameters['do_1'] * np.exp(-parameters['delUd_1'] / (Rg * T))
        ) / parameters['rho_s']
    if np.isnan(ratio2):
        ratio2 = 1.0 / (
            parameters['qsb_2'] * parameters['bo_2'] * np.exp(-parameters['delUb_2'] / (Rg * T))
            + parameters['qsd_2'] * parameters['do_2'] * np.exp(-parameters['delUd_2'] / (Rg * T))
        ) / parameters['rho_s']

    k1 = k01 * ratio1
    k2 = k02 * ratio2

    if parameters.get('processType') == 'Resin' or parameters.get('amine', False):
        k01_ref = parameters['LDF'] / np.exp(-38.87e3 / (8.314 * 303))
        k1 = k01_ref * np.exp(-38.87e3 / (8.314 * T))

    return k1, k2
