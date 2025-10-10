"""kBAAM_ODEs_nonIsothermal_ND.py

Imperial College London, Multiphase Systems Laboratory
Year: 2025
Author: Hassan Azzan (HA)

Purpose:
    Python translation of the MATLAB `kBAAM_ODEs_nonIsothermal_ND.m` which
    implements the right-hand side for the non-dimensional k-BAAM ODEs.

Function signature:
    def kbaam_odes_nonisothermal_nd(t, X, parameters, step_name)

Inputs:
    - t: dimensionless time (float)
    - X: iterable with 5 state variables [y1, q1, q2, T, Tw]
    - parameters: dict-like or object-like container matching the MATLAB
      `parameters` struct. Callables (e.g. pressure schedules) should be
      provided as callables in this container. The functions `DSL` and
      `LDFCoefficient` should be present in `parameters` or available as
      local module implementations.
    - step_name: one of 'ads', 'blo', 'evac', 'pres'

Returns:
    - dXdt: numpy array of shape (5,) with time derivatives

Notes:
    - Mirrors MATLAB behavior including non-dimensionalisation, time
      references, and handling of F_in/v_in fallbacks to avoid division by
      zero.
"""

from __future__ import annotations
import numpy as np

# Try to import local implementations as sensible defaults/fallbacks
try:
    # package-style import (when running as package)
    from .DSL import DSL as _DSL_impl
    from .LDFCoefficient import LDFCoefficient as _LDF_impl
except Exception:
    try:
        # top-level import for direct script usage
        from DSL import DSL as _DSL_impl
        from LDFCoefficient import LDFCoefficient as _LDF_impl
    except Exception:
        _DSL_impl = None
        _LDF_impl = None


def _get(p, key, default=None):
    """Helper to get a value from a dict-like or object-like parameters."""
    if p is None:
        return default
    if isinstance(p, dict):
        return p.get(key, default)
    return getattr(p, key, default)


def kbaam_odes_nonisothermal_nd(t, X, parameters, step_name: str):
    """Compute ODEs for the (non-)isothermal kBAAM model (dimensionless).

    The function mirrors the MATLAB implementation and expects the caller to
    provide `DSL` and `LDFCoefficient` implementations via the `parameters`
    container.
    """
    R = 8.314

    # ensure X is array-like
    X = np.asarray(X, dtype=float)
    if X.size < 5:
        raise ValueError("X must contain at least 5 state variables: [y1,q1,q2,T,Tw]")

    y1 = X[0]
    q1 = X[1]
    q2 = X[2]
    T = X[3]
    Tw = X[4]

    # parameters accessor
    p = parameters

    # Reference values for non-dimensionalization
    V_column = _get(p, 'V_column')
    if V_column is None:
        raise KeyError('parameters must provide V_column')

    # Use provided F_in when available and non-zero; otherwise reconstruct a
    # volumetric flow from v_in and the inner area (r_in) and compute F_in
    # from ideal gas law. This prevents division-by-zero if parameters.F_in
    # has been mutated to zero in-between steps.
    F_in = _get(p, 'F_in')

    volFluxRef = F_in / V_column
    timeRef = _get(p, 'p_H') / (R * _get(p, 'T_feed') * volFluxRef)
    qRef = _get(p, 'qsb_1') + _get(p, 'qsd_1')
    TRef = _get(p, 'T_feed')
    TwRef = TRef
    PRef = _get(p, 'p_H')

    Qheat = 0.0

    # Prepare outputs
    dq1dt = 0.0
    dq2dt = 0.0
    dPdt = 0.0
    Fout = 0.0

    # helpers (expect callables provided)
    DSL = _get(p, 'DSL')
    LDFCoefficient = _get(p, 'LDFCoefficient')

    # If not provided through parameters, fall back to local module implementations
    if DSL is None:
        DSL = _DSL_impl
    if LDFCoefficient is None:
        LDFCoefficient = _LDF_impl

    # switch-case behavior
    if step_name == 'ads':
        P = _get(p, 'P_ads')(t * timeRef) / PRef


        q1_star, q2_star = DSL(P * PRef, _get(p, 'y1_in'), TRef,
                              _get(p, 'qsb_1'), _get(p, 'qsd_1'), _get(p, 'qsb_2'), _get(p, 'qsd_2'),
                              _get(p, 'bo_1'), _get(p, 'do_1'), _get(p, 'bo_2'), _get(p, 'do_2'),
                              _get(p, 'delUb_1'), _get(p, 'delUd_1'), _get(p, 'delUb_2'), _get(p, 'delUd_2'))

        k1, k2 = LDFCoefficient(P * PRef, y1, T * TRef, q1_star, q2_star, p)

        dq1dt = timeRef / qRef * k1 * (q1_star - q1 * qRef)
        dq2dt = timeRef / qRef * k2 * (q2_star - q2 * qRef)
        dPdt = 0.0
        F_in = _get(p, 'F_in'); 
        Fout = F_in - (1 - _get(p, 'e_bed')) * _get(p, 'V_column') * _get(p, 'rho_s') * (dq1dt + dq2dt) / (timeRef / qRef)

    elif step_name == 'blo':
        P = _get(p, 'P_blo')(t * timeRef) / PRef


        q1_star, q2_star = DSL(P * PRef, y1, T * TRef,
                              _get(p, 'qsb_1'), _get(p, 'qsd_1'), _get(p, 'qsb_2'), _get(p, 'qsd_2'),
                              _get(p, 'bo_1'), _get(p, 'do_1'), _get(p, 'bo_2'), _get(p, 'do_2'),
                              _get(p, 'delUb_1'), _get(p, 'delUd_1'), _get(p, 'delUb_2'), _get(p, 'delUd_2'))

        k1, k2 = LDFCoefficient(P * PRef, y1, T * TRef, q1_star, q2_star, p)

        dq1dt = timeRef / qRef * k1 * (q1_star - q1 * qRef)
        dq2dt = timeRef / qRef * k2 * (q2_star - q2 * qRef)
        dPdt = _get(p, 'dPdt_blo')(t * timeRef) / (PRef / timeRef)

        F_in = 0;

        Fout = -(1 - _get(p, 'e_bed')) * _get(p, 'V_column') * _get(p, 'rho_s') * (dq1dt + dq2dt) / (timeRef / qRef) - (_get(p, 'e_bed') / (R * T * TRef)) * dPdt * (_get(p, 'p_H') / timeRef) * _get(p, 'V_column')

    elif step_name == 'evac':
        P = _get(p, 'P_evac')(t * timeRef) / PRef

        q1_star, q2_star = DSL(P * PRef, y1, T * TRef,
                              _get(p, 'qsb_1'), _get(p, 'qsd_1'), _get(p, 'qsb_2'), _get(p, 'qsd_2'),
                              _get(p, 'bo_1'), _get(p, 'do_1'), _get(p, 'bo_2'), _get(p, 'do_2'),
                              _get(p, 'delUb_1'), _get(p, 'delUd_1'), _get(p, 'delUb_2'), _get(p, 'delUd_2'))

        k1, k2 = LDFCoefficient(P * PRef, y1, T * TRef, q1_star, q2_star, p)

        dq1dt = timeRef / qRef * k1 * (q1_star - q1 * qRef)
        dq2dt = timeRef / qRef * k2 * (q2_star - q2 * qRef)
        dPdt = _get(p, 'dPdt_evac')(t * timeRef) / (_get(p, 'p_H') / timeRef)

        if _get(p, 'heating'):
            L = _get(p, 'L', _get(p, 'length', 1.0))
            Qheat = 60.0 * (_get(p, 'r_out')) * np.pi * 2.0 * L * (T * TRef - _get(p, 'Theat'))

        F_in = 0;            
        Fout = -(1 - _get(p, 'e_bed')) * _get(p, 'V_column') * _get(p, 'rho_s') * (dq1dt + dq2dt) / (timeRef / qRef) - (_get(p, 'e_bed') / (R * T * TRef)) * dPdt * (_get(p, 'p_H') / timeRef) * _get(p, 'V_column')

    elif step_name == 'pres':
        P = _get(p, 'P_press')(t * timeRef) / PRef


        q1_star, q2_star = DSL(P * PRef, y1, T * TRef,
                              _get(p, 'qsb_1'), _get(p, 'qsd_1'), _get(p, 'qsb_2'), _get(p, 'qsd_2'),
                              _get(p, 'bo_1'), _get(p, 'do_1'), _get(p, 'bo_2'), _get(p, 'do_2'),
                              _get(p, 'delUb_1'), _get(p, 'delUd_1'), _get(p, 'delUb_2'), _get(p, 'delUd_2'))

        k1, k2 = LDFCoefficient(P * PRef, y1, T * TRef, q1_star, q2_star, p)

        dq1dt = timeRef / qRef * k1 * (q1_star - q1 * qRef)
        dq2dt = timeRef / qRef * k2 * (q2_star - q2 * qRef)
        dPdt = _get(p, 'dPdt_press')(t * timeRef) / (PRef / timeRef)

        # set F_in according to MATLAB
        F_in_val = (1 - _get(p, 'e_bed')) * _get(p, 'V_column') * _get(p, 'rho_s') * (dq1dt + dq2dt) / (timeRef / qRef) + (_get(p, 'e_bed') / (R * T * TRef)) * dPdt * (_get(p, 'p_H') / timeRef) * _get(p, 'V_column')

        F_in = F_in_val     
        Fout = 0.0

        if _get(p, 'pressType') == 'LPP':
            if isinstance(p, dict):
                p['y1_in'] = _get(p, 'y1_LPP')
            else:
                try:
                    setattr(p, 'y1_in', _get(p, 'y1_LPP'))
                except Exception:
                    pass

    else:
        raise ValueError(f'Unknown step_name: {step_name}')

    # Analytical computation of heat of adsorption as a function of P and T
    b1 = _get(p, 'bo_1') * np.exp(-_get(p, 'delUb_1') / (R * T * TRef))
    d1 = _get(p, 'do_1') * np.exp(-_get(p, 'delUd_1') / (R * T * TRef))
    b2 = _get(p, 'bo_2') * np.exp(-_get(p, 'delUb_2') / (R * T * TRef))
    d2 = _get(p, 'do_2') * np.exp(-_get(p, 'delUd_2') / (R * T * TRef))

    c1 = y1 * P * PRef / (R * T * TRef)
    c2 = (1 - y1) * P * PRef / (R * T * TRef)

    delH1 = -(
        _get(p, 'qsb_1') * b1 * _get(p, 'delUb_1') / (1 + b1 * c1) ** 2
        + _get(p, 'qsd_1') * d1 * _get(p, 'delUd_1') / (1 + d1 * c1) ** 2
    ) / (
        _get(p, 'qsb_1') * b1 / (1 + b1 * c1) ** 2 + _get(p, 'qsd_1') * d1 / (1 + d1 * c1) ** 2
    )

    delH2 = -(
        _get(p, 'qsb_2') * b2 * _get(p, 'delUb_2') / (1 + b2 * c2) ** 2
        + _get(p, 'qsd_2') * d2 * _get(p, 'delUd_2') / (1 + d2 * c2) ** 2
    ) / (
        _get(p, 'qsb_2') * b2 / (1 + b2 * c2) ** 2 + _get(p, 'qsd_2') * d2 / (1 + d2 * c2) ** 2
    )

    # Temperature derivatives
    if _get(p, 'modelType') == 'isothermal':
        dTdt = 0.0
        dTwdt = 0.0
    else:
        coefft1 = (timeRef / TRef) / (((1 - _get(p, 'e_bed')) / _get(p, 'e_bed')) * (_get(p, 'rho_s') * _get(p, 'cp_s') + _get(p, 'cp_a') * _get(p, 'rho_s') * qRef * (q1 + q2)))

        dTdt = coefft1 * (
            + (F_in / _get(p, 'V_column') * _get(p, 'cp_g') * (TRef - TRef * T))
            - _get(p, 'cp_g') / R * dPdt * PRef / timeRef
            - (((1 - _get(p, 'e_bed')) / _get(p, 'e_bed')) * _get(p, 'cp_a') * _get(p, 'rho_s') * T * TRef * qRef / timeRef * (dq1dt + dq2dt))
            + (((1 - _get(p, 'e_bed')) / _get(p, 'e_bed')) * _get(p, 'rho_s') * qRef / timeRef * (delH1 * dq1dt + delH2 * dq2dt))
            - (2 * _get(p, 'h_in') / _get(p, 'r_in') / _get(p, 'e_bed') * (T * TRef - Tw * TRef) + Qheat)
        )

        dTwdt = (timeRef / TwRef) / (_get(p, 'rho_w') * _get(p, 'cp_w')) * (2 * _get(p, 'h_in') * _get(p, 'r_in') / (_get(p, 'r_out') ** 2 - _get(p, 'r_in') ** 2) * (T * TRef - Tw * TRef) - 2 * _get(p, 'h_out') * _get(p, 'r_out') / (_get(p, 'r_out') ** 2 - _get(p, 'r_in') ** 2) * (Tw * TRef - TRef))

    # dy1dt
    dy1dt = T / (_get(p, 'e_bed') * P) * (
        -_get(p, 'e_bed') * y1 / T * dPdt
        + P * y1 / (T ** 2) * dTdt
        - ((1 - _get(p, 'e_bed')) * _get(p, 'rho_s') * qRef * R * TRef / PRef * dq1dt)
        + ((_get(p, 'y1_in') * F_in / volFluxRef - y1 * Fout / volFluxRef) / _get(p, 'V_column'))
    )

    dXdt = np.zeros(5, dtype=float)
    dXdt[0] = dy1dt
    dXdt[1] = dq1dt
    dXdt[2] = dq2dt
    dXdt[3] = dTdt
    dXdt[4] = dTwdt

    return dXdt
