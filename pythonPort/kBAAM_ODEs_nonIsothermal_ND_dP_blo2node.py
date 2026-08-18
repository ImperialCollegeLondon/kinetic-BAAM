"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Imperial College London, United Kingdom
Multiphase Systems Laboratory
Year:     2026
MATLAB:   R2024a
Authors:  Hassan Azzan (HA)

Purpose:
Evaluates the RHS of the dimensionless ODEs for the blowdown step of the non-isothermal
k-BAAM using TWO spatial nodes. P, T, and Tw are shared between the two nodes. Only
composition (y1) and adsorbed amounts (q1, q2) differ between nodes.

State vector X (10 states, all dimensionless):
    X[0]  = y1_1       mole fraction of CO2 in node 1 (feed end)   [-]
    X[1]  = q1_1/qRef  adsorbed CO2 in node 1                       [-]
    X[2]  = q2_1/qRef  adsorbed N2  in node 1                       [-]
    X[3]  = y1_2       mole fraction of CO2 in node 2 (product end) [-]
    X[4]  = q1_2/qRef  adsorbed CO2 in node 2                       [-]
    X[5]  = q2_2/qRef  adsorbed N2  in node 2                       [-]
    X[6]  = T/TRef     shared gas/solid temperature                  [-]
    X[7]  = Tw/TwRef   shared wall temperature                       [-]
    X[8]  = P1/PRef    node 1 (feed end) pressure                    [-]
    X[9]  = P2/PRef    node 2 (product end) pressure                 [-]

Flow model (nodes sized by adsorption loading fraction, independent pressures):
    Node 1 (feed end)    occupies fraction f = loadingFraction of column volume.
    Node 2 (product end) occupies fraction (1-f) of column volume.
    Darcy velocity node 1 -> node 2 (node 1 centre to node 2 centre = L/2):
        v_12  = (2/L) * darcyK * (P1 - P2) * PRef   [m/s]
    Darcy velocity node 2 -> outlet (node 2 centre to outlet = L2/2):
        v_out = (2/L2) * darcyK * (P2 - P_blo) * PRef   [m/s]

Last modified:
- 2026-05-01, HA: Initial creation (Python port)

Dependencies:
    - DSL.py / SSLSTA.py
    - LDFCoefficient.py

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"""
from __future__ import annotations

import numpy as np

try:
    from .DSL import DSL as _DSL_impl
    from .LDFCoefficient import LDFCoefficient as _LDF_impl
except Exception:
    try:
        from DSL import DSL as _DSL_impl
        from LDFCoefficient import LDFCoefficient as _LDF_impl
    except Exception:
        _DSL_impl = None
        _LDF_impl = None


def _get(p, key, default=None):
    if isinstance(p, dict):
        return p.get(key, default)
    return getattr(p, key, default)


def _eq(parameters, p_nd, y1, t_nd):
    """Return equilibrium loadings (q1_star, q2_star) at dimensionless (P, y1, T)."""
    pref = _get(parameters, 'PRef')
    tref = _get(parameters, 'TRef')
    if _get(parameters, 'SSLSTA', False):
        sslsta_fn = _get(parameters, 'SSLSTAFunction') or _get(parameters, 'SSLSTA_fn')
        if sslsta_fn is None:
            try:
                from SSLSTA import SSLSTA as sslsta_fn
            except Exception:
                sslsta_fn = None
        if sslsta_fn is None:
            raise KeyError("SSLSTA implementation required when parameters['SSLSTA'] is True.")
        return sslsta_fn(p_nd * pref, y1, t_nd * tref, parameters)
    dsl = _get(parameters, 'DSL') or _DSL_impl
    if dsl is None:
        raise KeyError('DSL implementation is required.')
    return dsl(
        p_nd * pref, y1, t_nd * tref,
        _get(parameters, 'qsb_1'), _get(parameters, 'qsd_1'),
        _get(parameters, 'qsb_2'), _get(parameters, 'qsd_2'),
        _get(parameters, 'bo_1'),  _get(parameters, 'do_1'),
        _get(parameters, 'bo_2'),  _get(parameters, 'do_2'),
        _get(parameters, 'delUb_1'), _get(parameters, 'delUd_1'),
        _get(parameters, 'delUb_2'), _get(parameters, 'delUd_2'),
    )


def _dsl_heats(parameters, y1, t_nd, p_nd):
    """Analytical DSL isosteric heat (dimensionless T and P inputs)."""
    r = 8.3145
    tref = _get(parameters, 'TRef')
    pref = _get(parameters, 'PRef')
    tdim = t_nd * tref
    c1 = y1 * p_nd * pref / (r * tdim)
    c2 = (1.0 - y1) * p_nd * pref / (r * tdim)
    b1 = _get(parameters, 'bo_1') * np.exp(-_get(parameters, 'delUb_1') / (r * tdim))
    d1 = _get(parameters, 'do_1') * np.exp(-_get(parameters, 'delUd_1') / (r * tdim))
    b2 = _get(parameters, 'bo_2') * np.exp(-_get(parameters, 'delUb_2') / (r * tdim))
    d2 = _get(parameters, 'do_2') * np.exp(-_get(parameters, 'delUd_2') / (r * tdim))
    num1 = -(
        _get(parameters, 'qsb_1') * b1 * _get(parameters, 'delUb_1') / (1.0 + b1 * c1) ** 2
        + _get(parameters, 'qsd_1') * d1 * _get(parameters, 'delUd_1') / (1.0 + d1 * c1) ** 2
    )
    den1 = (
        _get(parameters, 'qsb_1') * b1 / (1.0 + b1 * c1) ** 2
        + _get(parameters, 'qsd_1') * d1 / (1.0 + d1 * c1) ** 2
    )
    num2 = -(
        _get(parameters, 'qsb_2') * b2 * _get(parameters, 'delUb_2') / (1.0 + b2 * c2) ** 2
        + _get(parameters, 'qsd_2') * d2 * _get(parameters, 'delUd_2') / (1.0 + d2 * c2) ** 2
    )
    den2 = (
        _get(parameters, 'qsb_2') * b2 / (1.0 + b2 * c2) ** 2
        + _get(parameters, 'qsd_2') * d2 / (1.0 + d2 * c2) ** 2
    )
    return num1 / den1, num2 / den2


def kbaam_odes_nonisothermal_nd_dp_blo2node(t, x, parameters, step_name: str):
    """Evaluate dX/dt for the 10-state 2-node blowdown/evac ODE.

    Parameters
    ----------
    t : float
        Dimensionless time [-].
    x : array-like, length 10
        Dimensionless state vector (see module docstring).
    parameters : dict
        Process and adsorbent parameters (must contain all precomputed constants).
    step_name : str
        'blo' or 'evac'.

    Returns
    -------
    dXdt : ndarray, shape (10,)
    """
    x = np.asarray(x, dtype=float)

    y1_1 = np.clip(x[0], 0.0, 1.0)
    q1_1 = x[1]
    q2_1 = x[2]
    y1_2 = np.clip(x[3], 0.0, 1.0)
    q1_2 = x[4]
    q2_2 = x[5]
    T    = max(x[6], 1e-9)   # shared dimensionless temperature
    Tw   = x[7]               # shared dimensionless wall temperature
    P1   = max(x[8], 1e-9)   # node 1 dimensionless pressure
    P2   = max(x[9], 1e-9)   # node 2 dimensionless pressure

    R = 8.3145

    # Node-size fraction based on loading fraction
    if step_name == 'blo':
        f_node = max(0.05, min(0.95, _get(parameters, 'loadingFraction', 0.5)))
    else:
        f_node = 1.0 - max(0.03, min(0.97, _get(parameters, 'loadingFraction', 0.5)))

    # Unpack parameters
    tRef      = _get(parameters, 'timeRef')
    qRef      = _get(parameters, 'qRef')
    TRef      = _get(parameters, 'TRef')
    TwRef     = _get(parameters, 'TwRef')
    PRef      = _get(parameters, 'PRef')
    Ab        = _get(parameters, 'Ab')
    coeff_q   = _get(parameters, 'Ab_rhos_qRef_tRef')
    cp_a      = _get(parameters, 'cp_a')
    cp_g      = _get(parameters, 'cp_g')
    e         = _get(parameters, 'e_bed')
    V         = _get(parameters, 'V_column')
    L         = _get(parameters, 'L')
    A         = _get(parameters, 'A_in')
    tRef_qRef = _get(parameters, 'tRef_qRef')
    darcyK    = _get(parameters, 'darcyK')

    V_node1 = f_node * V
    V_node2 = (1.0 - f_node) * V

    P1_dim = P1 * PRef
    P2_dim = P2 * PRef
    T_dim  = T  * TRef

    # Prescribed outlet pressure
    if step_name == 'blo':
        P_out_nd = _get(parameters, 'P_blo')(t * tRef) / PRef
    else:
        P_out_nd = _get(parameters, 'P_evac')(t * tRef) / PRef

    # Darcy flows
    L2 = (1.0 - f_node) * L
    v_12  = (2.0 / L)  * darcyK * (P1 - P2)      * PRef   # node 1 -> node 2 [m/s]
    v_out = (2.0 / L2) * darcyK * (P2 - P_out_nd) * PRef   # node 2 -> outlet [m/s]

    F_12  = P1_dim * A * e / (R * T_dim) * v_12   # inter-node molar flow [mol/s]
    F_out = P2_dim * A * e / (R * T_dim) * v_out   # outlet molar flow [mol/s]
    F_in  = 0.0

    ldf = _get(parameters, 'LDFCoefficient') or _LDF_impl
    if ldf is None:
        raise KeyError('LDFCoefficient implementation is required.')

    # LDF kinetics — node 1
    q1_star1, q2_star1 = _eq(parameters, P1, y1_1, T)
    k1_1, k2_1 = ldf(P1_dim, y1_1, T_dim, q1_star1, q2_star1, parameters)
    dq1dt_1 = tRef_qRef * k1_1 * (q1_star1 - q1_1 * qRef)
    dq2dt_1 = tRef_qRef * k2_1 * (q2_star1 - q2_1 * qRef)

    # LDF kinetics — node 2
    q1_star2, q2_star2 = _eq(parameters, P2, y1_2, T)
    k1_2, k2_2 = ldf(P2_dim, y1_2, T_dim, q1_star2, q2_star2, parameters)
    dq1dt_2 = tRef_qRef * k1_2 * (q1_star2 - q1_2 * qRef)
    dq2dt_2 = tRef_qRef * k2_2 * (q2_star2 - q2_2 * qRef)

    # Assemble RHS
    f = np.zeros(10, dtype=float)

    # Node 1: y1 balance (no inflow, convective outflow F_12)
    f[0] = R * T_dim / P1_dim * (0.0 - y1_1 * F_12) / (e * V_node1)
    f[1] = dq1dt_1
    f[2] = dq2dt_1

    # Node 2: y1 balance (inflow F_12 at y1_1, outflow F_out at y1_2)
    f[3] = R * T_dim / P2_dim * (y1_1 * F_12 - y1_2 * F_out) / (e * V_node2)
    f[4] = dq1dt_2
    f[5] = dq2dt_2

    cpg_eV_total = cp_g / (V * e)
    two_hin = _get(parameters, 'two_hin_rin_e')
    Qheat = 0.0

    is_isothermal = _get(parameters, 'isIsothermal', _get(parameters, 'modelType') == 'isothermal')

    if not is_isothermal:
        # Shared energy balance (whole column, F_out exits at node-2 temperature)
        f[6] = cpg_eV_total * (0.0 - F_out * T_dim) - two_hin * (T_dim - Tw * TwRef)
        # Shared wall energy balance
        f[7] = _get(parameters, 'wall_prefactor') * (
            _get(parameters, 'wall_coeff1') * (T_dim - Tw * TwRef)
            - _get(parameters, 'wall_coeff2') * (Tw * TwRef - TRef)
            + Qheat
        )

    # Per-node overall material balance
    f[8] = R * T_dim * (F_in - F_12)  / (e * V_node1)
    f[9] = R * T_dim * (F_12 - F_out) / (e * V_node2)

    # Build 10x10 mass matrix and solve
    M = _build_mass_matrix_10(
        y1_1, q1_1, q2_1, y1_2, q1_2, q2_2, T, Tw, P1, P2,
        parameters, R, tRef, qRef, TRef, TwRef, PRef, Ab, coeff_q, cp_a, cp_g,
        f_node, V_node1, V_node2, is_isothermal,
    )
    dXdt = np.linalg.solve(M, f)
    return dXdt


def _build_mass_matrix_10(y1_1, q1_1, q2_1, y1_2, q1_2, q2_2, T, Tw, P1, P2,
                          parameters, R, tRef, qRef, TRef, TwRef, PRef, Ab, coeff_q,
                          cp_a, cp_g, f_node, V_node1, V_node2, is_isothermal):
    """Assemble the 10x10 mass matrix for the 2-node system."""
    P1_dim = P1 * PRef
    P2_dim = P2 * PRef
    T_dim  = T  * TRef
    RT_PP1 = R * T_dim / P1_dim
    RT_PP2 = R * T_dim / P2_dim

    M = np.zeros((10, 10), dtype=float)

    # Row 0: y1_1 species balance
    M[0, 0] = 1.0 / tRef
    M[0, 1] = RT_PP1 * coeff_q
    M[0, 6] = -y1_1 / (T * tRef)
    M[0, 8] =  y1_1 / (P1 * tRef)

    # Row 1: q1_1 LDF
    M[1, 1] = 1.0

    # Row 2: q2_1 LDF
    M[2, 2] = 1.0

    # Row 3: y1_2 species balance
    M[3, 3] = 1.0 / tRef
    M[3, 4] = RT_PP2 * coeff_q
    M[3, 6] = -y1_2 / (T * tRef)
    M[3, 9] =  y1_2 / (P2 * tRef)

    # Row 4: q1_2 LDF
    M[4, 4] = 1.0

    # Row 5: q2_2 LDF
    M[5, 5] = 1.0

    if not is_isothermal:
        # Volume-weighted heat capacity and heats of adsorption
        q1_avg = f_node * q1_1 + (1.0 - f_node) * q1_2
        q2_avg = f_node * q2_1 + (1.0 - f_node) * q2_2
        Ceff = Ab * (parameters['rho_s'] * parameters['cp_s']
                     + cp_a * parameters['rho_s'] * qRef * (q1_avg + q2_avg))

        if _get(parameters, 'SSLSTA', False):
            y1_avg   = f_node * y1_1 + (1.0 - f_node) * y1_2
            P_avg_nd = f_node * P1   + (1.0 - f_node) * P2
            try:
                from computeSSLSTAHeatBinaryBT import computeSSLSTAHeatBinaryBT
                delH1_1, delH2_1 = computeSSLSTAHeatBinaryBT(
                    P_avg_nd * PRef / 1e5, y1_avg, T_dim,
                    [parameters['SSLSTA1'], parameters['SSLSTA2']],
                )
            except Exception:
                delH1_1 = -_get(parameters, 'delUb_1', 0.0)
                delH2_1 = -_get(parameters, 'delUb_2', 0.0)
            delH1_2 = delH1_1
            delH2_2 = delH2_1
        else:
            delH1_1, delH2_1 = _dsl_heats(parameters, y1_1, T, P1)
            delH1_2, delH2_2 = _dsl_heats(parameters, y1_2, T, P2)

        if _get(parameters, 'isResin', False):
            delH1_1 = -_get(parameters, 'delUb_1', 0.0)
            delH2_1 = -_get(parameters, 'delUb_2', 0.0)
            delH1_2 = delH1_1
            delH2_2 = delH2_1

        # Row 6: shared energy balance
        M[6, 1] = -coeff_q * f_node       * (delH1_1 - cp_a * T_dim)
        M[6, 2] = -coeff_q * f_node       * (delH2_1 - cp_a * T_dim)
        M[6, 4] = -coeff_q * (1.0-f_node) * (delH1_2 - cp_a * T_dim)
        M[6, 5] = -coeff_q * (1.0-f_node) * (delH2_2 - cp_a * T_dim)
        M[6, 6] = Ceff * TRef / tRef
        M[6, 8] = f_node       * cp_g / R * PRef / tRef
        M[6, 9] = (1.0-f_node) * cp_g / R * PRef / tRef

    # Row 7: wall energy balance (always identity in M — rhs handles dynamics)
    M[7, 7] = 1.0

    # Row 8: node 1 overall material balance
    M[8, 1] = R * T_dim * coeff_q
    M[8, 2] = R * T_dim * coeff_q
    M[8, 6] = -P1_dim / (T * tRef)
    M[8, 8] = PRef / tRef

    # Row 9: node 2 overall material balance
    M[9, 4] = R * T_dim * coeff_q
    M[9, 5] = R * T_dim * coeff_q
    M[9, 6] = -P2_dim / (T * tRef)
    M[9, 9] = PRef / tRef

    # Isothermal simplification: remove T couplings, replace T/Tw rows with identity
    if is_isothermal:
        M[0, 6] = 0.0
        M[3, 6] = 0.0
        M[6, :] = 0.0;  M[6, 6] = 1.0
        # Row 7 already has only M[7,7]=1 (no other entries), so stays as identity
        M[8, 6] = 0.0
        M[9, 6] = 0.0

    return M
