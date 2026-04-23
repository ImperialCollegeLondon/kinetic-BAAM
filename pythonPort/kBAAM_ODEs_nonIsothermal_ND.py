"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Imperial College London, United Kingdom
Multiphase Systems Laboratory
Year:     2025
MATLAB:   R2024a
Authors:  Hassan Azzan (HA), Ayca Yilmaz (AY)

Purpose:
Evaluates the RHS of the dimensionless ODEs for the non-isothermal k-BAAM with pressure
drop. Takes dimensionless time t, state vector X, parameters struct, cycle step name,
and returns dXdt. Handles four PSA steps (ads, blo, evac, pres) via a switch-style
branch. Supports both isothermal and non-isothermal models.

The DAE system is written in semi-explicit mass-matrix form M*dX/dt = f. The function
solves this directly by computing M^-1 f at each call for the 6-state system.

State vector X (all dimensionless):
    X(1) = y1       mole fraction of component 1 (CO2) [-]
    X(2) = q1/qRef  adsorbed amount of component 1     [-]
    X(3) = q2/qRef  adsorbed amount of component 2     [-]
    X(4) = T/TRef   gas and solid temperature          [-]
    X(5) = Tw/TwRef wall temperature                   [-]
    X(6) = P/PRef   total pressure                     [-]

Equilibrium model: DSL (default) or SSLSTA (parameters.SSLSTA = 1)
Kinetics: Linear Driving Force (LDF) via LDFCoefficient.m
Flow: Darcy's law with linear pressure profile along bed length
Pressure drop: Ergun-based Darcy permeability precomputed in Outputs

Last modified:
- 2026-04-23, HA: Add section headers, inline comments, precomputed constants,
                                    state clamping, dimensional caching, fewer LDF calls in ads step,
                                    remove buildMassMatrix / mode-switching infrastructure
- 2025-10-09, HA: Add wall energy balance
- 2025-09-21, HA: Initial creation

Input arguments:
    - t: dimensionless time [-]
    - X: 6x1 vector of dimensionless state variables
    - parameters: struct of adsorbent properties and process parameters
    - stepName: current cycle step ('ads', 'blo', 'evac', 'pres')

Output arguments:
    - dXdt: 6x1 vector of time derivatives

Dependencies:
    - DSL.m
    - SSLSTA.m
    - LDFCoefficient.m
    - computeDSLHeatUnary.m
    - computeSSLSTAHeatBinaryBT.m

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


def _set(p, key, value):
    if isinstance(p, dict):
        p[key] = value
    else:
        setattr(p, key, value)


def _dsl_heats(parameters, y1, t_nd, p_nd):
    """Analytical DSL heat expression consistent with MATLAB path for DSL mode."""
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

    del_h1 = num1 / den1
    del_h2 = num2 / den2
    return del_h1, del_h2


def _equilibrium_loadings(parameters, p_nd, y1, t_nd):
    dsl = _get(parameters, 'DSL') or _DSL_impl
    if dsl is None:
        raise KeyError('DSL implementation is required in parameters or imports.')
    pref = _get(parameters, 'PRef')
    tref = _get(parameters, 'TRef')
    return dsl(
        p_nd * pref,
        y1,
        t_nd * tref,
        _get(parameters, 'qsb_1'),
        _get(parameters, 'qsd_1'),
        _get(parameters, 'qsb_2'),
        _get(parameters, 'qsd_2'),
        _get(parameters, 'bo_1'),
        _get(parameters, 'do_1'),
        _get(parameters, 'bo_2'),
        _get(parameters, 'do_2'),
        _get(parameters, 'delUb_1'),
        _get(parameters, 'delUd_1'),
        _get(parameters, 'delUb_2'),
        _get(parameters, 'delUd_2'),
    )


def kbaam_odes_nonisothermal_nd(t, x, parameters, step_name: str):
    """Evaluate dX/dt by solving M\f for MATLAB-equivalent pressure-drop system."""
    x = np.asarray(x, dtype=float)
    if x.size < 6:
        raise ValueError('State must be length 6: [y1,q1,q2,T,Tw,P].')

    y1 = np.clip(x[0], 0.0, 1.0)
    q1 = x[1]
    q2 = x[2]
    t_nd = max(x[3], 1e-9)
    tw_nd = x[4]
    p_nd = max(x[5], 1e-9)

    r = 8.3145

    tref = _get(parameters, 'timeRef')
    qref = _get(parameters, 'qRef')
    t_ref = _get(parameters, 'TRef')
    tw_ref = _get(parameters, 'TwRef')
    p_ref = _get(parameters, 'PRef')
    ab = _get(parameters, 'Ab')
    coeff_q = _get(parameters, 'Ab_rhos_qRef_tRef')
    cp_a = _get(parameters, 'cp_a')
    cp_g = _get(parameters, 'cp_g')
    ebed = _get(parameters, 'e_bed')
    vcol = _get(parameters, 'V_column')
    col_l = _get(parameters, 'L')
    area = _get(parameters, 'A_in')
    t_ref_qref = _get(parameters, 'tRef_qRef')
    darcy_k = _get(parameters, 'darcyK')
    cpg_e_v = _get(parameters, 'cpg_eV')
    two_hin = _get(parameters, 'two_hin_rin_e')

    p_dim = p_nd * p_ref
    t_dim = t_nd * t_ref

    qheat = 0.0

    ldf = _get(parameters, 'LDFCoefficient') or _LDF_impl
    if ldf is None:
        raise KeyError('LDFCoefficient implementation is required in parameters or imports.')

    if step_name == 'ads':
        p_out = _get(parameters, 'P_ads')(t * tref) / p_ref

        if _get(parameters, 'cCSTR', False):
            q1_star, q2_star = _equilibrium_loadings(parameters, p_nd, y1, t_nd)
            k1, k2 = ldf(p_dim, y1, t_dim, q1_star, q2_star, parameters)
            dq1dt = t_ref_qref * k1 * (q1_star - q1 * qref)
            dq2dt = t_ref_qref * k2 * (q2_star - q2 * qref)
        else:
            q1_in, q2_in = _equilibrium_loadings(parameters, p_nd, _get(parameters, 'y1_in'), t_nd)
            q1_star, q2_star = _equilibrium_loadings(parameters, p_nd, y1, t_nd)
            k1_in, _ = ldf(p_dim, y1, t_dim, q1_in, q2_in, parameters)
            _, k2 = ldf(p_dim, y1, t_dim, q1_star, q2_star, parameters)
            dq1dt = t_ref_qref * k1_in * (q1_in - q1 * qref)
            dq2dt = t_ref_qref * k2 * (q2_star - q2 * qref)

        f_in = _get(parameters, 'volFlowin') * (2.0 * p_nd - p_out) * p_ref / (r * np.mean([t_ref, t_dim]))
        y1_in = _get(parameters, 'y1_in')
        v_out = (2.0 / col_l) * darcy_k * (p_nd - p_out) * p_ref
        f_out = p_out * p_ref * area * ebed / (r * t_nd * t_ref) * v_out

    elif step_name == 'blo':
        p_out = _get(parameters, 'P_blo')(t * tref) / p_ref
        q1_star, q2_star = _equilibrium_loadings(parameters, p_nd, y1, t_nd)
        k1, k2 = ldf(
            p_dim,
            _get(parameters, 'y1init', y1),
            t_dim,
            _get(parameters, 'q1init', q1_star),
            _get(parameters, 'q2init', q2_star),
            parameters,
        )
        dq1dt = t_ref_qref * k1 * (q1_star - q1 * qref)
        dq2dt = t_ref_qref * k2 * (q2_star - q2 * qref)
        f_in = 0.0
        y1_in = 0.0
        v_out = (2.0 / col_l) * darcy_k * (p_nd - p_out) * p_ref
        f_out = p_dim * area * ebed / (r * t_dim) * v_out

    elif step_name == 'evac':
        p_out = _get(parameters, 'P_evac')(t * tref) / p_ref
        q1_star, q2_star = _equilibrium_loadings(parameters, p_nd, y1, t_nd)
        k1, k2 = ldf(p_dim, y1, t_dim, q1_star, q2_star, parameters)
        dq1dt = t_ref_qref * k1 * (q1_star - q1 * qref)
        dq2dt = t_ref_qref * k2 * (q2_star - q2 * qref)
        if _get(parameters, 'heating', False) and t_dim < _get(parameters, 'Theat', 0.0):
            qheat = (
                _get(parameters, 'heatPowerDensity', 0.0)
                * (_get(parameters, 'Theat', 0.0) - t_dim)
                / max(_get(parameters, 'Theat', 1.0) - t_ref, 1e-12)
                / max(_get(parameters, 'r_out') - _get(parameters, 'r_in'), 1e-12)
            )
        f_in = 0.0
        y1_in = 0.0
        v_out = (2.0 / col_l) * darcy_k * (p_nd - p_out) * p_ref
        f_out = p_dim * area * ebed / (r * t_dim) * v_out

    elif step_name == 'pres':
        p_out = _get(parameters, 'P_press')(t * tref) / p_ref
        q1_star, q2_star = _equilibrium_loadings(parameters, p_nd, y1, t_nd)
        k1, k2 = ldf(p_dim, y1, t_dim, q1_star, q2_star, parameters)
        dq1dt = t_ref_qref * k1 * (q1_star - q1 * qref)
        dq2dt = t_ref_qref * k2 * (q2_star - q2 * qref)
        v_in = (2.0 / col_l) * darcy_k * (p_out - p_nd) * p_ref
        f_in = p_out * p_ref * area * ebed / (r * t_ref) * v_in
        f_out = 0.0
        if _get(parameters, 'pressType') == 'LPP':
            y1_in = _get(parameters, 'y1_LPP', _get(parameters, 'y1_in'))
        else:
            y1_in = _get(parameters, 'y1_in')
    else:
        raise ValueError(f'Unknown step_name: {step_name}')

    f = np.zeros(6, dtype=float)
    f[0] = r * t_dim / p_dim * (y1_in * f_in - y1 * f_out) / (ebed * vcol)
    f[1] = dq1dt
    f[2] = dq2dt

    is_isothermal = _get(parameters, 'isIsothermal', _get(parameters, 'modelType') == 'isothermal')
    if not is_isothermal:
        f[3] = cpg_e_v * (f_in * t_ref - f_out * t_dim) - two_hin * (t_dim - tw_nd * tw_ref)
        f[4] = _get(parameters, 'wall_prefactor') * (
            _get(parameters, 'wall_coeff1') * (t_dim - tw_nd * tw_ref)
            - _get(parameters, 'wall_coeff2') * (tw_nd * tw_ref - t_ref)
            + qheat
        )

    f[5] = r * t_dim * (f_in - f_out) / (ebed * vcol)

    m = np.zeros((6, 6), dtype=float)
    rt_pp = r * t_dim / p_dim

    m[0, 0] = 1.0 / tref
    m[0, 1] = rt_pp * coeff_q
    m[0, 3] = -y1 / (t_nd * tref)
    m[0, 5] = y1 / (p_nd * tref)

    m[1, 1] = 1.0
    m[2, 2] = 1.0

    if not is_isothermal:
        ceff = ab * (_get(parameters, 'rho_s') * _get(parameters, 'cp_s') + cp_a * _get(parameters, 'rho_s') * qref * (q1 + q2))
        del_h1, del_h2 = _dsl_heats(parameters, y1, t_nd, p_nd)
        is_resin = _get(parameters, 'isResin', _get(parameters, 'processType') in ('Resin', 'ResinSens'))
        if is_resin:
            del_h1 = -_get(parameters, 'delUb_1')
            del_h2 = -_get(parameters, 'delUb_2')

        m[3, 1] = -coeff_q * (del_h1 - cp_a * t_dim)
        m[3, 2] = -coeff_q * (del_h2 - cp_a * t_dim)
        m[3, 3] = ceff * t_ref / tref
        m[3, 5] = cp_g / r * p_ref / tref
    else:
        m[3, 3] = 1.0

    m[4, 4] = 1.0

    m[5, 1] = r * t_dim * coeff_q
    m[5, 2] = r * t_dim * coeff_q
    m[5, 3] = -p_dim / (t_nd * tref)
    m[5, 5] = p_ref / tref

    if is_isothermal:
        m[0, 3] = 0.0
        m[5, 3] = 0.0

    d_x = np.linalg.solve(m, f)
    return d_x
