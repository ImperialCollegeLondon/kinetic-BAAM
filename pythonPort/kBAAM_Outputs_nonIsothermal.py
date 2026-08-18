"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Imperial College London, United Kingdom
Multiphase Systems Laboratory
Year:     2025
MATLAB:   R2024a
Authors:  Hassan Azzan (HA)

Purpose:
Simulates a 4-step PSA/PVSA cycle (adsorption, blowdown, evacuation, pressurisation)
using the non-isothermal 0-D kinetic batch adsorber analogue model (k-BAAM) with
pressure drop, and returns key performance indicators (KPIs): CO2 purity, recovery,
productivity, and specific energy consumption.

The cycle is iterated until cyclic steady state (CSS) is reached (max 200 cycles).
CSS is detected by convergence of the bed-averaged temperature and loading between
successive cycles. Step pressures follow exponential profiles parameterised by lambda.

ODE solver: ode15s in MATLAB; solve_ivp in this Python port.
State vector (dimensionless): [y1, q1/qRef, q2/qRef, T/TRef, Tw/TwRef, P/PRef]
Reference scales: timeRef = V/Qin, TRef = T_feed, PRef = p_H, qRef = qsb_1+qsd_1

Step-invariant composite constants (e.g. Ab, darcyK, cpg_eV, wall coefficients)
are precomputed once before the cycle loop and stored in the parameters dict to
avoid redundant arithmetic inside the ODE function.

Supports: VSA, PVSA, TVSA, PTVSA, multiple adsorbent types (DSL / SSLSTA isotherms),
isothermal and non-isothermal modes, LPP pressurisation, CSTR mixing assumption,
and breakthrough test mode.

Last modified:
- 2026-04-23, HA: Add section headers, inline comments, precompute composite ODE
                                    constants, remove dead variables, deduplicate rawData write,
                                    move odeset outside cycle loop, remove redundant DSL call,
                                    remove useMassMatrix flag (mass-matrix mode retired)
- 2025-12-18, HA: Add total material balance and pressure drop
- 2025-10-09, HA: Add wall energy balance
- 2025-10-08, HA: Add reverse engineering optimisation method
- 2025-09-17, HA: Initial creation

Input arguments:
    - parameters: struct of adsorbent properties and process parameters
    - thetaIn: optional theta vector of optimisation decision variables

Output arguments:
    - KPIs: process indicators or objective values, depending on output mode

Dependencies:
    - kBAAM_ODEs_nonIsothermal_ND_dP.m
    - DSL.m
    - SSLSTA.m
    - LDFCoefficient.m
    - computeDSLHeatUnary.m
    - computeSSLSTAHeatBinaryBT.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"""
from __future__ import annotations
import numpy as np
from typing import Any, Dict, Sequence, Optional
import warnings

from scipy.integrate import solve_ivp

from kBAAM_ODEs_nonIsothermal_ND import kbaam_odes_nonisothermal_nd
from kBAAM_ODEs_nonIsothermal_ND_dP_blo2node import kbaam_odes_nonisothermal_nd_dp_blo2node
from createParameters import create_parameters
import pdb
try:
    import matplotlib_inline as _mpl_inline
except ImportError:
    _mpl_inline = None
import os
from datetime import datetime
def kbaam_outputs_nonisothermal(parameters: Dict[str, Any], thetaIn: Optional[Sequence[float]] = None, raise_on_error: bool = False, solver_method: str = 'BDF'):
    """Run cyclic steady-state simulation and return KPIs.

    Parameters follows the dict returned by create_parameters(). If theta is
    provided, it will be used to override selected parameters depending on
    processType (matches MATLAB behaviour).

    Returns:
        KPIs: 2-element list/array with the KPIs (matching MATLAB output)
    """
    if _mpl_inline is not None:
        try:
            _mpl_inline.backend_inline.set_matplotlib_formats('png', 'jpeg')
        except Exception:
            pass
    
    startTime = datetime.now()
    

    try:
        # Handle theta input similarly to the MATLAB script
        if thetaIn is not None:
            # If parameters.outputType == 'opt' the MATLAB code multiplies
            # theta by xRef. We skip that unless xRef present.
            thetavec = np.asarray(thetaIn)
        else:
            thetavec = None

        if thetavec is not None:
            pt = parameters.get('processType', 'PVSA')
            if pt == 'VSA':
                parameters['v_in'] = thetavec[0]
                parameters['p_I'] = thetavec[1]
                parameters['t_ads'] = thetavec[2]
                parameters['t_blo'] = thetavec[3]
                parameters['t_evac'] = thetavec[4]
                parameters['p_L'] = thetavec[5]
                parameters['p_H'] = 101325.0
            elif pt == 'PVSA':
                parameters['v_in'] = thetavec[0]
                parameters['p_I'] = thetavec[1]
                parameters['t_ads'] = thetavec[2]
                parameters['t_blo'] = thetavec[3]
                parameters['t_evac'] = thetavec[4]
                if len(thetavec) > 6:
                    parameters['p_L'] = thetavec[5]
                    parameters['p_H'] = thetavec[6]
                else:
                    parameters['p_H'] = thetavec[5]
            elif pt == 'AdsorbentVSA':
                parameters['v_in'] = 0.4
                parameters['p_I'] = thetavec[0]
                parameters['t_ads'] = thetavec[1]
                parameters['t_blo'] = thetavec[2]
                parameters['t_evac'] = 800.0
                parameters['qsb_1'] = thetavec[3]
                parameters['qsb_2'] = thetavec[3]
                parameters['qsd_1'] = 0.0
                parameters['qsd_2'] = 0.0
                parameters['bo_1'] = thetavec[4]
                parameters['bo_2'] = thetavec[4]
                parameters['do_1'] = 0.0
                parameters['do_2'] = 0.0
                parameters['delUb_1'] = -thetavec[5]
                parameters['delUb_2'] = -thetavec[6]
                parameters['delUd_1'] = 0.0
                parameters['delUd_2'] = 0.0
                if parameters.get('outputType') == 'opt':
                    parameters['p_I'] = 10 ** parameters['p_I']
            elif pt == 'AdsorbentPVSA':
                parameters['v_in'] = 0.4
                parameters['p_I'] = thetavec[0]
                parameters['t_ads'] = thetavec[1]
                parameters['t_blo'] = thetavec[2]
                parameters['t_evac'] = 2000.0
                parameters['qsb_1'] = thetavec[3]
                parameters['qsb_2'] = thetavec[3]
                parameters['qsd_1'] = 0.0
                parameters['qsd_2'] = 0.0
                parameters['bo_1'] = thetavec[4]
                parameters['bo_2'] = thetavec[4]
                parameters['do_1'] = 0.0
                parameters['do_2'] = 0.0
                parameters['delUb_1'] = -thetavec[5]
                parameters['delUb_2'] = -thetavec[6]
                parameters['delUd_1'] = 0.0
                parameters['delUd_2'] = 0.0
                parameters['p_H'] = thetavec[7]
                if parameters.get('outputType') == 'opt':
                    parameters['p_I'] = 10 ** parameters['p_I']
                    parameters['p_H'] = 10 ** parameters['p_H']
            elif pt == 'AdsorbentVSAb0':
                parameters['v_in'] = 0.4
                parameters['p_I'] = thetavec[0]
                parameters['t_ads'] = thetavec[1]
                parameters['t_blo'] = thetavec[2]
                parameters['t_evac'] = 800.0
                parameters['qsb_1'] = thetavec[3]
                parameters['qsb_2'] = thetavec[3]
                parameters['qsd_1'] = 0.0
                parameters['qsd_2'] = 0.0
                parameters['bo_1'] = thetavec[4]
                parameters['bo_2'] = thetavec[5]
                parameters['do_1'] = 0.0
                parameters['do_2'] = 0.0
                parameters['delUb_1'] = -thetavec[6]
                parameters['delUb_2'] = -thetavec[6]
                parameters['delUd_1'] = 0.0
                parameters['delUd_2'] = 0.0
                if parameters.get('outputType') == 'opt':
                    parameters['p_I'] = 10 ** parameters['p_I']
            elif pt == 'AdsorbentPVSAb0':
                parameters['v_in'] = 0.4
                parameters['p_I'] = thetavec[0]
                parameters['t_ads'] = thetavec[1]
                parameters['t_blo'] = thetavec[2]
                parameters['t_evac'] = 2000.0
                parameters['qsb_1'] = thetavec[3]
                parameters['qsb_2'] = thetavec[3]
                parameters['qsd_1'] = 0.0
                parameters['qsd_2'] = 0.0
                parameters['bo_1'] = thetavec[4]
                parameters['bo_2'] = thetavec[5]
                parameters['do_1'] = 0.0
                parameters['do_2'] = 0.0
                parameters['delUb_1'] = -thetavec[6]
                parameters['delUb_2'] = -thetavec[6]
                parameters['delUd_1'] = 0.0
                parameters['delUd_2'] = 0.0
                parameters['p_H'] = thetavec[7]
                if parameters.get('outputType') == 'opt':
                    parameters['p_I'] = 10 ** parameters['p_I']
                    parameters['p_H'] = 10 ** parameters['p_H']
            elif pt == 'Resin':
                parameters['v_in'] = thetavec[5]
                parameters['y1_in'] = 0.0004
                parameters['p_I'] = thetavec[0]
                parameters['t_ads'] = thetavec[1]
                parameters['t_blo'] = thetavec[2]
                parameters['t_evac'] = thetavec[3]
                parameters['p_L'] = thetavec[4]
                parameters['qsb_2'] = 5.3446e-02 * parameters['qsb_1'] / 2.38
                parameters['qsd_1'] = 0.0
                parameters['qsd_2'] = 0.0
                parameters['bo_2'] = 1.0137e-05
                parameters['do_1'] = 0.0
                parameters['do_2'] = 0.0
                parameters['delUb_1'] = -100e3
                parameters['delUb_2'] = -1.3912e4
                parameters['delUd_1'] = 0.0
                parameters['delUd_2'] = 0.0
                parameters['heating'] = True
                parameters['Theat'] = 373.0
                parameters['rho_s'] = 1123.0
                parameters['cp_s'] = 1300.0
                if parameters.get('outputType') == 'opt':
                    parameters['p_I'] = 10 ** parameters['p_I']
                    parameters['p_L'] = 10 ** parameters['p_L']
                    parameters['p_H'] = 1e5
                if not parameters.get('fixResins', False):
                    parameters['LDF'] = 0.132 * np.exp(-2.076 * parameters['qsb_1'])
                    parameters['bo_1'] = 1.06e-16 * parameters['qsb_1']

        Rg = 8.3145

        dt = 0.04 if parameters.get('processType') in ('Resin', 'ResinSens') else 0.05
        t_ads = np.arange(0.0, parameters['t_ads'] + dt, dt)
        t_blo = np.arange(0.0, parameters['t_blo'] + dt, dt)
        t_evac = np.arange(0.0, parameters['t_evac'] + dt, dt)
        t_press = np.arange(0.0, parameters['t_press'] + dt, dt)

        qRef = parameters['qsb_1'] + parameters['qsd_1']
        parameters['qRef'] = qRef
        parameters['refVals'] = np.array([1.0, qRef, qRef, parameters['T_feed'], parameters['T_feed'], parameters['p_H']])

        if parameters.get('heating', False):
            parameters['heatPowerDensity'] = 5e3
            parameters['h_in'] = 8.6 * 3.0
            parameters['Theat'] = parameters.get('Theat', 0.0) + 6.0
        else:
            parameters['Theat'] = 0.0

        # geometry
        parameters['Lbyr'] = parameters.get('Lbyr', 7.0)
        parameters['r_in'] = (parameters['V_column'] / (parameters['Lbyr'] * np.pi)) ** (1.0 / 3.0)
        parameters['r_out'] = parameters['r_in'] + 0.0175
        parameters['v_in'] = parameters['v_in'] / parameters['e_bed']
        parameters['L'] = parameters['Lbyr'] * parameters['r_in']
        parameters['A_in'] = parameters['r_in'] ** 2 * np.pi
        parameters['volFlowin'] = parameters['v_in'] * parameters['A_in'] * parameters['e_bed']
        parameters['deltaP'] = 0.0
        parameters['F_in'] = parameters['volFlowin'] * (parameters['p_H'] + parameters['deltaP']) / (Rg * parameters['T_feed'])
        timeRef = parameters['V_column'] / parameters['volFlowin']
        parameters['timeRef'] = timeRef
        parameters['TRef'] = parameters['T_feed']
        parameters['TwRef'] = parameters['T_feed']
        parameters['PRef'] = parameters['p_H']

        # pressure schedules
        parameters['lambda'] = parameters.get('lambda', 0.5)

        # helper to ensure a time-function accepts scalar or array t and
        # returns a scalar for scalar input or an array shaped like t for array input
        def _ensure_vector_fn(fn_or_val):
            if callable(fn_or_val):
                fn = fn_or_val
                def wrapped(t):
                    t_arr = np.asarray(t)
                    # scalar input -> return scalar
                    if t_arr.shape == ():
                        return fn(float(t))
                    # array input -> try forwarding the array; if fn returns scalar,
                    # broadcast it to the input shape
                    out = fn(t_arr)
                    out_arr = np.asarray(out)
                    if out_arr.shape == ():
                        return np.full(t_arr.shape, out_arr.item())
                    return out_arr
                return wrapped
            else:
                # constant value -> return constant or array of constants
                val = fn_or_val
                def const_fn(t):
                    t_arr = np.asarray(t)
                    if t_arr.shape == ():
                        return val
                    return np.full(t_arr.shape, val)
                return const_fn

        # default schedule expressions (vectorized via numpy)
        P_ads_def = lambda t: parameters['p_H'] + parameters['deltaP']
        P_blo_def = lambda t: parameters['p_I'] + (parameters['P_initH'] - parameters['p_I']) * np.exp(-parameters['lambda'] * t)
        P_evac_def = lambda t: parameters['p_L'] + (parameters['p_I'] - parameters['p_L']) * np.exp(-parameters['lambda'] * t)
        P_press_def = lambda t: parameters['P_initH'] + (parameters['P_initR'] - parameters['P_initH']) * np.exp(-parameters['lambda'] * t)

        dPdt_blo_def = lambda t: -parameters['lambda'] * (parameters['P_initH'] - parameters['p_I']) * np.exp(-parameters['lambda'] * t)
        dPdt_evac_def = lambda t: -parameters['lambda'] * (parameters['P_initL'] - parameters['p_L']) * np.exp(-parameters['lambda'] * t)
        dPdt_press_def = lambda t: -parameters['lambda'] * (parameters['P_initR'] - parameters['P_initH']) * np.exp(-parameters['lambda'] * t)

        # attach defaults or use provided functions/values, ensuring vector behaviour
        parameters['P_ads'] = _ensure_vector_fn(parameters.get('P_ads', P_ads_def))
        parameters['P_initH'] = parameters['P_ads'](parameters['t_ads'])
        parameters['P_blo'] = _ensure_vector_fn(parameters.get('P_blo', P_blo_def))
        parameters['P_initL'] = parameters['P_blo'](parameters['t_blo'])
        parameters['P_evac'] = _ensure_vector_fn(parameters.get('P_evac', P_evac_def))
        parameters['P_initR'] = parameters['P_evac'](parameters['t_evac'])
        parameters['P_press'] = _ensure_vector_fn(parameters.get('P_press', P_press_def))

        parameters['dPdt_blo'] = _ensure_vector_fn(parameters.get('dPdt_blo', dPdt_blo_def))
        parameters['dPdt_evac'] = _ensure_vector_fn(parameters.get('dPdt_evac', dPdt_evac_def))
        parameters['dPdt_press'] = _ensure_vector_fn(parameters.get('dPdt_press', dPdt_press_def))
        
        
        # initial conditions
        y1Init = 0.99
        # call DSL to get q1Init q2Init
        DSL = parameters.get('DSL')
        if DSL is None:
            # fallback import
            try:
                from DSL import DSL as DSL
            except Exception:
                DSL = None
        if DSL is None:
            raise KeyError('DSL function must be provided in parameters')

        q1Init, q2Init = DSL(parameters['p_L'], y1Init, parameters['T_feed'], parameters['qsb_1'], parameters['qsd_1'], parameters['qsb_2'], parameters['qsd_2'], parameters['bo_1'], parameters['do_1'], parameters['bo_2'], parameters['do_2'], parameters['delUb_1'], parameters['delUd_1'], parameters['delUb_2'], parameters['delUd_2'])
        X0 = np.array([y1Init, q1Init, q2Init, parameters['T_feed'], parameters['T_feed'], parameters['p_L']]) / parameters['refVals']
        parameters['y1init'] = y1Init
        parameters['q1init'] = q1Init
        parameters['q2init'] = q2Init

        max_no_Cycles = 1 if parameters.get('testBT', False) else 200

        parameters['isIsothermal'] = parameters.get('modelType') == 'isothermal'
        parameters['isResin'] = parameters.get('processType') in ('Resin', 'ResinSens')
        parameters['Ab'] = (1.0 - parameters['e_bed']) / parameters['e_bed']
        parameters['darcyK'] = (4.0 / 150.0 / 1.72e-5) * (parameters['e_bed'] / (1.0 - parameters['e_bed'])) ** 2 * parameters['rp'] ** 2
        parameters['tRef_qRef'] = parameters['timeRef'] / parameters['qRef']
        parameters['Ab_rhos_qRef_tRef'] = parameters['Ab'] * parameters['rho_s'] * parameters['qRef'] / parameters['timeRef']
        parameters['inv_tRef'] = 1.0 / parameters['timeRef']
        parameters['PRef_tRef'] = parameters['PRef'] / parameters['timeRef']
        parameters['R_TRef'] = 8.3145 * parameters['TRef']
        parameters['cpg_eV'] = parameters['cp_g'] / (parameters['V_column'] * parameters['e_bed'])
        parameters['two_hin_rin_e'] = 2.0 * parameters['h_in'] / (parameters['r_in'] * parameters['e_bed'])
        parameters['wall_coeff1'] = 2.0 * parameters['h_in'] * parameters['r_in'] / (parameters['r_out'] ** 2 - parameters['r_in'] ** 2)
        parameters['wall_coeff2'] = 2.0 * parameters['h_out'] * parameters['r_out'] / (parameters['r_out'] ** 2 - parameters['r_in'] ** 2)
        parameters['wall_prefactor'] = parameters['timeRef'] / parameters['TwRef'] / (parameters['rho_w'] * parameters['cp_w'])

        temp_check = np.zeros(11)

        cycle = 1

        recovery_percentageValues = []
        purity_percentageValues = []
        productivity_Values = []
        SEC_Values = []

        warnings.filterwarnings('ignore')

        while cycle < max_no_Cycles and temp_check.mean() < 1 and parameters['p_H'] > parameters['p_I'] and parameters['p_I'] > parameters['p_L']:
            parameters['loadingFraction'] = 1
            parameters['P_press'] = _ensure_vector_fn(lambda t: parameters['P_initH'] + (parameters['P_initR'] - parameters['P_initH']) * np.exp(-parameters['lambda'] * np.asarray(t)))
            parameters['dPdt_press'] = _ensure_vector_fn(lambda t: -parameters['lambda'] * (parameters['P_initR'] - parameters['P_initH']) * np.exp(-parameters['lambda'] * np.asarray(t)))

            # solve pressurization
            t_eval_pres = t_press / timeRef
            sol4 = solve_ivp(lambda tt, xx: kbaam_odes_nonisothermal_nd(tt, xx, parameters, 'pres'), (t_eval_pres[0], t_eval_pres[-1]), X0, method=solver_method, t_eval=t_eval_pres, rtol=1e-5, atol=1e-5)
            t4 = sol4.t * timeRef
            X4 = sol4.y.T
            X4[X4 < 0] = 0
            X4[X4[:, 0] > 1, 0] = 1
            if parameters.get('testBT', False) and not parameters.get('testEvac', False):
                parameters['P_initH'] = parameters['p_H'] + 150.0 / 4.0 / parameters['rp'] ** 2 * ((1.0 - parameters['e_bed']) / parameters['e_bed']) ** 2 * 1.72e-5 * parameters['v_in'] * parameters['L']
                y1Init = 1e-6
                if not parameters.get('SSLSTA', False):
                    q1Init, q2Init = DSL(parameters['P_initH'], y1Init, parameters['T_feed'], parameters['qsb_1'], parameters['qsd_1'], parameters['qsb_2'], parameters['qsd_2'], parameters['bo_1'], parameters['do_1'], parameters['bo_2'], parameters['do_2'], parameters['delUb_1'], parameters['delUd_1'], parameters['delUb_2'], parameters['delUd_2'])
                else:
                    sslsta = parameters.get('SSLSTAFunction') or parameters.get('SSLSTA')
                    if not callable(sslsta):
                        raise KeyError('SSLSTA implementation is required when parameters[\'SSLSTA\'] is true.')
                    q1Init, q2Init = sslsta(parameters['P_initH'], y1Init, parameters['T_feed'], parameters)
                X0 = np.array([y1Init, q1Init, q2Init, parameters['T_feed'], parameters['T_feed'], parameters['P_initH']]) / parameters['refVals']
            else:
                X0 = X4[-1, :].copy()
            X4 = X4 * parameters['refVals']
            parameters['q1init'] = X0[1] * parameters['qRef']
            parameters['y1init'] = X0[0]
            parameters['q2init'] = X0[2] * parameters['qRef']

            # loading fraction after pressurization (used to size 2-node blowdown)
            _lf_q1max_p, _ = DSL(X4[-1, 5], parameters['y1_in'], X4[-1, 3], parameters['qsb_1'], parameters['qsd_1'], parameters['qsb_2'], parameters['qsd_2'], parameters['bo_1'], parameters['do_1'], parameters['bo_2'], parameters['do_2'], parameters['delUb_1'], parameters['delUd_1'], parameters['delUb_2'], parameters['delUd_2'])
            _lf_q10_p, _   = DSL(X4[-1, 5], X4[-1, 0], X4[-1, 3], parameters['qsb_1'], parameters['qsd_1'], parameters['qsb_2'], parameters['qsd_2'], parameters['bo_1'], parameters['do_1'], parameters['bo_2'], parameters['do_2'], parameters['delUb_1'], parameters['delUd_1'], parameters['delUb_2'], parameters['delUd_2'])
            _lf_denom_p = _lf_q1max_p - _lf_q10_p
            parameters['loadingFraction'] = (X4[-1, 1] - X4[0, 1]) / (_lf_denom_p if abs(_lf_denom_p) > 1e-12 else 1.0)

            # solve adsorption
            t_eval_ads = t_ads / timeRef
            sol1 = solve_ivp(lambda tt, xx: kbaam_odes_nonisothermal_nd(tt, xx, parameters, 'ads'), (t_eval_ads[0], t_eval_ads[-1]), X0, method=solver_method, t_eval=t_eval_ads, rtol=1e-5, atol=1e-5)
            t1 = sol1.t * timeRef
            X1 = sol1.y.T
            X1[X1 < 0] = 0
            X1[X1[:, 0] > 1, 0] = 1
            X0 = X1[-1, :].copy()
            X1 = X1 * parameters['refVals']

            if parameters.get('processType') in ('Resin', 'ResinSens') and t1[-1] < 0.95 * parameters['t_ads']:
                parameters['t_ads'] = t1[-1]
                t_ads = np.arange(0.0, parameters['t_ads'] + dt, dt)

            if parameters.get('pressureDrop', False):
                parameters['P_initH'] = X1[-1, 5]
                parameters['P_blo'] = _ensure_vector_fn(lambda t: parameters['p_I'] + (parameters['P_initH'] - parameters['p_I']) * np.exp(-parameters['lambda'] * np.asarray(t)))
            # compute inlet/outlet flows for adsorption step
            fin_ads = parameters['volFlowin'] * (2.0 * X1[:, 5] - parameters['p_H']) / (Rg * parameters['T_feed'])
            if parameters.get('pressureDrop', False):
                v_outA = (2.0 / parameters['L']) * parameters['darcyK'] * (X1[:, 5] - parameters['p_H'])
                Fout_ads = parameters['p_H'] * parameters['A_in'] * parameters['e_bed'] / (Rg * X1[:, 3]) * v_outA
                Fout_ads[Fout_ads < 0] = 0
            else:
                grad_q1_ads = np.gradient(X1[:, 1], dt)
                grad_q2_ads = np.gradient(X1[:, 2], dt)
                Fout_ads = parameters['F_in'] - (1 - parameters['e_bed']) * parameters['V_column'] * parameters['rho_s'] * (grad_q1_ads + grad_q2_ads)
                Fout_ads[Fout_ads < 0] = 0
            F_1_out_ads = Fout_ads * X1[:, 0]

            # integrate molar outputs
            mol_1_out_ads = np.trapz(F_1_out_ads, t1)
            moltot_out_ads = np.trapz(Fout_ads, t1)
            parameters['y1_LPP'] = mol_1_out_ads / moltot_out_ads if moltot_out_ads != 0 else 0.0
            parameters['y1init'] = X1[-1, 0]
            parameters['q1init'] = X1[-1, 1]
            parameters['q2init'] = X1[-1, 2]

            # loading fraction after adsorption (for 2-node blowdown sizing)
            _lf_q1max_a, _ = DSL(X1[-1, 5], parameters['y1_in'], X1[-1, 3], parameters['qsb_1'], parameters['qsd_1'], parameters['qsb_2'], parameters['qsd_2'], parameters['bo_1'], parameters['do_1'], parameters['bo_2'], parameters['do_2'], parameters['delUb_1'], parameters['delUd_1'], parameters['delUb_2'], parameters['delUd_2'])
            _lf_q10_a, _   = DSL(X4[-1, 5], X4[-1, 0], X4[-1, 3], parameters['qsb_1'], parameters['qsd_1'], parameters['qsb_2'], parameters['qsd_2'], parameters['bo_1'], parameters['do_1'], parameters['bo_2'], parameters['do_2'], parameters['delUb_1'], parameters['delUd_1'], parameters['delUb_2'], parameters['delUd_2'])
            _lf_denom_a = _lf_q1max_a - X1[0, 1]
            parameters['loadingFraction'] = (X1[-1, 1] - X1[0, 1]) / (_lf_denom_a if abs(_lf_denom_a) > 1e-12 else 1.0)
            f_w = max(0.05, min(0.95, parameters['loadingFraction']))

            # blowdown
            t_eval_blo = t_blo / timeRef
            if parameters.get('pressureDrop', False) and f_w < 0.95 and not parameters.get('amine', False):
                # 2-node blowdown
                refVals_blo = np.array([
                    1.0, parameters['qRef'], parameters['qRef'],   # node 1: y1_1, q1_1, q2_1
                    1.0, parameters['qRef'], parameters['qRef'],   # node 2: y1_2, q1_2, q2_2
                    parameters['TRef'], parameters['TwRef'],       # shared T, Tw
                    parameters['PRef'], parameters['PRef'],        # P1, P2
                ])
                P_ads_end = parameters['p_H']
                T_ads_end = X1[-1, 3]
                q1_n1, q2_n1 = DSL(P_ads_end, parameters['y1_in'], T_ads_end, parameters['qsb_1'], parameters['qsd_1'], parameters['qsb_2'], parameters['qsd_2'], parameters['bo_1'], parameters['do_1'], parameters['bo_2'], parameters['do_2'], parameters['delUb_1'], parameters['delUd_1'], parameters['delUb_2'], parameters['delUd_2'])
                y1_n2 = X1[-1, 0]
                q1_n2 = X1[ 0, 1]
                q2_n2 = X1[ 0, 2]
                X0_blo = np.array([
                    parameters['y1_in'], q1_n1, q2_n1,
                    y1_n2, q1_n2, q2_n2,
                    T_ads_end, X1[-1, 4],
                    P_ads_end, P_ads_end,
                ]) / refVals_blo
                sol2_10 = solve_ivp(
                    lambda tt, xx: kbaam_odes_nonisothermal_nd_dp_blo2node(tt, xx, parameters, 'blo'),
                    (t_eval_blo[0], t_eval_blo[-1]), X0_blo,
                    method=solver_method, t_eval=t_eval_blo, rtol=1e-5, atol=1e-5,
                )
                t2 = sol2_10.t * timeRef
                X2_10 = sol2_10.y.T
                X2_10[X2_10 < 0] = 0
                X2_10[X2_10[:, 0] > 1, 0] = 1
                X2_10[X2_10[:, 3] > 1, 3] = 1
                X2_10dim = X2_10 * refVals_blo
                f_blo = max(0.05, min(0.9, parameters['loadingFraction']))
                y1_blo_end = f_blo * X2_10dim[-1, 0] + (1 - f_blo) * X2_10dim[-1, 3]
                q1_blo_end = f_blo * X2_10dim[-1, 1] + (1 - f_blo) * X2_10dim[-1, 4]
                q2_blo_end = f_blo * X2_10dim[-1, 2] + (1 - f_blo) * X2_10dim[-1, 5]
                T_blo_end  = X2_10dim[-1, 6]
                Tw_blo_end = X2_10dim[-1, 7]
                P_blo_end  = f_blo * X2_10dim[-1, 8] + (1 - f_blo) * X2_10dim[-1, 9]
                X0 = np.array([y1_blo_end, q1_blo_end, q2_blo_end, T_blo_end, Tw_blo_end, P_blo_end]) / parameters['refVals']
                P_avg_blo = f_blo * X2_10dim[:, 8] + (1 - f_blo) * X2_10dim[:, 9]
                X2 = np.column_stack([
                    f_blo * X2_10dim[:, 0] + (1 - f_blo) * X2_10dim[:, 3],
                    f_blo * X2_10dim[:, 1] + (1 - f_blo) * X2_10dim[:, 4],
                    f_blo * X2_10dim[:, 2] + (1 - f_blo) * X2_10dim[:, 5],
                    X2_10dim[:, 6], X2_10dim[:, 7], P_avg_blo,
                ])
                y1_bd_out = X2_10dim[:, 3]  # node-2 (product end) composition
            else:
                sol2 = solve_ivp(lambda tt, xx: kbaam_odes_nonisothermal_nd(tt, xx, parameters, 'blo'), (t_eval_blo[0], t_eval_blo[-1]), X0, method=solver_method, t_eval=t_eval_blo, rtol=1e-5, atol=1e-5)
                t2 = sol2.t * timeRef
                X2 = sol2.y.T
                X2[X2 < 0] = 0
                X2[X2[:, 0] > 1, 0] = 1
                X0 = X2[-1, :].copy()
                X2 = X2 * parameters['refVals']
                f_blo = 1.0
                y1_bd_out = X2[:, 0]
                X2_10dim = None

            # evacuation
            t_eval_evac = t_evac / timeRef
            sol3 = solve_ivp(lambda tt, xx: kbaam_odes_nonisothermal_nd(tt, xx, parameters, 'evac'), (t_eval_evac[0], t_eval_evac[-1]), X0, method=solver_method, t_eval=t_eval_evac, rtol=1e-5, atol=1e-5)
            t3 = sol3.t * timeRef
            X3 = sol3.y.T
            X3[X3 < 0] = 0
            X3[X3[:, 0] > 1, 0] = 1
            X0 = X3[-1, :].copy()
            X3 = X3 * parameters['refVals']

            if parameters.get('processType') in ('Resin', 'ResinSens') and t3[-1] < 0.95 * parameters['t_evac']:
                parameters['t_evac'] = t3[-1]
                t_evac = np.arange(0.0, parameters['t_evac'] + dt, dt)
                parameters['P_initR'] = parameters['P_evac'](parameters['t_evac'])
                cycle = 1

            # compute molar balances
            n_1_ads = (X1[-1, 0] * X1[-1, 5] * parameters['V_column'] * parameters['e_bed'] / (Rg * X1[-1, 3])) + X1[-1, 1] * parameters['V_column'] * (1 - parameters['e_bed']) * parameters['rho_s']
            n_1_bd = (X3[ 0, 0] * X3[ 0, 5] * parameters['V_column'] * parameters['e_bed'] / (Rg * X3[ 0, 3])) + X3[ 0, 1] * parameters['V_column'] * (1 - parameters['e_bed']) * parameters['rho_s']
            n_1_evac = (X3[-1, 0] * X3[-1, 5] * parameters['V_column'] * parameters['e_bed'] / (Rg * X3[-1, 3])) + X3[-1, 1] * parameters['V_column'] * (1 - parameters['e_bed']) * parameters['rho_s']
            n_2_bd = ((1 - X3[ 0, 0]) * X3[ 0, 5] * parameters['V_column'] * parameters['e_bed'] / (Rg * X3[ 0, 3])) + X3[ 0, 2] * parameters['V_column'] * (1 - parameters['e_bed']) * parameters['rho_s']
            n_2_evac = ((1 - X3[-1, 0]) * X3[-1, 5] * parameters['V_column'] * parameters['e_bed'] / (Rg * X3[-1, 3])) + X3[-1, 2] * parameters['V_column'] * (1 - parameters['e_bed']) * parameters['rho_s']
            n_1_pres = (X4[-1, 0] * X4[-1, 5] * parameters['V_column'] * parameters['e_bed'] / (Rg * X4[-1, 3])) + X4[-1, 1] * parameters['V_column'] * (1 - parameters['e_bed']) * parameters['rho_s']
            n_1_presInit = (X4[0, 0] * X4[0, 5] * parameters['V_column'] * parameters['e_bed'] / (Rg * X4[0, 3])) + X4[0, 1] * parameters['V_column'] * (1 - parameters['e_bed']) * parameters['rho_s']

            if parameters.get('pressType') == 'LPP':
                mole_lp_recycle = n_1_pres - n_1_presInit
                denom = n_1_ads - n_1_evac + mol_1_out_ads - max(0.0, mole_lp_recycle)
                recovery_percentage = 100.0 * (n_1_bd - n_1_evac) / (denom if denom != 0 else 1.0)
            else:
                denom = n_1_ads - n_1_evac + mol_1_out_ads
                recovery_percentage = 100.0 * (n_1_bd - n_1_evac) / (denom if denom != 0 else 1.0)
            purity_denom = n_1_bd - n_1_evac + max(0.0, n_2_bd - n_2_evac)
            purity_percentage = 100.0 * (n_1_bd - n_1_evac) / (purity_denom if purity_denom != 0 else 1.0)

            cycle_time = (parameters['t_ads'] + parameters['t_blo'] + parameters['t_evac'] + parameters['t_press'])
            productivity = (n_1_bd - n_1_evac) / (parameters['V_column'] * cycle_time if parameters['V_column'] * cycle_time != 0 else 1.0)

            # Energy calculation — step flowrates
            fin_ads = parameters['volFlowin'] * (2.0 * X1[:, 5] - parameters['p_H']) / (Rg * parameters['T_feed'])

            if parameters.get('pressureDrop', False) and f_w < 0.95 and not parameters.get('amine', False):
                # 2-node blowdown: outlet flow from node-2 (product end)
                L2_blo_e = (1 - f_blo) * parameters['L']
                v_outB = (2.0 / L2_blo_e) * parameters['darcyK'] * (X2_10dim[:, 9] - parameters['P_blo'](t2))
                Fout_bd = X2_10dim[:, 9] * parameters['A_in'] * parameters['e_bed'] / (Rg * X2_10dim[:, 6]) * v_outB
                Fout_bd[Fout_bd < 0] = 0
                v_outE = (2.0 / parameters['L']) * parameters['darcyK'] * (X3[:, 5] - parameters['P_evac'](t3))
                Fout_evac = X3[:, 5] * parameters['A_in'] * parameters['e_bed'] / (Rg * X3[:, 3]) * v_outE
                Fout_evac[Fout_evac < 0] = 0
                v_outP = -(2.0 / parameters['L']) * parameters['darcyK'] * (X4[:, 5] - parameters['P_press'](t4))
                Fin_pres = parameters['P_press'](t4) * parameters['A_in'] * parameters['e_bed'] / (Rg * parameters['T_feed']) * v_outP
            elif parameters.get('pressureDrop', False) and not parameters.get('amine', False):
                # single-node with pressure drop
                v_outB = (2.0 / parameters['L']) * parameters['darcyK'] * (X2[:, 5] - parameters['P_blo'](t2))
                Fout_bd = X2[:, 5] * parameters['A_in'] * parameters['e_bed'] / (Rg * X2[:, 3]) * v_outB
                Fout_bd[Fout_bd < 0] = 0
                v_outE = (2.0 / parameters['L']) * parameters['darcyK'] * (X3[:, 5] - parameters['P_evac'](t3))
                Fout_evac = X3[:, 5] * parameters['A_in'] * parameters['e_bed'] / (Rg * X3[:, 3]) * v_outE
                Fout_evac[Fout_evac < 0] = 0
                v_outP = -(2.0 / parameters['L']) * parameters['darcyK'] * (X4[:, 5] - parameters['P_press'](t4))
                Fin_pres = parameters['P_press'](t4) * parameters['A_in'] * parameters['e_bed'] / (Rg * parameters['T_feed']) * v_outP
            else:
                # no pressure drop: gradient-based material balance
                grad_X2_q1 = np.gradient(X2[:, 1], dt)
                grad_X2_q2 = np.gradient(X2[:, 2], dt)
                grad_X2_p  = np.gradient(X2[:, 5], dt)
                grad_X2_t  = np.gradient(X2[:, 3], dt)
                grad_X3_q1 = np.gradient(X3[:, 1], dt)
                grad_X3_q2 = np.gradient(X3[:, 2], dt)
                grad_X3_p  = np.gradient(X3[:, 5], dt)
                grad_X3_t  = np.gradient(X3[:, 3], dt)
                grad_X4_q1 = np.gradient(X4[:, 1], dt)
                grad_X4_q2 = np.gradient(X4[:, 2], dt)
                grad_X4_p  = np.gradient(X4[:, 5], dt)
                grad_X4_t  = np.gradient(X4[:, 3], dt)
                Fout_bd = 0 - (1 - parameters['e_bed']) * parameters['V_column'] * parameters['rho_s'] * (grad_X2_q1 + grad_X2_q2) - (parameters['e_bed'] / (Rg * X2[:, 3])) * grad_X2_p * parameters['V_column'] + (parameters['e_bed'] * X2[:, 5] / (Rg * X2[:, 3] ** 2)) * grad_X2_t * parameters['V_column']
                Fout_bd[Fout_bd < 0] = 0
                Fout_evac = 0 - (1 - parameters['e_bed']) * parameters['V_column'] * parameters['rho_s'] * (grad_X3_q1 + grad_X3_q2) - 0.0 * (parameters['e_bed'] / (Rg * X3[:, 3])) * grad_X3_p * parameters['V_column'] + (parameters['e_bed'] * X3[:, 5] / (Rg * X3[:, 3] ** 2)) * grad_X3_t * parameters['V_column']
                Fout_evac[Fout_evac < 0] = 0
                Fin_pres = (1 - parameters['e_bed']) * parameters['V_column'] * parameters['rho_s'] * (grad_X4_q1 + grad_X4_q2) + (parameters['e_bed'] / (Rg * X4[:, 3])) * grad_X4_p * parameters['V_column'] - (parameters['e_bed'] * X4[:, 5] / (Rg * X4[:, 3] ** 2)) * grad_X4_t * parameters['V_column']

            p_out_bd = parameters['P_blo'](t2)
            p_out_evac = parameters['P_evac'](t3)
            eta_bd = 0.8 * (19.55 * p_out_bd * 1e-5 / (1 + 19.55 * p_out_bd * 1e-5))
            eta_evac = 0.8 * (19.55 * p_out_evac * 1e-5 / (1 + 19.55 * p_out_evac * 1e-5))
            eta_press = 0.72
            eta_ads = 0.72
            p_atm = 101325.0

            if parameters.get('heating', False):
                heat_flag = X3[:, 3] < parameters['Theat']
                qheat = np.trapz(t3, parameters['heatPowerDensity'] * (parameters['Theat'] - X3[:, 3]) / (parameters['Theat'] - parameters['T_feed']) * parameters['r_out'] * 2.0 * parameters['L'] * heat_flag)
                ec_heat = qheat
            else:
                ec_heat = 0.0

            EC_BD = np.trapz((1.0 / eta_bd) * Fout_bd * Rg * X2[:, 3] * (1.4 / 0.4) * ((p_atm / np.minimum(p_atm, p_out_bd)) ** (0.4 / 1.4) - 1.0), t2)
            EC_EVAC = np.trapz((1.0 / eta_evac) * Fout_evac * Rg * X3[:, 3] * (1.4 / 0.4) * ((p_atm / np.minimum(p_atm, p_out_evac)) ** (0.4 / 1.4) - 1.0), t3)
            EC_PRES = np.trapz((1.0 / eta_press) * Fin_pres * Rg * parameters['T_feed'] * (1.4 / 0.4) * ((np.maximum(p_atm, parameters['P_press'](t4)) / p_atm) ** (0.4 / 1.4) - 1.0), t4)
            EC_ADS = np.trapz((1.0 / eta_ads) * fin_ads * Rg * parameters['T_feed'] * (1.4 / 0.4) * ((np.maximum(p_atm, 2.0 * X1[:, 5] - parameters['p_H']) / p_atm) ** (0.4 / 1.4) - 1.0), t1)

            sec_denom = (n_1_bd - n_1_evac) * 0.04401
            SEC = (EC_PRES + EC_ADS + EC_BD + EC_EVAC + ec_heat) / (sec_denom if sec_denom != 0 else 1.0) / 3600.0

            recovery_percentageValues.append(recovery_percentage)
            purity_percentageValues.append(purity_percentage)
            productivity_Values.append(productivity)
            SEC_Values.append(SEC)

            process_indicators = np.vstack([purity_percentageValues, recovery_percentageValues, productivity_Values, SEC_Values])

            # Match MATLAB CSS convergence: all 4 KPIs within 0.005%
            # for 11 consecutive cycles, referenced against cycle-5.
            if cycle > 11:
                for i in range(4):
                    for k in range(11):
                        prev = process_indicators[i, cycle - 6]
                        curr = process_indicators[i, cycle - k - 1]
                        if prev != 0 and abs(100 * (curr - prev) / prev) <= 0.005:
                            temp_check[k] = 1
                        else:
                            temp_check[k] = 0

            cycle += 1

        # Compute KPIs and penalties
        purity_percentage = purity_percentageValues[-1] if len(purity_percentageValues) else 0.0
        recovery_percentage = recovery_percentageValues[-1] if len(recovery_percentageValues) else 0.0
        productivity = productivity_Values[-1] if len(productivity_Values) else 0.0
        SEC = SEC_Values[-1] if len(SEC_Values) else 0.0

        phi_pen = [0.0, 0.0]
        phi_pen[0] = 0.80 * ((1.0 * max(0.0, (95 - purity_percentage))) ** 2) + (max(0.0, (90 - recovery_percentage)) ** 2)
        phi_pen[1] = 0.30 * (((1.0 * max(0.0, (95 - purity_percentage))) ** 2) + (max(0.0, (90 - recovery_percentage)) ** 2))
        
        # pdb.set_trace()
        # Only drop into the debugger when explicitly requested via parameters['debug']
        if parameters.get('debug', False):
            pdb.set_trace()

        if parameters.get('outputType') == 'plot' and 'process_indicators' in locals():
            KPIs = process_indicators.T
        elif parameters.get('OptType') == 'Const':
            KPIs = [(-productivity + phi_pen[0]), 10 * (SEC * 2.77778e-7 + phi_pen[1])]
        else:
            KPIs = [-recovery_percentage, -purity_percentage]

    except Exception:
        # On error return large KPIs like the original MATLAB implementation
        if raise_on_error:
            raise
        KPIs = [1e5, 1e5]
        SEC = 1e5
        
    
    #Python 3: 
    # print(datetime.now() - startTime)


    # Optionally plot (basic reproduction of MATLAB plotting)
    if parameters.get('outputType') == 'plot' and all(name in locals() for name in ('t1', 't2', 't3', 't4', 'X1', 'X2', 'X3', 'X4')):
        import matplotlib.pyplot as plt

        t_cycle = np.concatenate([t1, t2 + t1[-1], t3 + t1[-1] + t2[-1], t4 + t1[-1] + t2[-1] + t3[-1]])
        X_cycle = np.vstack([X1, X2, X3, X4])
        P_cycle = X_cycle[:, 5]

        fig, axs = plt.subplots(6, 1, figsize=(7, 10), sharex=True)
        axs[0].plot(t_cycle, P_cycle / 1e5, '-', color='#0B0', linewidth=2)
        axs[0].set_ylabel('P [bar]')
        axs[0].set_title('Column Pressure (Cyclic Steady State)')
        
        # q* and q1
        q1_starvals, q2_starvals = DSL(P_cycle, X_cycle[:, 0], X_cycle[:, 3], parameters['qsb_1'], parameters['qsd_1'], parameters['qsb_2'], parameters['qsd_2'], parameters['bo_1'], parameters['do_1'], parameters['bo_2'], parameters['do_2'], parameters['delUb_1'], parameters['delUd_1'], parameters['delUb_2'], parameters['delUd_2'])
        q1_starvals_ads, q2_starvals_ads = DSL(parameters['P_ads'](t_ads * timeRef), parameters['y1_in'], X1[:, 3], parameters['qsb_1'], parameters['qsd_1'], parameters['qsb_2'], parameters['qsd_2'], parameters['bo_1'], parameters['do_1'], parameters['bo_2'], parameters['do_2'], parameters['delUb_1'], parameters['delUd_1'], parameters['delUb_2'], parameters['delUd_2'])
        q1_starvals[:len(t_ads)] = q1_starvals_ads
        q2_starvals[:len(t_ads)] = q2_starvals_ads

        axs[1].plot(t_cycle, q1_starvals, 'k--', linewidth=2)
        axs[1].plot(t_cycle, X_cycle[:, 1], 'b-', linewidth=2)
        axs[1].set_ylabel('q_1 [mol/kg]')
        axs[1].set_title('Adsorbed amount of Component 1 (Cyclic Steady State)')

        axs[2].plot(t_cycle, q2_starvals, 'k--', linewidth=2)
        axs[2].plot(t_cycle, X_cycle[:, 2], 'r--', linewidth=2)
        axs[2].set_ylabel('q_1 [mol/kg]')
        axs[2].set_title('Adsorbed amount of Component 2 (Cyclic Steady State)')

        axs[3].plot(t_cycle, X_cycle[:, 3], 'm-', linewidth=2)
        axs[3].set_ylabel('T [K]')
        axs[3].set_title('Bed/Gas temperature (Cyclic Steady State)')

        axs[4].plot(t_cycle, X_cycle[:, 4], 'g-', linewidth=2)
        axs[4].set_ylabel('Tw [K]')
        axs[4].set_title('Column wall temperature (Cyclic Steady State)')

        axs[5].plot(t_cycle, X_cycle[:, 0], 'k-', linewidth=2)
        axs[5].set_ylabel('y_1')
        axs[5].set_ylim([-0.05, 1.05])
        axs[5].set_title('Mole fraction of Component 1 (Cyclic Steady State)')
        plt.xlabel('time [s]')
        plt.tight_layout()
        plot_path = parameters.get('plotPath')
        if plot_path:
            fig.savefig(plot_path, dpi=160, bbox_inches='tight')
        plt.show()

    # parameters.get('outputType') = 'opt'
    # If not plotting, save decision variables / results to rawData/<fileName>.txt
    if parameters.get('outputType') != 'plot':
        try:
            raw_dir = 'rawData'
            os.makedirs(raw_dir, exist_ok=True)
            file_path = os.path.join(raw_dir, f"{parameters.get('fileName','run')}.txt")

            # Determine theta to save: prefer the caller-provided thetaIn, else
            # try to reconstruct relevant elements for logging.
            if thetaIn is not None:
                theta_to_save = np.asarray(thetaIn)
            else:
                # Construct a fallback theta for non-Adsorbent cases
                pt = parameters.get('processType', 'PVSA')
                if pt == 'PVSA':
                    theta_to_save = np.array([parameters['v_in'], parameters['p_I'], parameters['t_ads'], parameters['t_blo'], parameters['t_evac'], parameters['p_H']])
                elif pt == 'VSA':
                    theta_to_save = np.array([parameters['v_in'], parameters['p_I'], parameters['t_ads'], parameters['t_blo'], parameters['t_evac']])
                else:
                    # For Adsorbent* variants try to build a 7-entry theta if possible
                    theta_to_save = None

            with open(file_path, 'a+') as fh:
                if parameters.get('processType','').startswith('Adsorbent') and theta_to_save is not None and theta_to_save.size >= 7:
                    # MATLAB writes p_H, p_I, p_L, F_in, t_ads, t_blo, t_evac, purity,recovery,productivity,SEC, theta(1:7)
                    vals = [parameters.get('p_H', np.nan), parameters.get('p_I', np.nan), parameters.get('p_L', np.nan), parameters.get('F_in', np.nan), parameters.get('t_ads', np.nan), parameters.get('t_blo', np.nan), parameters.get('t_evac', np.nan), purity_percentage, recovery_percentage, productivity, SEC]
                    # ensure we have at least 7 theta entries
                    theta_head = theta_to_save[:7]
                    line_vals = np.concatenate([np.asarray(vals, dtype=float), np.asarray(theta_head, dtype=float)])
                    fmt = ' '.join(['%12.9f'] * len(line_vals)) + '\n'
                    fh.write(fmt % tuple(line_vals.tolist()))
                else:
                    # Regular case: p_H, p_I, p_L, F_in, t_ads, t_blo, t_evac, purity,recovery,productivity,SEC, V_column, v_in
                    vals = [parameters.get('p_H', np.nan), parameters.get('p_I', np.nan), parameters.get('p_L', np.nan), parameters.get('F_in', np.nan), parameters.get('t_ads', np.nan), parameters.get('t_blo', np.nan), parameters.get('t_evac', np.nan), purity_percentage, recovery_percentage, productivity, SEC, parameters.get('V_column', np.nan), parameters.get('v_in', np.nan)]
                    fmt = ' '.join(['%12.9f'] * len(vals)) + '\n'
                    fh.write(fmt % tuple([float(x) for x in vals]))
        except Exception:
            # Do not let logging break the main result
            pass

    return KPIs


if __name__ == '__main__':
    # quick run using default parameters
    params = create_parameters()
    print('Running kBAAM outputs smoke run (debug mode)...')
    try:
        kpis = kbaam_outputs_nonisothermal(params, raise_on_error=True)
        print('KPIs:', kpis)
    except Exception as exc:
        import traceback

        print('kBAAM outputs raised an exception:')
        traceback.print_exc()
