"""kBAAM_Outputs_nonIsothermal.py

Imperial College London, Multiphase Systems Laboratory
Year: 2025
Author: Hassan Azzan (HA)

Purpose:
    Solve the non-isothermal kinetic batch adsorber analogue model (k-BAAM)
    and compute key performance indicators (KPIs) such as heavy product
    recovery and purity. This is a Python port of the MATLAB
    ``kBAAM_Outputs_nonIsothermal.m`` driver.

Function:
    ``kbaam_outputs_nonisothermal(parameters, thetaIn=None, ...) -> KPIs``

Inputs:
    - parameters: a dict-like parameters structure (see ``createParameters.py``).
    - thetaIn: optional decision vector to override process parameters

Outputs:
    - KPIs: list/array, typically [objective1, objective2] matching MATLAB

Dependencies:
    - kBAAM_ODEs_nonIsothermal_ND
    - DSL
    - LDFCoefficient

Notes:
    - The loader and default parameter handling attempt to mimic the MATLAB
      behaviour (pressure schedules, non-dimensionalisation, defaults,
      fallback DSL/LDF implementations).
"""
from __future__ import annotations
import numpy as np
from typing import Any, Dict, Sequence, Optional
import warnings

from scipy.integrate import solve_ivp

from kBAAM_ODEs_nonIsothermal_ND import kbaam_odes_nonisothermal_nd
from createParameters import create_parameters
import pdb
import matplotlib_inline;
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
    matplotlib_inline.backend_inline.set_matplotlib_formats('png', 'jpeg');
        
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
                if parameters.get('outputType') == 'opt':
                    parameters['p_I'] = 10 ** parameters['p_I']
            elif pt == 'PVSA':
                parameters['v_in'] = thetavec[0]
                parameters['p_I'] = thetavec[1]
                parameters['t_ads'] = thetavec[2]
                parameters['t_blo'] = thetavec[3]
                parameters['t_evac'] = thetavec[4]
                parameters['p_H'] = thetavec[5]
                if parameters.get('outputType') == 'opt':
                    parameters['p_I'] = 10 ** parameters['p_I']
                    parameters['p_H'] = 10 ** parameters['p_H']
            elif pt.startswith('Adsorbent'):
                # handle Adsorbent* variants similarly to MATLAB script (partial)
                parameters['v_in'] = 0.4
                parameters['p_I'] = thetavec[0]
                parameters['t_ads'] = thetavec[1]
                parameters['t_blo'] = thetavec[2]
                # other fields set below depending on specific variant
                # For brevity we mimic MATLAB assignments where straightforward

        Rg = 8.3145

        dt = 0.1
        t_ads = np.arange(0.0, parameters['t_ads'] + dt, dt)
        t_blo = np.arange(0.0, parameters['t_blo'] + dt, dt)
        t_evac = np.arange(0.0, parameters['t_evac'] + dt, dt)
        t_press = np.arange(0.0, parameters['t_press'] + dt, dt)

        qRef = parameters['qsb_1'] + parameters['qsd_1']
        parameters['refVals'] = np.array([1.0, qRef, qRef, parameters['T_feed'], parameters['T_feed']])

        # geometry
        parameters['Lbyr'] = 6
        parameters['r_in'] = (parameters['V_column'] / (parameters['Lbyr'] * np.pi)) ** (1.0 / 3.0)
        parameters['r_out'] = parameters.get('r_in', parameters['r_in']) + parameters.get('t_wall', 0.003)
        parameters['L'] = parameters['Lbyr'] * parameters['r_in']
        A_in = parameters['r_in'] ** 2 * np.pi
        volFlowin = parameters['v_in'] * A_in
        parameters['F_in'] = volFlowin * parameters['p_H'] / (Rg * parameters['T_feed'])
        volFluxRef = parameters['F_in'] / parameters['V_column']
        timeRef = parameters['p_H'] / (Rg * parameters['T_feed'] * volFluxRef)
        # pdb.set_trace()
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
        P_ads_def = lambda t: parameters['p_H']
        P_blo_def = lambda t: parameters['p_I'] + (parameters['p_H'] - parameters['p_I']) * np.exp(-parameters['lambda'] * t)
        P_evac_def = lambda t: parameters['p_L'] + (parameters.get('P_initL', parameters.get('p_I')) - parameters['p_L']) * np.exp(-parameters['lambda'] * t)
        P_press_def = lambda t: parameters['p_H'] + (parameters.get('P_initR', parameters.get('p_H')) - parameters['p_H']) * np.exp(-3 * t)

        dPdt_blo_def = lambda t: -parameters['lambda'] * (parameters['p_H'] - parameters['p_I']) * np.exp(-parameters['lambda'] * t)
        dPdt_evac_def = lambda t: -parameters['lambda'] * (parameters['p_I'] - parameters['p_L']) * np.exp(-parameters['lambda'] * t)
        dPdt_press_def = lambda t: -parameters['lambda'] * (parameters['p_L'] - parameters['p_H']) * np.exp(-parameters['lambda'] * t)

        # attach defaults or use provided functions/values, ensuring vector behaviour
        parameters['P_ads'] = _ensure_vector_fn(parameters.get('P_ads', P_ads_def))
        parameters['P_blo'] = _ensure_vector_fn(parameters.get('P_blo', P_blo_def))
        # P_initL depends on P_blo; call with scalar t_blo
        parameters['P_initL'] = parameters['P_blo'](parameters['t_blo'])
        parameters['P_evac'] = _ensure_vector_fn(parameters.get('P_evac', P_evac_def))
        # P_initR depends on P_evac; call with scalar t_evac
        parameters['P_initR'] = parameters['P_evac'](parameters['t_evac'])
        parameters['P_press'] = _ensure_vector_fn(parameters.get('P_press', P_press_def))

        parameters['dPdt_blo'] = _ensure_vector_fn(parameters.get('dPdt_blo', dPdt_blo_def))
        parameters['dPdt_evac'] = _ensure_vector_fn(parameters.get('dPdt_evac', dPdt_evac_def))
        parameters['dPdt_press'] = _ensure_vector_fn(parameters.get('dPdt_press', dPdt_press_def))
        
        
        # initial conditions
        y1Init = 0.01
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

        q1Init, q2Init = DSL(parameters['p_H'], y1Init, parameters['T_feed'], parameters['qsb_1'], parameters['qsd_1'], parameters['qsb_2'], parameters['qsd_2'], parameters['bo_1'], parameters['do_1'], parameters['bo_2'], parameters['do_2'], parameters['delUb_1'], parameters['delUd_1'], parameters['delUb_2'], parameters['delUd_2'])
        X0 = np.array([y1Init, q1Init, q2Init, parameters['T_feed'], parameters['T_feed']]) / parameters['refVals']

        max_no_Cycles = 100

        if parameters.get('heating'):
            parameters['Theat'] = 493
        else:
            parameters['Theat'] = 0

        temp_check = np.zeros(5)

        cycle = 1

        recovery_percentageValues = []
        purity_percentageValues = []
        productivity_Values = []
        SEC_Values = []

        warnings.filterwarnings('ignore')

        while cycle < max_no_Cycles and temp_check.mean() < 1:
            # solve adsorption
            t_eval_ads = t_ads / timeRef
            sol1 = solve_ivp(lambda tt, xx: kbaam_odes_nonisothermal_nd(tt, xx, parameters, 'ads'), (t_eval_ads[0], t_eval_ads[-1]), X0, method=solver_method, t_eval=t_eval_ads, rtol=1e-5, atol=1e-5)
            t1 = sol1.t * timeRef
            X1 = sol1.y.T
            X1[X1 < 0] = 0
            X1[X1[:, 0] > 1, 0] = 1
            X0 = X1[-1, :].copy()
            X1 = X1 * parameters['refVals']
            # compute outlet flows for adsorption
            grad_q1 = np.gradient(X1[:, 1], dt * timeRef)
            grad_q2 = np.gradient(X1[:, 2], dt * timeRef)
            Fout_ads = parameters['F_in'] - (1 - parameters['e_bed']) * parameters['V_column'] * parameters['rho_s'] * (grad_q1 + grad_q2)
            Fout_ads[Fout_ads < 0] = 0
            F_1_out_ads = Fout_ads * X1[:, 0]

            # integrate molar outputs
            mol_1_out_ads = np.trapz(F_1_out_ads, t1 * timeRef)
            moltot_out_ads = np.trapz(Fout_ads, t1 * timeRef)
            parameters['y1_LPP'] = mol_1_out_ads / moltot_out_ads if moltot_out_ads != 0 else 0.0

            # blowdown
            t_eval_blo = t_blo / timeRef
            sol2 = solve_ivp(lambda tt, xx: kbaam_odes_nonisothermal_nd(tt, xx, parameters, 'blo'), (t_eval_blo[0], t_eval_blo[-1]), X0, method=solver_method, t_eval=t_eval_blo, rtol=1e-5, atol=1e-5)
            t2 = sol2.t * timeRef
            X2 = sol2.y.T
            X2[X2 < 0] = 0
            X2[X2[:, 0] > 1, 0] = 1
            X0 = X2[-1, :].copy()
            X2 = X2 * parameters['refVals']

            # evacuation
            t_eval_evac = t_evac / timeRef
            sol3 = solve_ivp(lambda tt, xx: kbaam_odes_nonisothermal_nd(tt, xx, parameters, 'evac'), (t_eval_evac[0], t_eval_evac[-1]), X0, method=solver_method, t_eval=t_eval_evac, rtol=1e-5, atol=1e-5)
            t3 = sol3.t * timeRef
            X3 = sol3.y.T
            X3[X3 < 0] = 0
            X3[X3[:, 0] > 1, 0] = 1
            X0 = X3[-1, :].copy()
            X3 = X3 * parameters['refVals']

            if t3[-1] < 0.95 * parameters['t_evac']:
                parameters['t_evac'] = t3[-1]
                t_evac = np.arange(0.0, parameters['t_evac'] + dt, dt)
                parameters['P_initR'] = parameters['P_evac'](parameters['t_evac'])
                cycle = 1

            # pressurization
            t_eval_pres = t_press / timeRef
            sol4 = solve_ivp(lambda tt, xx: kbaam_odes_nonisothermal_nd(tt, xx, parameters, 'pres'), (t_eval_pres[0], t_eval_pres[-1]), X0, method=solver_method, t_eval=t_eval_pres, rtol=1e-5, atol=1e-5)
            t4 = sol4.t * timeRef
            X4 = sol4.y.T
            X4[X4 < 0] = 0
            X4[X4[:, 0] > 1, 0] = 1
            X0 = X4[-1, :].copy()
            X4 = X4 * parameters['refVals']

            # compute molar balances
            n_1_bd = (X2[-1, 0] * parameters['P_blo'](t2[-1]) * parameters['V_column'] * parameters['e_bed'] / (Rg * X2[-1, 3])) + X2[-1, 1] * parameters['V_column'] * (1 - parameters['e_bed']) * parameters['rho_s']
            n_1_evac = (X3[-1, 0] * parameters['P_evac'](t3[-1]) * parameters['V_column'] * parameters['e_bed'] / (Rg * X3[-1, 3])) + X3[-1, 1] * parameters['V_column'] * (1 - parameters['e_bed']) * parameters['rho_s']
            n_2_bd = ((1 - X2[-1, 0]) * parameters['P_blo'](t2[-1]) * parameters['V_column'] * parameters['e_bed'] / (Rg * X2[-1, 3])) + X2[-1, 2] * parameters['V_column'] * (1 - parameters['e_bed']) * parameters['rho_s']
            n_2_evac = ((1 - X3[-1, 0]) * parameters['P_evac'](t3[-1]) * parameters['V_column'] * parameters['e_bed'] / (Rg * X3[-1, 3])) + X3[-1, 2] * parameters['V_column'] * (1 - parameters['e_bed']) * parameters['rho_s']
            n_1_pres = (X4[-1, 0] * parameters['P_press'](t4[-1]) * parameters['V_column'] * parameters['e_bed'] / (Rg * X4[-1, 3])) + X4[-1, 1] * parameters['V_column'] * (1 - parameters['e_bed']) * parameters['rho_s']
            n_1_adsin = parameters['F_in'] * parameters['y1_in'] * parameters['t_ads']

            if parameters.get('pressType') == 'LPP':
                recovery_percentage = 100.0 * (n_1_bd - n_1_evac) / (n_1_adsin if n_1_adsin != 0 else 1.0)
            else:
                recovery_percentage = 100.0 * (n_1_bd - n_1_evac) / (n_1_pres - n_1_evac + n_1_adsin if (n_1_pres - n_1_evac + n_1_adsin) != 0 else 1.0)
            purity_percentage = 100.0 * (n_1_bd - n_1_evac) / ((n_1_bd - n_1_evac) + (n_2_bd - n_2_evac) if ((n_1_bd - n_1_evac) + (n_2_bd - n_2_evac)) != 0 else 1.0)

            cycle_time = (parameters['t_ads'] + parameters['t_blo'] + parameters['t_evac'] + parameters['t_press'])
            productivity = (n_1_bd - n_1_evac) / (parameters['V_column'] * cycle_time if parameters['V_column'] * cycle_time != 0 else 1.0)

            # Energy calculation
            grad_X2_q1 = np.gradient(X2[:, 1], dt)
            grad_X2_q2 = np.gradient(X2[:, 2], dt)
            Fout_bd = 0 - (1 - parameters['e_bed']) * parameters['V_column'] * parameters['rho_s'] * (grad_X2_q1 + grad_X2_q2) - (parameters['e_bed'] / (Rg * X2[:, 3])) * parameters['dPdt_blo'](t2) * parameters['V_column']
            Fout_bd[Fout_bd < 0] = 0
            vout_bd = Fout_bd * Rg * X2[:, 3] / (parameters['P_blo'](t2)) / A_in

            Fout_evac = 0 - (1 - parameters['e_bed']) * parameters['V_column'] * parameters['rho_s'] * (np.gradient(X3[:, 1], dt) + np.gradient(X3[:, 2], dt)) - (parameters['e_bed'] / (Rg * X3[:, 3])) * parameters['dPdt_evac'](t3) * parameters['V_column']
            vout_evac = Fout_evac * Rg * X3[:, 3] / (parameters['P_evac'](t3)) / A_in

            Fin_pres = (1 - parameters['e_bed']) * parameters['V_column'] * (np.gradient(X4[:, 1], dt) + np.gradient(X4[:, 2], dt)) + (parameters['e_bed'] / (Rg * X4[:, 3])) * parameters['dPdt_press'](t4) * parameters['V_column']
            vin_pres = Fin_pres * Rg * X4[:, 3] / (parameters['P_press'](t4)) / A_in

            eta_bd = 0.8 * (19.55 * parameters['P_blo'](t2) * 1e-5 / (1 + 19.55 * parameters['P_blo'](t2) * 1e-5))
            eta_evac = 0.8 * (19.55 * parameters['P_evac'](t3) * 1e-5 / (1 + 19.55 * parameters['P_evac'](t3) * 1e-5))
            eta_press = 0.8
            eta_ads = 0.8

            EC_BD = np.trapz(1.0 / eta_bd * vout_bd * A_in * 1.0 * parameters['P_blo'](t2) * (1.4 / 0.4) * ((1e5 / np.minimum(1e5, parameters['P_blo'](t2))) ** (0.4 / 1.4) - 1), t2)
            EC_EVAC = np.trapz(1.0 / eta_evac * vout_evac * A_in * 1.0 * parameters['P_evac'](t3) * (1.4 / 0.4) * ((1e5 / np.minimum(1e5, parameters['P_evac'](t3))) ** (0.4 / 1.4) - 1), t3)
            EC_PRES = np.trapz(1.0 / eta_press * vin_pres * A_in * 1.0 * parameters['P_press'](t4) * (1.4 / 0.4) * ((np.maximum(1e5, parameters['P_press'](t4)) / 1e5) ** (0.4 / 1.4) - 1), t4)
            EC_ADS = (t1[-1] * 1.0 / eta_ads * parameters['v_in'] * A_in * 1.0 * parameters['p_H'] * (1.4 / 0.4) * ((np.maximum(1e5, parameters['p_H']) / 1e5) ** (0.4 / 1.4) - 1))

            SEC = (EC_PRES + EC_ADS + EC_BD + EC_EVAC) / ((n_1_bd - n_1_evac) * 0.044 if (n_1_bd - n_1_evac) * 0.044 != 0 else 1.0)

            recovery_percentageValues.append(recovery_percentage)
            purity_percentageValues.append(purity_percentage)
            productivity_Values.append(productivity)
            SEC_Values.append(SEC)

            process_indicators = np.vstack([purity_percentageValues, recovery_percentageValues, productivity_Values, SEC_Values])

            # pdb.set_trace()
            if cycle > 5:
                for i in range(4):
                    for k in range(5):
                        prev = process_indicators[i, cycle - 5]
                        curr = process_indicators[i, cycle - k - 1]
                        if prev != 0 and abs(100 * (curr - prev) / prev) <= 0.5:
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

        if parameters.get('OptType') == 'Const':
            KPIs = [(-productivity + phi_pen[0]), 10 * (SEC * 2.77778e-7 + phi_pen[1])]
        else:
            KPIs = [-recovery_percentage, -purity_percentage]

    except Exception:
        # On error return large KPIs like the original MATLAB implementation
        KPIs = [1e5, 1e5]
        SEC = 1e5
        
    
    #Python 3: 
    # print(datetime.now() - startTime)


    # Optionally plot (basic reproduction of MATLAB plotting)
    if parameters.get('outputType') == 'plot':
        import matplotlib.pyplot as plt

        t_cycle = np.concatenate([t1, t2 + t1[-1], t3 + t1[-1] + t2[-1], t4 + t1[-1] + t2[-1] + t3[-1]])
        X_cycle = np.vstack([X1, X2, X3, X4])
        P1 = np.full_like(t1, parameters['P_ads'](t_ads))
        P2 = parameters['P_blo'](t2)
        P3 = parameters['P_evac'](t3)
        P4 = parameters['P_press'](t4)
        P_cycle = np.hstack([P1, P2, P3, P4])

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
