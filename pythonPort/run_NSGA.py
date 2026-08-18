"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Imperial College London, United Kingdom
Multiphase Systems Laboratory
Year:     2025
MATLAB:   R2024a
Authors:  Hassan Azzan (HA)

Purpose:
Function that takes parameters as inputs and carries out multiobjective
process optimization for the system defined in parameters.

Last modified:
- 2025-10-08, HA: Add reverse engineering method
- 2025-09-21, HA: Initial creation

Input arguments:
    - parameters: contains adsorbent properties and process parameters

Output arguments:
    - x, fval equivalents in the Python return value (X_pareto, F_pareto)

Dependencies:
    - kBAAM_Outputs_nonIsothermal.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"""

from __future__ import annotations
import numpy as np
from typing import Dict, Any, Tuple
from datetime import datetime
import os
import scipy.io

from kBAAM_Outputs_nonIsothermal import kbaam_outputs_nonisothermal
from createParameters import create_parameters


def _pareto_mask(costs: np.ndarray) -> np.ndarray:
    """Return boolean mask of Pareto-efficient points (minimization).

    costs: shape (n_points, n_objectives)
    """
    n_points = costs.shape[0]
    is_efficient = np.ones(n_points, dtype=bool)
    for i in range(n_points):
        if not is_efficient[i]:
            continue
        # any point that strictly dominates i -> mark i as not efficient
        better = np.all(costs <= costs[i], axis=1) & np.any(costs < costs[i], axis=1)
        better[i] = False
        if np.any(better):
            is_efficient[i] = False
        else:
            # remove points dominated by i
            dominated = np.all(costs >= costs[i], axis=1) & np.any(costs > costs[i], axis=1)
            is_efficient[dominated] = False
    return is_efficient


def run_nsga(parameters: Dict[str, Any], ngens: int = 120, pop_size: int = 300, use_pymoo: bool = True, n_cores: int | None = None) -> Tuple[np.ndarray, np.ndarray]:
    """Run NSGA-style optimisation for the provided `parameters`.

    Returns (X_pareto, F_pareto) arrays containing decision vectors and
    objective values on the discovered Pareto front.
    
    New parameter / inputs flag:
        - n_cores: number of worker processes to use for parallel evaluation
                   of candidate solutions when available. If None, the
                   function will look for ``parameters['n_cores']`` and
                   otherwise default to 1 (serial). To request 10 cores you
                   can either call ``run_nsga(..., n_cores=10)`` or set
                   ``parameters['n_cores']=10`` before calling.
    """
    parameters.setdefault('twoNode', 0)
    pt = parameters.get('processType', 'PVSA')

    # lb / ub / A / b and xRef match run_NSGA.m exactly
    if pt == 'PVSA':
        lb = np.array([0.2 * 0.37, 0.13e5,  30.0, 30.0, 30.0,  1e5])
        ub = np.array([2.0 * 0.37, 9e5,    300.0, 300.0, 300.0, 10e5])
        A  = np.array([[0, -1, 0, 0, 0,  0],
                       [0,  1, 0, 0, 0, -1]], dtype=float)
        b  = np.array([0.0, 0.0])
        parameters['xRef'] = np.ones(len(lb))
    elif pt == 'VSA':
        if parameters.get('amine', False):
            lb = np.array([0.2 * 0.37, 0.021e5, 100.0,   20.0,   30.0, 0.02e5])
            ub = np.array([3.0 * 0.37, 0.9e5,   3e4,   1200.0, 3.5e4,  0.5e5])
        else:
            lb = np.array([0.3 * 0.37, 0.021e5,  30.0,  30.0,  30.0, 0.02e5])
            ub = np.array([2.0 * 0.37, 0.9e5,   300.0, 300.0, 300.0,  0.5e5])
        A  = np.array([[0, -1, 0, 0, 0, 1]], dtype=float)
        b  = np.array([0.0])
        parameters['xRef'] = ub.copy()
    elif pt == 'AdsorbentVSA':
        lb = np.array([np.log10(0.02e5), 100,  40, 0.5, 1e-6, 10e3, 10e3])
        ub = np.array([np.log10(0.9e5),  3000, 300,  8, 1e-3, 45e3, 45e3])
        A  = np.array([[0, 0, 0, 0, 0, -1, 1]], dtype=float)
        b  = np.array([0.0])
        parameters['xRef'] = ub.copy()
    elif pt == 'AdsorbentPVSA':
        # 12-variable encoding: [v_in, p_I, t_ads, t_blo, t_evac, p_H,
        #                        qsb, log10(b01), log10(b02), delU1, delU2, rho_s]
        lb = np.array([0.2*0.37, 0.13e5, 30, 30, 30,  1e5,  2, -7, -7, 20e3,  5e3,  300])
        ub = np.array([2.0*0.37,  9e5,  300, 300, 300, 10e5, 10, -2, -2, 50e3, 30e3, 1600])
        A  = np.array([[0, -1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0],
                       [0,  1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0]], dtype=float)
        b  = np.array([0.0, 0.0])
        parameters['xRef'] = np.ones(len(lb))
        ngens    = 90
        pop_size = 300
    elif pt == 'AdsorbentPVSADSL':
        # 15-variable encoding: [v_in, p_I, t_ads, t_blo, t_evac, p_H, qsb, qsd,
        #                        log10(b01), log10(d01), log10(b02/d02),
        #                        delUb1, delUd1, delU2, rho_s]
        lb = np.array([0.2*0.37, 0.13e5, 30, 30, 30,  1e5,  1, 0.1, -7, -7, -7, 20e3, 20e3,  5e3,  600])
        ub = np.array([2.0*0.37,  9e5,  300, 300, 300, 10e5, 10,  10, -2, -2, -2, 50e3, 50e3, 20e3, 1600])
        A  = np.array([[0, -1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0,  1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=float)
        b  = np.array([0.0, 0.0])
        parameters['xRef'] = np.ones(len(lb))
    elif pt == 'AdsorbentVSAb0':
        lb = np.array([np.log10(0.02e5), 100,  40, 0.5, 1e-6, 10e3, 10e3])
        ub = np.array([np.log10(0.9e5),  3000, 300,  8, 1e-3, 45e3, 45e3])
        A  = np.array([[1, 0, 0, 0, -1, 1, 0]], dtype=float)
        b  = np.array([0.0])
        parameters['xRef'] = ub.copy()
    elif pt == 'AdsorbentPVSAb0':
        lb = np.array([np.log10(0.02e5), 100,  40, 0.5, 1e-6, 1e-6, 10e3, np.log10(1e5)])
        ub = np.array([np.log10(10e5),  3000, 300,  8, 1e-3, 1e-3, 45e3, np.log10(10e5)])
        A  = np.array([[1, 0, 0, 0, -1, 1, 0, -1]], dtype=float)
        b  = np.array([0.0])
        parameters['xRef'] = ub.copy()
    elif pt == 'Resin':
        lb = np.array([np.log10(0.021e5), 1000.0,   20.0,   500.0, np.log10(0.05e5), 0.005])
        ub = np.array([np.log10(0.975e5), 30000.0, 1000.0, 45000.0, np.log10(0.2e5),  1.0])
        A  = np.array([[-1, 0, 0, 0, 1, 0]], dtype=float)
        b  = np.array([0.0])
        parameters['xRef'] = ub.copy()
    else:
        raise ValueError(f"Unsupported processType: {pt}")

    n_vars = lb.size

    # filename for saving
    timestamp = datetime.now().strftime('%d%m%y%H%M')
    fname = f"{parameters.get('adsorbentName','ads')}_{pt}_{parameters.get('pressType','FP')}_{parameters.get('OptType','Unc')}_{parameters.get('modelType','model')}_{timestamp}"
    parameters['fileName'] = fname
    parameters['outputType'] = 'opt'

    # Resolve n_cores from explicit argument or parameters dict
    if n_cores is None:
        try:
            n_cores = int(parameters.get('n_cores', 1))
        except Exception:
            print('n_cores = 1')
            n_cores = 1

    def _run_one(params_copy, _ngens, _pop_size, _fname):
        """Run a single NSGA-II optimisation and save results."""
        if not use_pymoo:
            raise RuntimeError('use_pymoo must be True')
        from pymoo.algorithms.moo.nsga2 import NSGA2
        from pymoo.operators.sampling.lhs import LHS
        from pymoo.optimize import minimize
        from pymoo.core.problem import ElementwiseProblem, StarmapParallelization
        from multiprocessing.pool import ThreadPool

        class _Problem(ElementwiseProblem):
            def __init__(self):
                super().__init__(n_var=n_vars, n_obj=2, n_constr=A.shape[0],
                                 xl=lb, xu=ub, elementwise_runner=runner)

            def _evaluate(self, x, out, *args, **kwargs):
                kpis = kbaam_outputs_nonisothermal(params_copy.copy(), thetaIn=x.tolist())
                out['F'] = np.array(kpis, dtype=float)
                if A.size > 0:
                    out['G'] = np.atleast_1d(A.dot(np.asarray(x)) - b)

        pool   = ThreadPool(max(1, n_cores))
        runner = StarmapParallelization(pool.starmap)
        problem = _Problem()

        algo_kwargs = dict(pop_size=_pop_size, sampling=LHS(smooth=True, iterations=_pop_size))
        try:
            algo = NSGA2(**algo_kwargs, n_jobs=n_cores) if n_cores and n_cores > 1 else NSGA2(**algo_kwargs)
        except TypeError:
            algo = NSGA2(**algo_kwargs)

        try:
            print('running parallel')
            res = minimize(problem, algo, ('n_gen', _ngens), verbose=True, n_proc=n_cores)
        except TypeError:
            res = minimize(problem, algo, ('n_gen', _ngens), verbose=True)

        pool.close()
        os.makedirs('matFiles', exist_ok=True)
        scipy.io.savemat(os.path.join('matFiles', f"{_fname}.mat"),
                         {'x': res.X, 'fval': res.F})
        return res.X, res.F

    # ── Resin: iterate over adsorbent variants (mirrors MATLAB loops) ─────────
    if pt == 'Resin':
        if not parameters.get('fixResins', False):
            qsb_vals = np.arange(0.4, 3.6, 0.2)
            _ngens, _pop = 50, 120
            all_X, all_F = [], []
            for qsb in qsb_vals:
                p = parameters.copy()
                p['qsb_1'] = qsb
                p['fileName'] = fname
                X, F = _run_one(p, _ngens, _pop, f"{fname}_qsb{qsb:.2f}")
                all_X.append(X); all_F.append(F)
            return np.vstack(all_X), np.vstack(all_F)
        else:
            resin_vals = np.array([
                [0.454357067, 4.06469e-17, 0.0399],
                [0.51399776,  6.11738e-17, 0.034815],
                [0.984989599, 1.58041e-16, 0.019913],
                [1.336757211, 1.65654e-16, 0.01002],
                [0.781131901, 8.65318e-17, 0.0247],
                [1.174854923, 1.15484e-16, 0.0199],
                [2.37,        2.16419e-16, 0.0007],
            ])
            _ngens, _pop = 20, 300
            all_X, all_F = [], []
            for row in resin_vals:
                p = parameters.copy()
                p['qsb_1'], p['bo_1'], p['LDF'] = row[0], row[1], row[2]
                p['fileName'] = fname
                X, F = _run_one(p, _ngens, _pop, f"{fname}_qsb{row[0]:.3f}")
                all_X.append(X); all_F.append(F)
            return np.vstack(all_X), np.vstack(all_F)

    # ── All other process types ───────────────────────────────────────────────
    return _run_one(parameters, ngens, pop_size, fname)

if __name__ == '__main__':
    params = create_parameters()
    print('Starting run_nsga ...')
    Xp, Fp = run_nsga(params, ngens=120, pop_size=300, use_pymoo=True, n_cores=1)
    print('Pareto front size:', Xp.shape[0])
    print('Saved results to matFiles/')
