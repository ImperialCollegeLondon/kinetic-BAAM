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


def run_nsga(parameters: Dict[str, Any], ngens: int = 90, pop_size: int = 10, use_pymoo: bool = True, n_cores: int | None = None) -> Tuple[np.ndarray, np.ndarray]:
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
    pt = parameters.get('processType', 'PVSA')

    # Build lb/ub exactly as in the MATLAB script (these are in the same
    # encoding expected by the outputs code when parameters['outputType']=='opt')
    if pt == 'PVSA':
        lb = np.array([0.3 * 0.37, 0.13e5, 40.0, 30.0, 30.0, 1e5])
        ub = np.array([3.0 * 0.37, 5e5, 300.0, 300.0, 300.0, 10e5])
        A = np.array([[0, -1, 0, 0, 0, 0], [0, 1, 0, 0, 0, -1]], dtype=float)
        b = np.array([0.0, 0.0])
    elif pt == 'VSA':
        if parameters.get('amine', False):
            lb = np.array([0.3 * 0.37, 0.021e5, 100.0, 20.0, 30.0, 0.02e5])
            ub = np.array([3.0 * 0.37, 0.9e5, 3e4, 1200.0, 3.5e4, 0.5e5])
        else:
            lb = np.array([0.3 * 0.37, 0.021e5, 30.0, 30.0, 30.0, 0.02e5])
            ub = np.array([3.0 * 0.37, 0.9e5, 200.0, 300.0, 300.0, 0.5e5])
        A = np.array([[0, -1, 0, 0, 0, 1]], dtype=float)
        b = np.array([0.0])
    elif pt == 'AdsorbentVSA':
        lb = np.array([np.log10(0.02e5), 100, 40, 0.5, 1e-6, 10e3, 10e3])
        ub = np.array([np.log10(0.9e5), 3000, 300, 8, 1e-3, 45e3, 45e3])
        A = np.array([[0, 0, 0, 0, 0, -1, 1]])
        b = np.array([0.0])
    elif pt == 'AdsorbentPVSA':
        lb = np.array([np.log10(0.02e5), 100, 40, 0.5, 1e-6, 10e3, 10e3, np.log10(1e5)])
        ub = np.array([np.log10(10e5), 3000, 300, 8, 1e-3, 45e3, 45e3, np.log10(10e5)])
        A = np.array([[1, 0, 0, 0, 0, -1, 1, -1]])
        b = np.array([0.0])
    elif pt == 'AdsorbentVSAb0':
        lb = np.array([np.log10(0.02e5), 100, 40, 0.5, 1e-6, 10e3, 10e3])
        ub = np.array([np.log10(0.9e5), 3000, 300, 8, 1e-3, 45e3, 45e3])
        A = np.array([[1, 0, 0, 0, -1, 1, 0]])
        b = np.array([0.0])
    elif pt == 'AdsorbentPVSAb0':
        lb = np.array([np.log10(0.02e5), 100, 40, 0.5, 1e-6, 1e-6, 10e3, np.log10(1e5)])
        ub = np.array([np.log10(10e5), 3000, 300, 8, 1e-3, 1e-3, 45e3, np.log10(10e5)])
        A = np.array([[1, 0, 0, 0, -1, 1, 0, -1]])
        b = np.array([0.0])
    elif pt == 'Resin':
        lb = np.array([np.log10(0.021e5), 1000.0, 20.0, 500.0, np.log10(0.05e5), 0.005])
        ub = np.array([np.log10(0.975e5), 30000.0, 1000.0, 45000.0, np.log10(0.2e5), 1.0])
        A = np.array([[-1, 0, 0, 0, 1, 0]], dtype=float)
        b = np.array([0.0])
    else:
        raise ValueError(f"Unsupported processType: {pt}")

    n_vars = lb.size
    A = np.atleast_2d(A) if 'A' in locals() else np.zeros((0, n_vars))
    b = np.asarray(b) if 'b' in locals() else np.array([])

    # filename for saving
    timestamp = datetime.now().strftime('%d%m%y%H%M')
    fname = f"{parameters.get('adsorbentName','ads')}_{pt}_{parameters.get('pressType','FP')}_{parameters.get('OptType','Unc')}_{parameters.get('modelType','model')}_{timestamp}"
    parameters['fileName'] = fname
    parameters['outputType'] = 'opt'
    parameters['xRef'] = np.ones_like(ub, dtype=float)

    # Resolve n_cores from explicit argument or parameters dict
    if n_cores is None:
        try:
            n_cores = int(parameters.get('n_cores', 1))
        except Exception:
            print('n_cores = 1')
            n_cores = 1

    # Try using pymoo if requested
    if use_pymoo:
        # try:
        from pymoo.algorithms.moo.nsga2 import NSGA2
        from pymoo.operators.sampling.lhs import LHS
        # from pymoo.operators.mutation.bitflip import BitflipMutation
        # from pymoo.operators.crossover.hux import HUX        
        from pymoo.optimize import minimize
        from pymoo.core.problem import ElementwiseProblem
        from pymoo.algorithms.soo.nonconvex.de import DE
        from pymoo.core.problem import StarmapParallelization
        from multiprocessing.pool import ThreadPool
        # import multiprocessing
        # from multiprocessing.pool import ThreadPool
        # from multiprocessing import Pool
        # from multiprocessing.pool import ThreadPool
        # from pymoo.core.problem import StarmapParallelization
        # from pymoo.algorithms.soo.nonconvex.ga import GA
        # from pymoo.optimize import minimize
        class _Problem(ElementwiseProblem):
            def __init__(self):
                super().__init__(n_var=n_vars, n_obj=2, n_constr=A.shape[0], xl=lb, xu=ub,elementwise_runner=runner)

            def _evaluate(self, x, out, *args, **kwargs):
                # x is a 1D array in decision space (same encoding as MATLAB)
                kpis = kbaam_outputs_nonisothermal(parameters.copy(), thetaIn=x.tolist())
                # objective is same as MATLAB return (2-element list)
                out['F'] = np.array(kpis, dtype=float)
                # linear inequality constraints A x <= b -> G(x) = A x - b <= 0
                if A.size > 0:
                    g = A.dot(np.asarray(x)) - b
                    out['G'] = np.atleast_1d(g)
                    
        N_THREADS = n_cores
        pool = ThreadPool(N_THREADS)
        runner = StarmapParallelization(pool.starmap)
        problem = _Problem()        
                    
        # Configure algorithm. Many pymoo components accept `n_jobs` or
        # `n_processors` in newer versions; pass n_cores where applicable.
        algo_kwargs = dict(pop_size=pop_size, sampling=LHS(smooth=True, iterations=pop_size))
        try:
            # some pymoo versions accept `n_jobs` on the algorithm
            algo = NSGA2(**algo_kwargs, n_jobs=n_cores) if n_cores and n_cores > 1 else NSGA2(**algo_kwargs)
        except TypeError:
            # fallback if the NSGA2 constructor doesn't accept n_jobs
            algo = NSGA2(**algo_kwargs)

        # minimize supports `n_proc` via `minimize(..., n_processes=...)` in some
        # environments; attempt to pass it and otherwise call without it.
        try:
            print('running parallel')
            res = minimize(problem, algo, ('n_gen', ngens), verbose=True, n_proc=n_cores)
        except TypeError:
            res = minimize(problem, algo, ('n_gen', ngens), verbose=True)
            
        pool.close()
        X = res.X
        F = res.F

        # Save results
        os.makedirs('matFiles', exist_ok=True)
        scipy.io.savemat(os.path.join('matFiles', f"{fname}.mat"), {'x': X, 'fval': F})
        return X, F
        # except Exception:
        #     print('pymoo not available or failed; falling back to random search.')

    # # Random-search fallback: draw R = pop_size * ngens samples uniformly
    # R = int(pop_size * max(1, ngens))
    # # Latin hypercube sampling gives better coverage; use scipy.stats.qmc
    # try:
    #     from scipy.stats import qmc
    #     sampler = qmc.LatinHypercube(d=n_vars)
    #     X_unit = sampler.random(n=R)
    # except Exception:
    #     X_unit = np.random.rand(R, n_vars)

    # X_samples = lb + X_unit * (ub - lb)

    # # Apply linear inequality constraints A x <= b to sampled candidates
    # if A.size > 0:
    #     # A shape (m, n_vars), X_samples shape (R, n_vars)
    #     lhs = (A @ X_samples.T).T  # shape (R, m)
    #     # ensure b has shape (m,)
    #     b_vec = np.asarray(b).reshape(-1)
    #     if b_vec.size == 1:
    #         mask_valid = lhs[:, 0] <= b_vec.item()
    #     else:
    #         mask_valid = np.all(lhs <= b_vec[None, :], axis=1)
    #     X_samples = X_samples[mask_valid]

    # F_samples = np.zeros((len(X_samples), 2), dtype=float)
    # # Evaluate samples in parallel if requested
    # if n_cores and n_cores > 1:

    #     def _eval_theta(theta_row):
    #         theta = theta_row.tolist()
    #         try:
    #             kpis = kbaam_outputs_nonisothermal(parameters.copy(), thetaIn=theta)
    #             return np.array(kpis, dtype=float)
    #         except Exception as e:
    #             # return large penalty values on failure to keep candidate out
    #             return np.array([1e6, 1e6], dtype=float)

    #     with Pool(processes=n_cores) as p:
    #         results = p.map(_eval_theta, [X_samples[i, :] for i in range(len(X_samples))])
    #     F_samples = np.vstack(results)
    # else:
    #     for i in range(len(X_samples)):
    #         theta = X_samples[i, :]
    #         kpis = kbaam_outputs_nonisothermal(parameters.copy(), thetaIn=theta.tolist())
    #         F_samples[i, :] = np.array(kpis, dtype=float)
    #         if (i + 1) % 100 == 0:
    #             print(f'Evaluated {i+1}/{len(X_samples)} samples')

    # # Determine Pareto front (minimization)
    # mask = _pareto_mask(F_samples)
    # X_pareto = X_samples[mask]
    # F_pareto = F_samples[mask]

    # # save results
    # os.makedirs('matFiles', exist_ok=True)
    # scipy.io.savemat(os.path.join('matFiles', f"{fname}.mat"), {'x': X_samples, 'fval': F_samples, 'x_pareto': X_pareto, 'f_pareto': F_pareto})

    # return X_pareto, F_pareto


if __name__ == '__main__':
    params = create_parameters()
    print('Starting run_nsga (random-search fallback if pymoo not installed)...')
    # Example: use 4 cores for the random-search fallback
    Xp, Fp = run_nsga(params, ngens=10, pop_size=20, use_pymoo=False, n_cores=4)
    print('Pareto front size:', Xp.shape[0])
    print('Saved results to matFiles/')
