"""createParameters.py

Imperial College London, Multiphase Systems Laboratory
Year: 2025
Author: Hassan Azzan (HA)

Purpose:
        Create and return a parameters dictionary equivalent to the MATLAB
        ``createParameters`` function. The returned dict mirrors the MATLAB
        ``parameters`` struct (adsorbent properties, cycle parameters, models
        settings) and includes lightweight DSL and LDFCoefficient fallbacks so
        the translated code can run smoke tests immediately.

Functions:
        - create_parameters() -> dict
        - save_parameters_mat(parameters, folder='AdsorbentFiles') (available when
            running the module directly) — saves a MATLAB .mat file named by
            ``parameters['adsorbentName']`` into ``AdsorbentFiles/``.

Notes:
        - Units and semantics follow the MATLAB implementation; see the original
            MATLAB header comments in ``createParameters.m`` for detailed field
            descriptions.
"""

from __future__ import annotations
from typing import Dict, Any
import warnings
import numpy as np
import os
import scipy.io


def create_parameters() -> Dict[str, Any]:
    """Return parameters dict matching MATLAB `createParameters`.

    The returned dict mirrors the MATLAB `createParameters` struct. See the
    README or function docstring comments for notes on units and expected
    semantics. Defaults below are chosen to match the original MATLAB script
    so the translated code can run smoke tests immediately.
    """

    parameters: Dict[str, Any] = {}

    # Adsorbent and bed properties
    parameters['adsorbentName'] = "Z13X_AW_2022"
    parameters['qsb_1'] = 3.09
    parameters['qsd_1'] = 2.54
    parameters['bo_1'] = 8.65e-7
    parameters['do_1'] = 2.63e-8
    parameters['delUb_1'] = -3.66e4
    parameters['delUd_1'] = -3.57e4
    parameters['qsb_2'] = 5.84
    parameters['qsd_2'] = 0.0
    parameters['bo_2'] = 2.50e-7
    parameters['do_2'] = 0.0
    parameters['delUb_2'] = -1.58e4
    parameters['delUd_2'] = 0.0

    parameters['cp_g'] = 30.7
    parameters['cp_a'] = parameters['cp_g']
    parameters['cp_w'] = 502.0
    parameters['cp_s'] = 1070.0
    parameters['rho_s'] = 1130.0
    parameters['rho_w'] = 7800.0
    parameters['rp'] = 1e-3
    parameters['V_column'] = 0.066
    parameters['t_wall'] = 0.003
    parameters['e_bed'] = 0.37
    parameters['h_in'] = 8.6
    parameters['h_out'] = 2.5
    parameters['epsilon_p'] = 0.35
    parameters['Dm'] = 1.6e-5
    parameters['tau'] = 3.0

    # Column geometry and wall properties (defaults added to match ODE use)
    parameters['r_in'] = 0.02
    parameters['r_out'] = 0.03
    parameters['L'] = 1.0
    parameters['Theat'] = 0

    # Cycle properties
    parameters['v_in'] = 0.5
    parameters['y1_in'] = 0.15
    parameters['T_feed'] = 298.0
    parameters['p_H'] = 1.00e5
    parameters['p_I'] = 0.30e5
    parameters['p_L'] = 0.01e5
    parameters['t_ads'] = 300.0
    parameters['t_blo'] = 100.0
    parameters['t_evac'] = 300.0
    parameters['t_press'] = 20.0
    parameters['heating'] = False
    parameters['pressType'] = "FP"
    parameters['processType'] = "PVSA"
    parameters['lambda'] = 0.5

    # Model properties
    parameters['modelType'] = "nonisothermal"
    parameters['OptType'] = "Unc"
    parameters['outputType'] = "plot"
    # parameters['outputType'] = "opt"

    # Derived / required fields used by the Python ODE implementation
    # F_in may be computed by the outputs routine from v_in; leave None by
    # default so the outputs code computes it consistently.
    parameters['F_in'] = None
    parameters['y1_LPP'] = parameters['y1_in']

    # Default simple pressure schedules and derivatives (vector-aware)
    # These default functions accept scalar or numpy-array t and return a
    # scalar or array respectively.
    parameters['P_ads'] = lambda tt: (np.asarray(parameters['p_H']) if np.asarray(tt).shape == () else np.full(np.shape(np.asarray(tt)), parameters['p_H']))
    parameters['P_blo'] = lambda tt: (parameters['p_I'] + (parameters['p_H'] - parameters['p_I']) * np.exp(-parameters['lambda'] * np.asarray(tt)))
    parameters['P_evac'] = lambda tt: (parameters['p_L'] + (parameters.get('P_initL', parameters['p_I']) - parameters['p_L']) * np.exp(-parameters['lambda'] * np.asarray(tt)))
    parameters['P_press'] = lambda tt: (parameters['p_H'] + (parameters.get('P_initR', parameters['p_H']) - parameters['p_H']) * np.exp(-3 * np.asarray(tt)))

    parameters['dPdt_blo'] = lambda tt: -parameters['lambda'] * (parameters['p_H'] - parameters['p_I']) * np.exp(-parameters['lambda'] * np.asarray(tt))
    parameters['dPdt_evac'] = lambda tt: -parameters['lambda'] * (parameters['p_I'] - parameters['p_L']) * np.exp(-parameters['lambda'] * np.asarray(tt))
    parameters['dPdt_press'] = lambda tt: -3.0 * (parameters.get('P_initR', parameters['p_H']) - parameters['p_H']) * np.exp(-3 * np.asarray(tt))

    # Attach DSL and LDFCoefficient implementations from local modules so
    # the returned parameters dict is immediately usable. If you prefer to
    # override these, assign different callables to these keys.
    _DSL = None
    _LDFCoefficient = None

    try:
        # try top-level imports first
        from DSL import DSL as _DSL  # type: ignore
        from LDFCoefficient import LDFCoefficient as _LDFCoefficient  # type: ignore
    except Exception:
        try:
            # package-style import (if this module is used as a package)
            from .DSL import DSL as _DSL  # type: ignore
            from .LDFCoefficient import LDFCoefficient as _LDFCoefficient  # type: ignore
        except Exception:
            _DSL = None
            _LDFCoefficient = None

    # Provide lightweight fallbacks when the real implementations are not
    # importable. These are intentionally simple and only aim to keep the
    # translated code runnable for smoke tests; they are NOT scientifically
    # accurate replacements.
    if _DSL is None:
        warnings.warn("Could not import DSL; using lightweight fallback (for testing only).", UserWarning)

        def _DSL(P, y1, T, qsb_1, qsd_1, qsb_2, qsd_2, bo_1, do_1, bo_2, do_2, delUb_1, delUd_1, delUb_2, delUd_2):
            P_arr = np.asarray(P)
            y1_arr = np.asarray(y1)
            # scalar case
            if P_arr.shape == () and y1_arr.shape == ():
                Pnorm = float(P_arr) / 1e5
                q1 = qsb_1 + qsd_1 * 0.5 * (1.0 + float(y1_arr)) * (1.0 + Pnorm)
                q2 = qsb_2 + qsd_2 * 0.5 * (1.0 - float(y1_arr)) * (1.0 + Pnorm)
                return q1, q2
            # broadcast y1 and P to same shape
            if P_arr.shape != y1_arr.shape:
                try:
                    P_b = np.broadcast_to(P_arr, y1_arr.shape)
                except Exception:
                    P_b = np.full(y1_arr.shape, P_arr.item() if P_arr.shape == () else P_arr)
            else:
                P_b = P_arr
            Pnorm = P_b / 1e5
            q1 = qsb_1 + qsd_1 * 0.5 * (1.0 + y1_arr) * (1.0 + Pnorm)
            q2 = qsb_2 + qsd_2 * 0.5 * (1.0 - y1_arr) * (1.0 + Pnorm)
            return q1, q2

    if _LDFCoefficient is None:
        warnings.warn("Could not import LDFCoefficient; using lightweight fallback (for testing only).", UserWarning)

        def _LDFCoefficient(P, y1, T, q1_star, q2_star, parameters_local):
            T_arr = np.asarray(T)
            k1_base = parameters_local.get('bo_1', 1e-6)
            k2_base = parameters_local.get('bo_2', 1e-6)
            scale = 1.0 + 0.01 * (T_arr - 300.0)
            k1 = k1_base * scale
            k2 = k2_base * scale
            return k1, k2

    parameters['DSL'] = _DSL
    parameters['LDFCoefficient'] = _LDFCoefficient

    return parameters


if __name__ == '__main__':
    # Simple smoke creation
    params = create_parameters()
    print('Created parameters with keys:', sorted(params.keys()))
    # Save the parameters dict to AdsorbentFiles/<adsorbentName>.mat
    def save_parameters_mat(parameters: Dict[str, Any], folder: str = 'AdsorbentFiles', overwrite: bool = False) -> str:
        """Save the parameters dict to a .mat file named by parameters['adsorbentName'].

        Returns the path to the saved file.
        """
        ads_name = parameters.get('adsorbentName')
        if not ads_name:
            raise ValueError("parameters must contain 'adsorbentName' to determine output filename")

        if not os.path.exists(folder):
            os.makedirs(folder, exist_ok=True)

        if not ads_name.endswith('.mat'):
            filename = f"{ads_name}.mat"
        else:
            filename = ads_name

        path = os.path.join(folder, filename)
        if os.path.exists(path) and not overwrite:
            raise FileExistsError(f"File already exists: {path} (set overwrite=True to replace)")

        # Prepare a copy of parameters suitable for saving to .mat.
        def _prune_for_mat(obj):
            # dict -> recurse
            if isinstance(obj, dict):
                out = {}
                for kk, vv in obj.items():
                    # skip callables and module/class objects
                    if callable(vv):
                        continue
                    # skip bound methods / types
                    if isinstance(vv, type) or hasattr(vv, '__module__') and vv.__module__ == 'builtins' and isinstance(vv, object) and False:
                        pass
                    pr = _prune_for_mat(vv)
                    out[kk] = pr
                return out
            # None -> empty array
            if obj is None:
                return np.array([])
            # numpy scalar -> python scalar
            try:
                if isinstance(obj, np.generic):
                    return obj.item()
            except Exception:
                pass
            # numpy arrays, lists, tuples -> keep as-is (savemat handles these)
            if isinstance(obj, (np.ndarray, list, tuple)):
                return obj
            # primitives
            if isinstance(obj, (int, float, str, bool)):
                return obj
            # fallback: try to convert to string
            try:
                return str(obj)
            except Exception:
                return None

        pruned = _prune_for_mat(parameters)

        # Save under the MATLAB variable name 'parameters' so MATLAB loads a struct
        scipy.io.savemat(path, {'parameters': pruned})
        return path

    try:
        saved = save_parameters_mat(params)
        print(f"Saved parameters to: {saved}")
    except FileExistsError as e:
        print(str(e))