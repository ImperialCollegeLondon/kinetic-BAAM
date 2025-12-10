"""auxillaryFunctionskBAAM.py

Helper utilities for loading adsorbent and parameter MAT files.

These helpers mirror the small utility behaviour used in the MATLAB
driver to locate and load parameter files from the repository's
``AdsorbentFiles/`` folder.
"""

from __future__ import annotations
import os
from typing import Dict, Any
import scipy.io
import sys

def load_adsorbent_mat(mat_name: str | None = None) -> Dict[str, Any] | None:
    """Load an adsorbent MAT-file from `AdsorbentFiles/`.

    Parameters
    ----------
    mat_name:
        The filename (e.g. 'Z13X_AW_2022.mat') or a base name without the
        '.mat' extension. If ``None`` the function returns ``None``.

    Returns
    -------
    dict or None
        The loaded mat contents (with ``simplify_cells=True``) or ``None`` if
        ``mat_name`` is ``None``.

    Raises
    ------
    FileNotFoundError
        If the requested file does not exist in the `AdsorbentFiles/` folder.
    RuntimeError
        If loading the MAT file fails for some reason.
    """
    if mat_name is None:
        return None
    
    # Absolute path to the repository root (parent of this pythonPort directory)
    REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    PYTHONPORT = os.path.join(REPO_ROOT, 'pythonPort')

    # Prepend pythonPort and the repo root to sys.path if not already present
    # so imports like `from createParameters import create_parameters` work when
    # running scripts from the project root.
    for p in (PYTHONPORT, REPO_ROOT):
        if p and p not in sys.path:
            sys.path.insert(0, p)

    # allow callers to pass a base name without the .mat extension
    if not mat_name.endswith('.mat'):
        mat_name = f"{mat_name}.mat"

    ads_dir = 'AdsorbentFiles'
    mat_path = os.path.join(ads_dir, mat_name)

    if not os.path.exists(mat_path):
        # try to provide a helpful list of available files
        try:
            available = sorted(os.listdir(ads_dir))
        except Exception:
            available = []
        raise FileNotFoundError(
            f"Requested adsorbent file '{mat_name}' not found in '{ads_dir}'. "
            f"Available files: {available}"
        )

    try:
        data = scipy.io.loadmat(mat_path, simplify_cells=True)
        return data
    except Exception as e:
        raise RuntimeError(f"Failed to load MAT file '{mat_path}': {e}")


def load_parameters_mat(mat_name: str | None = None) -> Dict[str, Any] | None:
    """Load a parameters MAT-file from `matFiles/` or `AdsorbentFiles/`.

    Tries `matFiles/` first, then `AdsorbentFiles/`. Returns the loaded dict or
    raises FileNotFoundError / RuntimeError analogous to `load_adsorbent_mat`.
    """
    
    
    ads_dir = 'AdsorbentFiles'
    mat_path = os.path.join(ads_dir, mat_name)

    if mat_name is None:
        return None

    if not mat_name.endswith('.mat'):
        mat_name = f"{mat_name}.mat"

    # only look in AdsorbentFiles/ (adsorbent parameter MATs are stored there)
    folder = 'AdsorbentFiles'
    path = os.path.join(folder, mat_name)
    if os.path.exists(path):
        try:
            raw = scipy.io.loadmat(path, simplify_cells=True)
            return _normalize_mat_params(raw)
        except Exception as e:
            raise RuntimeError(f"Failed to load MAT file '{path}': {e}")

    raise FileNotFoundError(f"Parameters MAT file '{mat_name}' not found in '{folder}/'.")


# helpers to convert scipy.loadmat outputs into plain Python types
def _convert_matobj(obj):
    """Recursively convert scipy.loadmat objects into Python-friendly types."""
    # dicts (already mapping-like)
    if isinstance(obj, dict):
        return _normalize_mat_params(obj)
    # numpy arrays and scalars
    try:
        import numpy as _np
        if isinstance(obj, _np.ndarray):
            if obj.dtype == object:
                try:
                    return [_convert_matobj(v) for v in obj.tolist()]
                except Exception:
                    return obj.tolist()
            if obj.size == 1:
                try:
                    return obj.ravel()[0].item()
                except Exception:
                    return obj.ravel()[0]
            return obj
    except Exception:
        pass
    return obj


def _normalize_mat_params(raw: Dict[str, Any]) -> Dict[str, Any]:
    """Strip MATLAB metadata keys and convert values to plain Python objects/dicts."""
    out: Dict[str, Any] = {}
    for k, v in raw.items():
        if k.startswith('__'):
            continue
        out[k] = _convert_matobj(v)
    return out
