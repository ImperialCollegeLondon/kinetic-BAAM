"""Compatibility shim to expose `pythonPort` modules at the repository root.

Importing this module will prepend the `pythonPort` folder to sys.path and
re-export top-level modules so existing scripts that do e.g. ``from createParameters import ...``
continue to work.

Use: import kinetic_baam_shim  # ensures pythonPort is in sys.path
"""
from __future__ import annotations
import os
import sys

HERE = os.path.dirname(__file__)
PYTHONPORT = os.path.join(HERE, 'pythonPort')
if PYTHONPORT not in sys.path:
    sys.path.insert(0, PYTHONPORT)

# Re-export common modules for convenience
try:
    from createParameters import create_parameters  # type: ignore
    from run_NSGA import run_nsga  # type: ignore
    from kBAAM_Outputs_nonIsothermal import kbaam_outputs_nonisothermal  # type: ignore
    from auxillaryFunctionskBAAM import load_parameters_mat  # type: ignore
except Exception:
    # Ignore import errors here; modules can still be imported directly from pythonPort
    pass
