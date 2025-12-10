from __future__ import annotations
import sys
import os

# Absolute path to the repository root (parent of this pythonPort directory)
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
PYTHONPORT = os.path.join(REPO_ROOT, 'pythonPort')

# Prepend pythonPort and the repo root to sys.path if not already present
# so imports like `from createParameters import create_parameters` work when
# running scripts from the project root.
for p in (PYTHONPORT, REPO_ROOT):
    if p and p not in sys.path:
        sys.path.insert(0, p)

from auxillaryFunctionskBAAM import load_parameters_mat
from run_NSGA import run_nsga

# load (basename without .mat is OK)
p = load_parameters_mat('Z13X_AW_2022')
parameters = p["parameters"]

parameters["rp"] = 1e-3;
parameters["p_L"]= 0.02e5;
parameters["V_column"] = 0.0661;
parameters["outputType"] = "opt";
parameters["pressType"] = "FP";


###################################

parameters["adsorbentName"] = "Z13X";
parameters["modelType"] = "nonisothermal";
parameters["OptType"] = "Const";
parameters["processType"] = "PVSA";
# Request 10 cores for parallel evaluation (optional)
parameters['n_cores'] = 8

ngens = 90
pop_size = 200
Xp, Fp = run_nsga(parameters, ngens=ngens, pop_size=pop_size)
