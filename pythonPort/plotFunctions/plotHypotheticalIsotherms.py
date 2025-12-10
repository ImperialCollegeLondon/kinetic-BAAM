from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from typing import Sequence
import git
import sys

try:
    from ..DSL import DSL
except Exception:
    from pythonPort.DSL import DSL  # type: ignore


def plot_hypothetical_isotherms(Pvals: np.ndarray, theta: Sequence[float], T: float = 298.0):
    """Plot hypothetical isotherms using DSL for component 1 and 2.

    Parameters
    ----------
    Pvals : array of pressures [Pa]
    theta : sequence encoding [unused slots..., qsb, bo, delUb, delUb2] consistent with MATLAB script
    T : temperature [K]
    """
    # Absolute path to the repository root (parent of this pythonPort directory)
    repo = git.Repo('.', search_parent_directories=True)
    REPO_ROOT = repo.working_tree_dir
    # PYTHONPORT = os.path.join(REPO_ROOT, 'pythonPort')
    # pdb.set_trace()
    # Prepend pythonPort and the repo root to sys.path if not already present
    # so imports like `from createParameters import create_parameters` work when
    # running scripts from the project root.
    for p in (REPO_ROOT):
        if p and p not in sys.path:
            sys.path.append(p)
    # Interpret theta similar to MATLAB snippet (adapt as needed)
    # Here we assume: theta[-6:] = [v_in?, p_I?, t_ads?, t_blo?, t_evac?, p_H?] in main, but for hypo we used
    # qsb (both comps) at index 4, bo at 5, delUb1 at 6, delUb2 at 7 in the MATLAB snippet
    qsb = float(theta[4])
    bo = float(theta[5])
    delUb1 = float(theta[6])
    delUb2 = float(theta[7])

    params = dict(
        qsb_1=qsb, qsb_2=qsb,
        qsd_1=0.0, qsd_2=0.0,
        bo_1=bo, bo_2=bo,
        do_1=0.0, do_2=0.0,
        delUb_1=-delUb1, delUd_1=0.0,
        delUb_2=-delUb2, delUd_2=0.0,
    )

    q1, _ = DSL(Pvals, 1.0, T, **params)
    _, q2 = DSL(Pvals, 0.0, T, **params)

    fig, axs = plt.subplots(1, 2, figsize=(8, 3), constrained_layout=True)
    axs[0].plot(Pvals / 1e5, q1, color='r', lw=2)
    axs[0].set_xlabel('P [bar]')
    axs[0].set_ylabel('q_CO2 [mol/kg]')

    axs[1].plot(Pvals / 1e5, q2, color='b', lw=2)
    axs[1].set_xlabel('P [bar]')
    axs[1].set_ylabel('q_N2 [mol/kg]')
    return fig, axs
