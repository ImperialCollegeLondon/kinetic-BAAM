from __future__ import annotations
import os
import numpy as np
import matplotlib.pyplot as plt
from typing import List
from find_pareto_frontier import find_pareto_frontier
import os
import sys
import git
import pdb

def get_git_root(path):

        git_repo = git.Repo(path, search_parent_directories=True)
        git_root = git_repo.git.rev_parse("--show-toplevel")
        print(git_root)


def plot_energy_vs_productivity(files: List[str], labels: List[str] | None = None):
    """Plot Energy (kWh/tonne) vs Productivity (mol/m^3/s) from rawData text outputs.

    Only points satisfying purity>95 and recovery>90 are considered for Pareto set.
    """
    # Absolute path to the repository root (parent of this pythonPort directory)

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

    if labels is None:
        labels = [os.path.splitext(os.path.basename(f))[0] for f in files]
    

    color_vals = ['b', 'r', 'm', 'none']
    marker_vals = ['s', 's', 'o', 'o']
    
    
    plt.figure()
    for i, path in enumerate(files):
        # pdb.set_trace()
        try:
            data = np.loadtxt(os.path.join(REPO_ROOT,'rawData/',path))
        except Exception:
            data = np.loadtxt(os.path.join(REPO_ROOT,'pythonPort/rawData/',path))
        
        
        mask = (data[:, 8-1] > 95) & (data[:, 9-1] > 90)
        ab = data[mask]
        input_pairs = []
        for row in ab:
            input_pairs.append([1.0 / max(row[10-1], 1e-12), row[11-1]])
        input_arr = np.asarray(input_pairs)
        flag, pareto = find_pareto_frontier(input_arr)
        outPareto = np.column_stack([1.0 / np.maximum(pareto[:, 0], 1e-12), pareto[:, 1]]) if pareto.size else np.empty((0, 2))
        if outPareto.size:
            plt.scatter(outPareto[:, 0], outPareto[:, 1] * 2.77778e-7 * 1e3, s=50, marker=marker_vals[i % len(marker_vals)], edgecolors='none', linewidths=0.8, facecolors=color_vals[i % len(color_vals)], label=labels[i])
        plt.box(True)
        plt.grid(True, linewidth=1)
        plt.ylabel('E_T [kWh/tonne]')
        plt.xlabel('CO2 Productivity [mol/m3/s]')
        plt.ylim([0, 1400])
        plt.legend(loc='lower right')
        plt.tight_layout()
    plt.show()
