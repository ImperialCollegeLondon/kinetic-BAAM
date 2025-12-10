from __future__ import annotations
import os
import numpy as np
import matplotlib.pyplot as plt
from typing import List
from find_pareto_frontier import find_pareto_frontier
import git
import sys

def plot_purity_recovery_from_txt(files: List[str], labels: List[str] | None = None):
    """Plot Recovery vs Purity from rawData text outputs.

    Each file is expected to be a whitespace-delimited text with columns as in
    the Outputs driver logging: purity in col 8, recovery in col 9, productivity col 10, SEC col 11.
    
    
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
    if labels is None:
        labels = [os.path.splitext(os.path.basename(f))[0] for f in files]

    color_vals = ['b', 'k', 'none', 'none']
    edge_vals = ['b', 'k', 'b', 'k']
    marker_vals = ['s', 's', 'o', 'o']

    plt.figure()
    for i, path in enumerate(files):
        try:
            data = np.loadtxt(os.path.join(REPO_ROOT,'rawData/',path))
        except Exception:
            data = np.loadtxt(os.path.join(REPO_ROOT,'pythonPort/rawData/',path))
        
        # Clamp/cleanup to MATLAB logic
        data[:, 10-1] = data[:, 10-1] * 0.044 * 3600 * 24
        input_pairs = []
        for row in data:
            row[8-1] = max(0, row[8-1])
            row[9-1] = max(0, row[9-1])
            if row[8-1] > 100 or row[9-1] > 100:
                row[8-1] = 0
                row[9-1] = 0
            input_pairs.append([1.0 / max(row[8-1], 1e-12), 1.0 / max(row[9-1], 1e-12)])
        input_arr = np.asarray(input_pairs)
        flag, pareto = find_pareto_frontier(input_arr)
        outPareto = np.column_stack([1.0 / np.maximum(pareto[:, 0], 1e-12), 1.0 / np.maximum(pareto[:, 1], 1e-12)]) if pareto.size else np.empty((0, 2))
        plt.scatter(data[flag, 8-1], data[flag, 9-1], s=50, marker=marker_vals[i % len(marker_vals)], edgecolors=edge_vals[i % len(edge_vals)], linewidths=0.8, facecolors=color_vals[i % len(color_vals)], label=labels[i])
    plt.box(True)
    plt.grid(True, linewidth=1)
    plt.xlabel('CO2 Purity [%]')
    plt.ylabel('CO2 Recovery [%]')
    plt.ylim([60, 100])
    plt.xlim([80, 100])
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.show()
