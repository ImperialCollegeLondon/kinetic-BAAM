from __future__ import annotations
"""Quick example plots for kinetic-BAAM Python port.

This script mirrors the MATLAB `kBAAM_example.m` at a high level by
showing how to generate common plots from saved results.

It expects text result files written by the Outputs driver to `rawData/`.
Adjust the `files` lists below to point at your actual runs.
"""
import os
import numpy as np

from .plotPurityRecoveryOpt import plot_purity_recovery_from_txt
from .plotEnProdOpt import plot_energy_vs_productivity
from .plotHypotheticalIsotherms import plot_hypothetical_isotherms


def main():
    raw = 'rawData'
    files = []
    if os.path.isdir(raw):
        # pick a few newest .txt logs if present
        txts = [os.path.join(raw, f) for f in os.listdir(raw) if f.endswith('.txt')]
        files = sorted(txts, key=os.path.getmtime, reverse=True)[:2]

    if files:
        print('Plotting Purity vs Recovery for:', files)
        plot_purity_recovery_from_txt(files)
        print('Plotting Energy vs Productivity for:', files)
        plot_energy_vs_productivity(files)
    else:
        print('No rawData/*.txt files found; skipping KPI plots.')

    # Hypothetical isotherms demo (requires DSL shape-consistent theta)
    # Placeholder theta (indices 4..7 used; adapt based on your optimisation encoding)
    theta = [0, 0, 0, 0, 3.0, 1e-6, 2e4, 1e4]
    Pvals = np.linspace(0, 1e5, 200)
    print('Plotting hypothetical isotherms...')
    plot_hypothetical_isotherms(Pvals, theta)


if __name__ == '__main__':
    main()
