"""Smoke test mirroring kBAAM_example.m exactly.

theta = [v_in=0.3*0.37, p_I=0.4e5, t_ads=110, t_blo=40, t_evac=100, p_L=0.02e5, p_H=3e5]
Adsorbent: 13XH_T  (PVSA, pressureDrop=True, non-isothermal)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import scipy.io
from createParameters import create_parameters
from DSL import DSL
from LDFCoefficient import LDFCoefficient
from kBAAM_Outputs_nonIsothermal import kbaam_outputs_nonisothermal

# Load 13XH_T adsorbent (mirrors MATLAB load('13XH_T.mat'))
mat = scipy.io.loadmat(
    os.path.join(HERE, '..', 'AdsorbentFiles', '13XH_T.mat'),
    simplify_cells=True,
)
mat_p = mat['parameters']

parameters = create_parameters()

for k in ['qsb_1','qsd_1','bo_1','do_1','delUb_1','delUd_1',
          'qsb_2','qsd_2','bo_2','do_2','delUb_2','delUd_2',
          'rho_s','cp_g','cp_a','cp_w','cp_s','rho_w',
          'rp','V_column','t_wall','e_bed','h_in','h_out',
          'epsilon_p','Dm','tau','y1_in','T_feed']:
    if k in mat_p and isinstance(mat_p[k], (int, float)):
        parameters[k] = float(mat_p[k])

parameters.update(dict(
    outputType='plot', processType='PVSA', OptType='Unc',
    Lbyr=7.0, equilibrium=False, cCSTR=False,
    testBT=False, testEvac=False, normPlot=False,
    forwardEvac=False, pressureDrop=True, SSLSTA=False,
    t_press=20.0, heating=False, amine=False,
    DSL=DSL, LDFCoefficient=LDFCoefficient,
))

theta = [0.3 * 0.37, 0.4e5, 110.0, 40.0, 100.0, 0.02e5, 3e5]

print("Running kBAAM_example smoke test (13XH_T, PVSA, pressureDrop=True) ...")
print(f"  theta = {theta}")

import time
t0 = time.time()
KPIs_raw = kbaam_outputs_nonisothermal(parameters, thetaIn=theta,
                                        raise_on_error=True, solver_method='BDF')
elapsed = time.time() - t0
print(f"\nSimulation time: {elapsed:.1f} s")

# outputType='plot' returns process_indicators.T  shape (n_cycles, 4)
# columns: [purity, recovery, productivity, SEC]
KPIs = np.asarray(KPIs_raw)
if KPIs.ndim == 2:
    final        = KPIs[-1]
    purity       = final[0]
    recovery     = final[1]
    productivity = final[2]
    SEC          = final[3]
    n_cycles     = KPIs.shape[0]
else:
    recovery, purity = -float(KPIs[0]), -float(KPIs[1])
    productivity = SEC = float('nan')
    n_cycles = 1

print("\n" + "="*52)
print("  FINAL CYCLIC STEADY-STATE KPIs")
print("="*52)
print(f"  CO2 Purity      : {purity:.2f} %")
print(f"  CO2 Recovery    : {recovery:.2f} %")
print(f"  Productivity    : {productivity:.4g} mol/m\u00b3/s")
print(f"  SEC             : {SEC:.4g} kWh/tonne CO2")
print(f"  CSS reached in  : {n_cycles} cycles")
print("="*52)

# ── build convergence plot ────────────────────────────────────────────────────
out_path = os.path.join(HERE, 'smoke_kbaam_example_output.png')

if KPIs.ndim == 2 and KPIs.shape[0] > 1:
    cycles = np.arange(1, n_cycles + 1)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle(
        'k-BAAM: 13XH_T PVSA — Cyclic Convergence\n'
        r'$\theta$=[v_in=0.111 m/s, p_I=0.4 bar, t_ads=110 s, t_blo=40 s, '
        r't_evac=100 s, p_L=0.02 bar, p_H=3 bar]',
        fontsize=10,
    )

    labels   = ['CO\u2082 Purity [%]', 'CO\u2082 Recovery [%]',
                'Productivity [mol/m\u00b3/s]', 'SEC [kWh/tonne CO\u2082]']
    colors   = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    col_idx  = [0, 1, 2, 3]

    for ax, lbl, col, ci in zip(axes.flat, labels, colors, col_idx):
        ax.plot(cycles, KPIs[:, ci], color=col, lw=1.8)
        # mark the final (CSS) value
        ax.axhline(KPIs[-1, ci], color=col, ls='--', lw=0.8, alpha=0.6)
        ax.set_xlabel('Cycle')
        ax.set_ylabel(lbl)
        ax.set_title(f'CSS: {KPIs[-1, ci]:.3g}')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f'\nConvergence plot saved to: {out_path}')
    plt.close(fig)
else:
    print('(only one cycle — no convergence plot generated)')
