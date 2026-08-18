"""
k-BAAM Outputs GUI
==================
Tkinter GUI for running kBAAM_Outputs_nonIsothermal (Python port).

Launch:
    python kbaam_gui.py

Tabs:
  Inputs  — Basic (adsorbent / process type / theta / options)
          — Advanced (column geometry, heat transfer, feed)
  Outputs — Convergence (4-panel KPI vs cycle)
          — CSS Profiles (6-panel last-cycle state profiles)
"""
from __future__ import annotations
import os, sys, threading, traceback
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import scipy.io
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from createParameters import create_parameters
from DSL import DSL
from LDFCoefficient import LDFCoefficient
from kBAAM_Outputs_nonIsothermal import kbaam_outputs_nonisothermal

ADS_DIR = os.path.join(HERE, '..', 'AdsorbentFiles')
ADSORBENT_KEYS = [
    'qsb_1','qsd_1','bo_1','do_1','delUb_1','delUd_1',
    'qsb_2','qsd_2','bo_2','do_2','delUb_2','delUd_2',
    'rho_s','cp_g','cp_a','cp_w','cp_s','rho_w',
    'rp','V_column','t_wall','e_bed','h_in','h_out',
    'epsilon_p','Dm','tau','y1_in','T_feed',
]

# (label, param key, default, unit)
ADVANCED_FIELDS = [
    ('Column volume',        'V_column', 0.066,  'm\u00b3'),
    ('Particle radius',      'rp',       1e-3,   'm'),
    ('L / r  ratio',         'Lbyr',     7.0,    '\u2014'),
    ('Bed void fraction',    'e_bed',    0.37,   '\u2014'),
    ('Feed CO\u2082 mol frac', 'y1_in', 0.15,   '\u2014'),
    ('Feed temperature',     'T_feed',   298.0,  'K'),
    ('Inner HTC',            'h_in',     8.6,    'W/m\u00b2K'),
    ('Outer HTC',            'h_out',    2.5,    'W/m\u00b2K'),
    ('Solid density',        'rho_s',    1130.0, 'kg/m\u00b3'),
    ('Solid heat capacity',  'cp_s',     1070.0, 'J/kgK'),
    ('Wall heat capacity',   'cp_w',     502.0,  'J/kgK'),
    ('Wall density',         'rho_w',    7800.0, 'kg/m\u00b3'),
]

THETA_FIELDS = {
    'PVSA': [
        ('v_in',   0.111,  'm/s'),
        ('p_I',    0.4e5,  'Pa'),
        ('t_ads',  110.0,  's'),
        ('t_blo',  40.0,   's'),
        ('t_evac', 100.0,  's'),
        ('p_L',    0.02e5, 'Pa'),
        ('p_H',    3e5,    'Pa'),
    ],
    'VSA': [
        ('v_in',   0.3,    'm/s'),
        ('p_I',    0.3e5,  'Pa'),
        ('t_ads',  200.0,  's'),
        ('t_blo',  80.0,   's'),
        ('t_evac', 200.0,  's'),
        ('p_L',    0.02e5, 'Pa'),
    ],
}

# (title, unit, colour)
STATE_META = [
    ('CO\u2082 mole fraction  y\u2081', '\u2014',     '#1f77b4'),
    ('CO\u2082 loading  q\u2081',       'mol/kg',   '#ff7f0e'),
    ('N\u2082 loading  q\u2082',        'mol/kg',   '#2ca02c'),
    ('Temperature  T',                  'K',        '#d62728'),
    ('Wall temperature  T\u1d64',       'K',        '#9467bd'),
    ('Pressure  P',                     'Pa',       '#8c564b'),
]


def _load_adsorbent(mat_name: str) -> dict:
    path = os.path.join(ADS_DIR, mat_name)
    d = scipy.io.loadmat(path, simplify_cells=True)
    mat_p = d.get('parameters', {})
    params = create_parameters()
    for k in ADSORBENT_KEYS:
        if k in mat_p and isinstance(mat_p[k], (int, float)):
            params[k] = float(mat_p[k])
    params['DSL'] = DSL
    params['LDFCoefficient'] = LDFCoefficient
    return params


def _list_adsorbents() -> list:
    return sorted(f for f in os.listdir(ADS_DIR)
                  if f.endswith('.mat') and not f.startswith('.'))


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('k-BAAM Outputs')
        self.resizable(True, True)
        self._sim_thread = None
        self._profile_store: dict = {}
        self._build_ui()
        self._on_pt_change()

    # ── top-level layout ──────────────────────────────────────────────────────
    def _build_ui(self):
        self.columnconfigure(0, weight=0, minsize=285)
        self.columnconfigure(1, weight=1, minsize=560)
        self.rowconfigure(0, weight=1)
        left  = ttk.Frame(self, padding=6)
        right = ttk.Frame(self, padding=6)
        left.grid(row=0, column=0, sticky='nsew')
        right.grid(row=0, column=1, sticky='nsew')
        right.rowconfigure(1, weight=1)
        right.columnconfigure(0, weight=1)
        self._build_left(left)
        self._build_right(right)

    # ── LEFT: input notebook + run button ────────────────────────────────────
    def _build_left(self, parent):
        parent.rowconfigure(0, weight=1)
        parent.columnconfigure(0, weight=1)

        nb = ttk.Notebook(parent)
        nb.grid(row=0, column=0, sticky='nsew')

        basic_f  = ttk.Frame(nb, padding=8)
        adv_f    = ttk.Frame(nb, padding=8)
        create_f = ttk.Frame(nb, padding=0)
        opt_f    = ttk.Frame(nb, padding=8)
        nb.add(basic_f,  text='  Basic  ')
        nb.add(adv_f,    text='  Advanced  ')
        nb.add(create_f, text='  New Sorbent  ')
        nb.add(opt_f,    text='  Optimise  ')

        self._build_basic(basic_f)
        self._build_advanced(adv_f)
        self._build_create(create_f)
        self._build_opt_input(opt_f)

        btn_f = ttk.Frame(parent, padding=(6, 4))
        btn_f.grid(row=1, column=0, sticky='ew')
        btn_f.columnconfigure(0, weight=1)
        self._run_btn = ttk.Button(btn_f, text='\u25b6  Run', command=self._on_run)
        self._run_btn.grid(row=0, column=0, sticky='ew', ipady=5)
        self._status_var = tk.StringVar(value='Ready')
        ttk.Label(btn_f, textvariable=self._status_var,
                  foreground='gray').grid(row=1, column=0, sticky='w', pady=(3, 0))

    def _build_basic(self, f):
        f.columnconfigure(1, weight=1)
        row = 0

        ttk.Label(f, text='Adsorbent', font=('', 9, 'bold')).grid(
            row=row, column=0, columnspan=2, sticky='w'); row += 1
        self._ads_var = tk.StringVar()
        ads_cb = ttk.Combobox(f, textvariable=self._ads_var, state='readonly', width=30)
        ads_cb['values'] = _list_adsorbents()
        if '13XH_T.mat' in ads_cb['values']:
            ads_cb.current(ads_cb['values'].index('13XH_T.mat'))
        elif ads_cb['values']:
            ads_cb.current(0)
        ads_cb.grid(row=row, column=0, columnspan=2, sticky='ew', pady=(0, 6)); row += 1

        ttk.Label(f, text='Process type', font=('', 9, 'bold')).grid(
            row=row, column=0, columnspan=2, sticky='w'); row += 1
        self._pt_var = tk.StringVar(value='PVSA')
        pt_cb = ttk.Combobox(f, textvariable=self._pt_var, state='readonly',
                              width=30, values=['PVSA', 'VSA'])
        pt_cb.grid(row=row, column=0, columnspan=2, sticky='ew', pady=(0, 6)); row += 1
        pt_cb.bind('<<ComboboxSelected>>', lambda _: self._on_pt_change())

        ttk.Label(f, text='\u03b8  Parameters', font=('', 9, 'bold')).grid(
            row=row, column=0, columnspan=2, sticky='w'); row += 1
        self._theta_frame = ttk.Frame(f)
        self._theta_frame.grid(row=row, column=0, columnspan=2, sticky='ew',
                               pady=(0, 6)); row += 1
        self._theta_frame.columnconfigure(1, weight=1)
        self._theta_entries = {}

        ttk.Label(f, text='Options', font=('', 9, 'bold')).grid(
            row=row, column=0, columnspan=2, sticky='w'); row += 1
        self._opt_pd      = tk.BooleanVar(value=True)
        self._opt_heating = tk.BooleanVar(value=False)
        self._opt_cCSTR   = tk.BooleanVar(value=False)
        for lbl, var in [('Pressure drop', self._opt_pd),
                         ('Heating (TVSA)', self._opt_heating),
                         ('CSTR mixing',    self._opt_cCSTR)]:
            ttk.Checkbutton(f, text=lbl, variable=var).grid(
                row=row, column=0, columnspan=2, sticky='w'); row += 1

    def _build_advanced(self, f):
        f.columnconfigure(2, weight=1)
        ttk.Label(f, text='Column geometry & physical properties',
                  font=('', 9, 'bold')).grid(
            row=0, column=0, columnspan=3, sticky='w', pady=(0, 6))
        self._adv_vars = {}
        for r, (lbl, key, default, unit) in enumerate(ADVANCED_FIELDS, start=1):
            ttk.Label(f, text=lbl, width=22, anchor='w').grid(
                row=r, column=0, sticky='w', pady=1)
            var = tk.StringVar(value=str(default))
            ttk.Entry(f, textvariable=var, width=12).grid(
                row=r, column=1, padx=4, pady=1, sticky='ew')
            ttk.Label(f, text=unit, foreground='gray', width=8).grid(
                row=r, column=2, sticky='w')
            self._adv_vars[key] = var

    def _build_create(self, f):
        """Scrollable form to define and save a new adsorbent .mat file."""
        f.rowconfigure(0, weight=1); f.columnconfigure(0, weight=1)

        canvas = tk.Canvas(f, borderwidth=0, highlightthickness=0)
        sb = ttk.Scrollbar(f, orient='vertical', command=canvas.yview)
        canvas.configure(yscrollcommand=sb.set)
        canvas.grid(row=0, column=0, sticky='nsew')
        sb.grid(row=0, column=1, sticky='ns')

        inner = ttk.Frame(canvas, padding=8)
        win_id = canvas.create_window((0, 0), window=inner, anchor='nw')
        inner.bind('<Configure>',
                   lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.bind('<Configure>',
                    lambda e: canvas.itemconfig(win_id, width=e.width))
        # mouse-wheel scrolling
        def _scroll(e): canvas.yview_scroll(int(-1*(e.delta/120)), 'units')
        canvas.bind_all('<MouseWheel>', _scroll)

        inner.columnconfigure(2, weight=1)
        r = 0

        def _section(label):
            nonlocal r
            if r: ttk.Separator(inner, orient='horizontal').grid(
                row=r, column=0, columnspan=3, sticky='ew', pady=(6, 2)); r += 1
            ttk.Label(inner, text=label, font=('', 9, 'bold')).grid(
                row=r, column=0, columnspan=3, sticky='w', pady=(0, 3)); r += 1

        def _row(label, key, default, unit):
            nonlocal r
            ttk.Label(inner, text=label, width=20, anchor='w').grid(
                row=r, column=0, sticky='w', pady=1)
            var = tk.StringVar(value=str(default))
            ttk.Entry(inner, textvariable=var, width=12).grid(
                row=r, column=1, padx=4, pady=1, sticky='ew')
            ttk.Label(inner, text=unit, foreground='gray').grid(
                row=r, column=2, sticky='w')
            self._create_vars[key] = var; r += 1

        self._create_vars: dict[str, tk.StringVar] = {}
        self._create_name_var = tk.StringVar(value='NewSorbent')

        ttk.Label(inner, text='File name', font=('', 9, 'bold')).grid(
            row=r, column=0, columnspan=3, sticky='w'); r += 1
        ttk.Entry(inner, textvariable=self._create_name_var, width=28).grid(
            row=r, column=0, columnspan=2, sticky='ew', pady=(0, 4)); r += 1
        ttk.Label(inner, text='.mat', foreground='gray').grid(
            row=r-1, column=2, sticky='w')

        _section('Component 1 (CO\u2082) — DSL isotherm')
        _row('q_sb,1  (b-site saturation)', 'qsb_1',   3.0,      'mol/kg')
        _row('q_sd,1  (d-site saturation)', 'qsd_1',   3.0,      'mol/kg')
        _row('b\u2080,1',                   'bo_1',    3e-7,     'm\u00b3/mol')
        _row('d\u2080,1',                   'do_1',    1e-6,     'm\u00b3/mol')
        _row('\u0394U_b,1',                 'delUb_1', -30000.0, 'J/mol')
        _row('\u0394U_d,1',                 'delUd_1', -35000.0, 'J/mol')

        _section('Component 2 (N\u2082) — DSL isotherm')
        _row('q_sb,2',                      'qsb_2',   3.0,      'mol/kg')
        _row('q_sd,2',                      'qsd_2',   3.0,      'mol/kg')
        _row('b\u2080,2',                   'bo_2',    8e-7,     'm\u00b3/mol')
        _row('d\u2080,2',                   'do_2',    5e-6,     'm\u00b3/mol')
        _row('\u0394U_b,2',                 'delUb_2', -13000.0, 'J/mol')
        _row('\u0394U_d,2',                 'delUd_2', -16000.0, 'J/mol')

        _section('Physical properties')
        _row('Solid density',               'rho_s',    1130.0,  'kg/m\u00b3')
        _row('Solid heat capacity',         'cp_s',     1070.0,  'J/kgK')
        _row('Particle radius',             'rp',       1e-3,    'm')
        _row('Particle porosity',           'epsilon_p', 0.35,   '\u2014')
        _row('Molecular diffusivity',       'Dm',       1.6e-5,  'm\u00b2/s')
        _row('Tortuosity',                  'tau',      3.0,     '\u2014')

        _section('Bed & wall properties')
        _row('Bed void fraction',           'e_bed',    0.37,    '\u2014')
        _row('Column volume',               'V_column', 0.066,   'm\u00b3')
        _row('Inner HTC',                   'h_in',     8.6,     'W/m\u00b2K')
        _row('Outer HTC',                   'h_out',    2.5,     'W/m\u00b2K')
        _row('Wall heat capacity',          'cp_w',     502.0,   'J/kgK')
        _row('Wall density',                'rho_w',    7800.0,  'kg/m\u00b3')

        _section('Feed conditions')
        _row('Feed CO\u2082 mole fraction', 'y1_in',    0.15,    '\u2014')
        _row('Feed temperature',            'T_feed',   298.0,   'K')

        r += 1
        save_btn = ttk.Button(inner, text='\U0001f4be  Save to AdsorbentFiles/',
                              command=self._on_save_adsorbent)
        save_btn.grid(row=r, column=0, columnspan=3, sticky='ew',
                      ipady=4, pady=(6, 2)); r += 1
        self._create_status_var = tk.StringVar(value='')
        ttk.Label(inner, textvariable=self._create_status_var,
                  foreground='gray').grid(row=r, column=0, columnspan=3, sticky='w')

    def _on_save_adsorbent(self):
        import math
        name = self._create_name_var.get().strip()
        if not name:
            messagebox.showerror('Save error', 'Please enter a file name.')
            return
        params = {}
        for key, var in self._create_vars.items():
            try:
                params[key] = float(var.get())
            except ValueError:
                messagebox.showerror('Save error', f'Invalid value for {key}.')
                return
        # also store name and scalar default fields expected by the outputs code
        params['adsorbentName'] = name
        params['cp_g']   = 30.7
        params['cp_a']   = 30.7
        params['t_wall'] = 0.003
        params['Lbyr']   = 7.0
        import numpy as np
        mat_params = {k: np.array([[v]]) if isinstance(v, float) else v
                      for k, v in params.items()}
        out_path = os.path.join(ADS_DIR, f'{name}.mat')
        scipy.io.savemat(out_path, {'parameters': mat_params})
        # refresh adsorbent dropdown
        new_list = _list_adsorbents()
        self._ads_var.set(f'{name}.mat')
        self._create_status_var.set(f'Saved: {name}.mat')
        self.after(3000, lambda: self._create_status_var.set(''))

    def _on_pt_change(self):
        for w in self._theta_frame.winfo_children():
            w.destroy()
        self._theta_entries.clear()
        for r, (name, default, unit) in enumerate(THETA_FIELDS[self._pt_var.get()]):
            ttk.Label(self._theta_frame, text=f'{name}  [{unit}]',
                      width=18).grid(row=r, column=0, sticky='w', pady=1)
            var = tk.StringVar(value=str(default))
            ttk.Entry(self._theta_frame, textvariable=var, width=13).grid(
                row=r, column=1, sticky='ew', padx=(4, 0), pady=1)
            self._theta_entries[name] = var

    # ── RIGHT: KPI bar + output notebook ─────────────────────────────────────
    def _build_right(self, parent):
        kpi_frame = ttk.LabelFrame(parent, text='Cyclic Steady-State KPIs', padding=6)
        kpi_frame.grid(row=0, column=0, sticky='ew', pady=(0, 6))
        for i in range(4):
            kpi_frame.columnconfigure(i, weight=1)

        headers = ['CO\u2082 Purity', 'CO\u2082 Recovery', 'Productivity', 'SEC']
        units   = ['%', '%', 'mol/m\u00b3/s', 'kWh/t CO\u2082']
        self._kpi_vars = []
        for col, (h, u) in enumerate(zip(headers, units)):
            ttk.Label(kpi_frame, text=h, font=('', 8, 'bold'),
                      anchor='center').grid(row=0, column=col, padx=6)
            ttk.Label(kpi_frame, text=u, foreground='gray', font=('', 7),
                      anchor='center').grid(row=1, column=col, padx=6)
            v = tk.StringVar(value='\u2014')
            ttk.Label(kpi_frame, textvariable=v, font=('', 14),
                      anchor='center').grid(row=2, column=col, padx=6, pady=2)
            self._kpi_vars.append(v)
        self._cycles_var = tk.StringVar(value='')
        ttk.Label(kpi_frame, textvariable=self._cycles_var,
                  foreground='gray', font=('', 7)).grid(
            row=3, column=0, columnspan=4, sticky='e', padx=6)

        out_nb = ttk.Notebook(parent)
        out_nb.grid(row=1, column=0, sticky='nsew')
        parent.rowconfigure(1, weight=1)

        conv_f    = ttk.Frame(out_nb)
        profile_f = ttk.Frame(out_nb)
        pareto_f  = ttk.Frame(out_nb)
        out_nb.add(conv_f,    text='  Convergence  ')
        out_nb.add(profile_f, text='  CSS Profiles  ')
        out_nb.add(pareto_f,  text='  Pareto Front  ')

        self._build_conv_tab(conv_f)
        self._build_profile_tab(profile_f)
        self._build_pareto_tab(pareto_f)

        # ── log box ───────────────────────────────────────────────────────────
        log_hdr = ttk.Frame(parent)
        log_hdr.grid(row=2, column=0, sticky='ew', pady=(6, 0))
        log_hdr.columnconfigure(0, weight=1)
        ttk.Label(log_hdr, text='Console log', font=('', 8, 'bold')).grid(
            row=0, column=0, sticky='w')
        ttk.Button(log_hdr, text='Clear', width=6,
                   command=lambda: (self._log.configure(state='normal'),
                                    self._log.delete('1.0', 'end'),
                                    self._log.configure(state='disabled'))
                   ).grid(row=0, column=1, sticky='e')

        from tkinter.scrolledtext import ScrolledText
        self._log = ScrolledText(parent, height=5, state='disabled',
                                 font=('Courier', 8), wrap='word',
                                 background='#1e1e1e', foreground='#d4d4d4',
                                 insertbackground='white',
                                 selectbackground='#264f78')
        self._log.grid(row=3, column=0, sticky='ew', pady=(2, 0))

        # redirect stdout/stderr to the log widget
        import sys as _sys
        _orig_out, _orig_err = _sys.__stdout__, _sys.__stderr__
        _log_w = self._log

        class _Redir:
            def __init__(self, orig): self._orig = orig
            def write(self, s):
                if s:
                    self._orig.write(s)
                    try:
                        _log_w.after(0, _Redir._append, _log_w, s)
                    except Exception:
                        pass
            @staticmethod
            def _append(w, s):
                w.configure(state='normal')
                w.insert('end', s)
                w.see('end')
                w.configure(state='disabled')
            def flush(self): self._orig.flush()

        _sys.stdout = _Redir(_orig_out)
        _sys.stderr = _Redir(_orig_err)

        def _restore():
            _sys.stdout = _orig_out
            _sys.stderr = _orig_err
        self.protocol('WM_DELETE_WINDOW', lambda: (_restore(), self.destroy()))

    def _build_conv_tab(self, f):
        f.rowconfigure(0, weight=1); f.columnconfigure(0, weight=1)
        self._fig_conv, self._axes_conv = plt.subplots(2, 2, figsize=(6, 4.2))
        self._fig_conv.tight_layout(pad=2.2)
        canvas = FigureCanvasTkAgg(self._fig_conv, master=f)
        canvas.get_tk_widget().grid(row=0, column=0, sticky='nsew')
        canvas.draw()
        self._canvas_conv = canvas

    def _build_profile_tab(self, f):
        f.rowconfigure(0, weight=1); f.columnconfigure(0, weight=1)
        self._fig_prof, self._axes_prof = plt.subplots(3, 2, figsize=(6, 5.4))
        self._fig_prof.tight_layout(pad=2.4)
        canvas = FigureCanvasTkAgg(self._fig_prof, master=f)
        canvas.get_tk_widget().grid(row=0, column=0, sticky='nsew')
        canvas.draw()
        self._canvas_prof = canvas

    def _build_pareto_tab(self, f):
        f.rowconfigure(0, weight=1); f.columnconfigure(0, weight=1)
        self._fig_pareto, self._ax_pareto = plt.subplots(1, 1, figsize=(5.5, 4.5))
        self._fig_pareto.tight_layout(pad=2.5)
        canvas = FigureCanvasTkAgg(self._fig_pareto, master=f)
        canvas.get_tk_widget().grid(row=0, column=0, sticky='nsew')
        canvas.draw()
        self._canvas_pareto = canvas

    def _build_opt_input(self, f):
        f.columnconfigure(1, weight=1)
        r = 0

        ttk.Label(f, text='Optimisation settings', font=('', 9, 'bold')).grid(
            row=r, column=0, columnspan=3, sticky='w', pady=(0, 6)); r += 1

        def _row(lbl, var, width=10):
            ttk.Label(f, text=lbl, width=18, anchor='w').grid(
                row=r, column=0, sticky='w', pady=2)
            ttk.Entry(f, textvariable=var, width=width).grid(
                row=r, column=1, sticky='ew', padx=(4, 0), pady=2)

        ttk.Label(f, text='Process type', anchor='w').grid(
            row=r, column=0, sticky='w', pady=2)
        self._opt_pt_var = tk.StringVar(value='PVSA')
        ttk.Combobox(f, textvariable=self._opt_pt_var, state='readonly', width=18,
                     values=['PVSA', 'VSA', 'AdsorbentPVSA', 'AdsorbentVSA',
                             'AdsorbentPVSADSL', 'AdsorbentVSAb0',
                             'AdsorbentPVSAb0', 'Resin']).grid(
            row=r, column=1, sticky='ew', padx=(4, 0), pady=2); r += 1

        self._opt_ngens_var   = tk.StringVar(value='90')
        self._opt_popsize_var = tk.StringVar(value='200')
        self._opt_ncores_var  = tk.StringVar(value='4')

        _row('Generations',  self._opt_ngens_var);   r += 1
        _row('Pop. size',    self._opt_popsize_var);  r += 1
        _row('CPU cores',    self._opt_ncores_var);   r += 1

        ttk.Separator(f, orient='horizontal').grid(
            row=r, column=0, columnspan=2, sticky='ew', pady=6); r += 1

        ttk.Label(f, text='Objective type', anchor='w').grid(
            row=r, column=0, sticky='w', pady=2)
        self._opt_OptType_var = tk.StringVar(value='Const')
        ttk.Combobox(f, textvariable=self._opt_OptType_var, state='readonly',
                     values=['Const', 'Unc'], width=12).grid(
            row=r, column=1, sticky='ew', padx=(4, 0), pady=2); r += 1

        ttk.Label(f, text='Press. type', anchor='w').grid(
            row=r, column=0, sticky='w', pady=2)
        self._opt_pressType_var = tk.StringVar(value='FP')
        ttk.Combobox(f, textvariable=self._opt_pressType_var, state='readonly',
                     values=['FP', 'LPP'], width=12).grid(
            row=r, column=1, sticky='ew', padx=(4, 0), pady=2); r += 1

        # p_L is a fixed parameter for PVSA and a decision variable for VSA
        ttk.Label(f, text='p_L  [Pa]', anchor='w').grid(
            row=r, column=0, sticky='w', pady=2)
        self._opt_pL_var = tk.StringVar(value='2000.0')
        ttk.Entry(f, textvariable=self._opt_pL_var, width=12).grid(
            row=r, column=1, sticky='ew', padx=(4, 0), pady=2)
        ttk.Label(f, text='fixed / lb for VSA', foreground='gray',
                  font=('', 7)).grid(row=r, column=2, sticky='w', padx=4); r += 1

        ttk.Separator(f, orient='horizontal').grid(
            row=r, column=0, columnspan=2, sticky='ew', pady=6); r += 1

        # run / stop buttons
        self._opt_run_btn  = ttk.Button(f, text='\u25b6  Run Optimisation',
                                        command=self._on_run_opt)
        self._opt_stop_btn = ttk.Button(f, text='\u25a0  Stop',
                                        command=self._on_stop_opt,
                                        state='disabled')
        self._opt_run_btn.grid( row=r, column=0, columnspan=2, sticky='ew',
                                ipady=4, pady=(0, 3)); r += 1
        self._opt_stop_btn.grid(row=r, column=0, columnspan=2, sticky='ew',
                                ipady=2); r += 1

        # progress
        self._opt_prog_var   = tk.StringVar(value='')
        self._opt_status_var = tk.StringVar(value='')
        ttk.Label(f, textvariable=self._opt_prog_var,
                  font=('', 8)).grid(row=r+1, column=0, columnspan=2,
                                     sticky='w', pady=(6, 0))
        ttk.Label(f, textvariable=self._opt_status_var,
                  foreground='gray', font=('', 7)).grid(
            row=r+2, column=0, columnspan=2, sticky='w')

    # ── Optimisation run / stop ───────────────────────────────────────────────
    def _on_run_opt(self):
        if getattr(self, '_opt_thread', None) and self._opt_thread.is_alive():
            return
        self._opt_stop_event = threading.Event()
        self._opt_store: dict = {}
        self._opt_run_btn.state(['disabled'])
        self._opt_stop_btn.state(['!disabled'])
        self._opt_prog_var.set('Starting …')
        self._opt_status_var.set('')
        self._ax_pareto.clear()
        self._canvas_pareto.draw()
        self._opt_thread = threading.Thread(target=self._run_opt_thread, daemon=True)
        self._opt_thread.start()
        self.after(800, self._poll_opt_progress)

    def _on_stop_opt(self):
        if hasattr(self, '_opt_stop_event'):
            self._opt_stop_event.set()
        self._opt_status_var.set('Stopping …')

    def _run_opt_thread(self):
        try:
            import numpy as np
            from run_NSGA import run_nsga
            from pymoo.core.callback import Callback

            params = _load_adsorbent(self._ads_var.get())
            for key, var in self._adv_vars.items():
                try: params[key] = float(var.get())
                except ValueError: pass

            pt = self._opt_pt_var.get()
            params.update(dict(
                outputType='opt', processType=pt,
                OptType=self._opt_OptType_var.get(),
                pressType=self._opt_pressType_var.get(),
                modelType='nonisothermal',
                Lbyr=float(self._adv_vars['Lbyr'].get()),
                p_L=float(self._opt_pL_var.get()),
                pressureDrop=self._opt_pd.get(),
                heating=self._opt_heating.get(),
                cCSTR=self._opt_cCSTR.get(),
                SSLSTA=False, equilibrium=False,
                amine=False, testBT=False, testEvac=False,
                normPlot=False, forwardEvac=False,
                t_press=20.0, rigid=True, plot0D=False,
                n_cores=int(self._opt_ncores_var.get()),
                DSL=DSL, LDFCoefficient=LDFCoefficient,
            ))

            store = self._opt_store
            stop_ev = self._opt_stop_event

            class _Cb(Callback):
                def notify(self, algo):
                    store['gen']    = algo.n_gen
                    store['n_eval'] = algo.evaluator.n_eval
                    if algo.pop is not None:
                        store['F_pop'] = algo.pop.get('F')
                    if stop_ev.is_set():
                        algo.termination.force_termination = True

            import time
            t0 = time.time()
            ngens    = int(self._opt_ngens_var.get())
            pop_size = int(self._opt_popsize_var.get())
            X, F = run_nsga(params, ngens=ngens, pop_size=pop_size,
                            use_pymoo=True,
                            n_cores=int(self._opt_ncores_var.get()))
            elapsed = time.time() - t0
            store['X_final'] = X
            store['F_final'] = F
            store['elapsed'] = elapsed
            store['done'] = True

        except Exception:
            self._opt_store['error'] = traceback.format_exc()
            self._opt_store['done'] = True

    def _poll_opt_progress(self):
        store = self._opt_store
        if store.get('error'):
            self._opt_prog_var.set('Error — see details')
            self._opt_status_var.set(store['error'][:120])
            self._opt_run_btn.state(['!disabled'])
            self._opt_stop_btn.state(['disabled'])
            return
        gen    = store.get('gen', 0)
        n_eval = store.get('n_eval', 0)
        ngens  = int(self._opt_ngens_var.get())
        self._opt_prog_var.set(
            f'Generation {gen} / {ngens}  |  Evaluations: {n_eval}')
        if store.get('F_pop') is not None:
            self._plot_pareto_live(store['F_pop'], gen, ngens)
        if store.get('done'):
            self._opt_run_btn.state(['!disabled'])
            self._opt_stop_btn.state(['disabled'])
            F = store.get('F_final')
            elapsed = store.get('elapsed', 0)
            if F is not None:
                self._opt_status_var.set(
                    f'Done — {len(F)} Pareto points  |  {elapsed:.0f} s')
                self._plot_pareto_final(F, store.get('X_final'))
            return
        self.after(800, self._poll_opt_progress)

    def _plot_pareto_live(self, F, gen, ngens):
        ax = self._ax_pareto
        ax.clear()
        ax.scatter(F[:, 0], F[:, 1], c='#1f77b4', s=18, alpha=0.7,
                   edgecolors='none')
        ax.set_xlabel('Objective 1', fontsize=8)
        ax.set_ylabel('Objective 2', fontsize=8)
        ax.set_title(f'Population — gen {gen}/{ngens}', fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.25)
        self._fig_pareto.tight_layout()
        self._canvas_pareto.draw()

    def _plot_pareto_final(self, F, X):
        opt_type = self._opt_OptType_var.get()
        ax = self._ax_pareto
        ax.clear()

        # label axes based on OptType
        if opt_type == 'Const':
            xlabel = 'Productivity penalty  (−productivity + φ)'
            ylabel = 'SEC penalty  (SEC·conv + φ)'
        else:
            xlabel = '−Recovery [%]'
            ylabel = '−Purity [%]'

        ax.scatter(F[:, 0], F[:, 1], c='#d62728', s=28, zorder=3,
                   edgecolors='none', label='Pareto front')
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ads = self._ads_var.get().replace('.mat', '')
        ax.set_title(
            f'Pareto front — {ads} {self._pt_var.get()}\n'
            f'{len(F)} points  |  OptType: {opt_type}',
            fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.25)
        self._fig_pareto.tight_layout()
        self._canvas_pareto.draw()

    # ── Run ───────────────────────────────────────────────────────────────────
    def _on_run(self):
        if self._sim_thread and self._sim_thread.is_alive():
            return
        self._run_btn.state(['disabled'])
        self._status_var.set('Running \u2026')
        for v in self._kpi_vars: v.set('\u2026')
        self._cycles_var.set('')
        self._profile_store.clear()
        self._sim_thread = threading.Thread(target=self._run_sim, daemon=True)
        self._sim_thread.start()

    def _run_sim(self):
        try:
            params = _load_adsorbent(self._ads_var.get())
            # apply advanced overrides
            for key, var in self._adv_vars.items():
                try: params[key] = float(var.get())
                except ValueError: pass
            pt = self._pt_var.get()
            params.update(dict(
                outputType='plot', processType=pt, OptType='Unc',
                Lbyr=float(self._adv_vars['Lbyr'].get()),
                t_press=20.0,
                pressureDrop=self._opt_pd.get(),
                heating=self._opt_heating.get(),
                cCSTR=self._opt_cCSTR.get(),
                SSLSTA=False, equilibrium=False,
                amine=False, testBT=False, testEvac=False,
                normPlot=False, forwardEvac=False,
            ))
            theta = [float(v.get()) for v in self._theta_entries.values()]
            import time
            t0 = time.time()
            KPIs_raw = kbaam_outputs_nonisothermal(
                params, thetaIn=theta, raise_on_error=True,
                solver_method='BDF', _profile_store=self._profile_store)
            elapsed = time.time() - t0

            KPIs = np.asarray(KPIs_raw)
            if KPIs.ndim == 2 and KPIs.shape[0] > 0:
                final = KPIs[-1]
                purity, recovery, prod, sec = final
                n_cycles = KPIs.shape[0]
            else:
                recovery, purity = -float(KPIs[0]), -float(KPIs[1])
                prod = sec = float('nan'); n_cycles = 1; KPIs = None

            self.after(0, self._update_results,
                       purity, recovery, prod, sec, n_cycles, KPIs, elapsed, theta, pt)
        except Exception:
            self.after(0, self._on_error, traceback.format_exc())

    # ── Result update ─────────────────────────────────────────────────────────
    def _update_results(self, purity, recovery, prod, sec,
                        n_cycles, KPIs, elapsed, theta, pt):
        self._kpi_vars[0].set(f'{purity:.2f}')
        self._kpi_vars[1].set(f'{recovery:.2f}')
        self._kpi_vars[2].set(f'{prod:.4g}')
        self._kpi_vars[3].set(f'{sec:.4g}')
        self._cycles_var.set(f'CSS in {n_cycles} cycles  |  {elapsed:.1f} s')
        self._status_var.set(f'Done  ({elapsed:.1f} s)')
        self._run_btn.state(['!disabled'])
        ads = self._ads_var.get().replace('.mat', '')
        theta_str = ', '.join(f'{v:.3g}' for v in theta)
        suptitle = f'{ads} \u2014 {pt}\n\u03b8 = [{theta_str}]'
        if KPIs is not None and KPIs.ndim == 2 and KPIs.shape[0] > 1:
            self._plot_convergence(KPIs, [purity, recovery, prod, sec], suptitle)
        if self._profile_store.get('X_cycle') is not None:
            self._plot_profiles(suptitle)

    def _plot_convergence(self, KPIs, css_vals, suptitle):
        cycles  = np.arange(1, KPIs.shape[0] + 1)
        titles  = ['CO\u2082 Purity [%]', 'CO\u2082 Recovery [%]',
                   'Productivity [mol/m\u00b3/s]', 'SEC [kWh/t CO\u2082]']
        colors  = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
        for ax, col, ci, title, css in zip(
                self._axes_conv.flat, colors, range(4), titles, css_vals):
            ax.clear()
            ax.plot(cycles, KPIs[:, ci], color=col, lw=1.6)
            ax.axhline(css, color=col, ls='--', lw=0.8, alpha=0.5)
            ax.set_title(f'{title}\nCSS: {css:.3g}', fontsize=7.5)
            ax.set_xlabel('Cycle', fontsize=7)
            ax.tick_params(labelsize=7)
            ax.grid(True, alpha=0.25)
        self._fig_conv.suptitle(suptitle, fontsize=8)
        self._fig_conv.tight_layout(rect=[0, 0, 1, 0.90])
        self._canvas_conv.draw()

    def _plot_profiles(self, suptitle):
        t   = self._profile_store['t_cycle']
        X   = self._profile_store['X_cycle']
        t_s = self._profile_store['t_steps']   # (t_ads_end, t_blo_end, t_evac_end)
        step_labels = ['ads', 'blo', 'evac', 'pres']
        step_colors = ['#4daf4a', '#ff7f00', '#377eb8', '#984ea3']
        bounds = [0.0, t_s[0], t_s[1], t_s[2], t[-1]]
        for ax, (title, unit, col), ci in zip(
                self._axes_prof.flat, STATE_META, range(6)):
            ax.clear()
            ax.plot(t, X[:, ci], color=col, lw=1.4)
            for i, (ta, tb) in enumerate(zip(bounds, bounds[1:])):
                ax.axvspan(ta, tb, alpha=0.07, color=step_colors[i])
            for tx in t_s:
                ax.axvline(tx, color='gray', lw=0.6, ls='--')
            ax.set_title(title, fontsize=7.5)
            ax.set_ylabel(unit, fontsize=7)
            ax.set_xlabel('time [s]', fontsize=7)
            ax.tick_params(labelsize=7)
            ax.grid(True, alpha=0.2)
            # round T and Tw y-limits outward to nearest 5 K
            if ci in (3, 4):
                import math
                ylo = math.floor(X[:, ci].min() / 5) * 5
                yhi = math.ceil(X[:, ci].max()  / 5) * 5
                ax.set_ylim(ylo, yhi)
        from matplotlib.patches import Patch
        patches = [Patch(color=c, alpha=0.25, label=l)
                   for c, l in zip(step_colors, step_labels)]
        self._axes_prof.flat[0].legend(handles=patches, fontsize=6,
                                       loc='upper right', framealpha=0.7)
        self._fig_prof.suptitle(suptitle, fontsize=8)
        self._fig_prof.tight_layout(rect=[0, 0, 1, 0.92])
        self._canvas_prof.draw()

    def _on_error(self, msg):
        self._status_var.set('Error \u2014 see console')
        self._run_btn.state(['!disabled'])
        for v in self._kpi_vars: v.set('\u2014')
        messagebox.showerror('Simulation error', msg[:1600])


if __name__ == '__main__':
    App().mainloop()
