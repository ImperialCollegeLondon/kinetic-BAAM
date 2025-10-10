"""Simple smoke test for the translated kBAAM modules.

This script creates parameters, ensures DSL/LDF functions are attached,
and runs a short integration of the ODEs to confirm there are no runtime
errors. It prints the status and final state vector.
"""
from __future__ import annotations
import traceback
import numpy as np

from createParameters import create_parameters
from kBAAM_ODEs_nonIsothermal_ND import kbaam_odes_nonisothermal_nd


def main():
    try:
        params = create_parameters()

        # Ensure DSL and LDFCoefficient are present (create_parameters wires them
        # where possible). If still None, try importing directly.
        if params.get('DSL') is None:
            try:
                from DSL import DSL as DSL_impl
                params['DSL'] = DSL_impl
            except Exception:
                pass

        if params.get('LDFCoefficient') is None:
            try:
                from LDFCoefficient import LDFCoefficient as LDF_impl
                params['LDFCoefficient'] = LDF_impl
            except Exception:
                pass

        # Initial (dimensionless) state: [y1, q1, q2, T, Tw]
        y0 = np.array([params['y1_in'], 0.0, 0.0, 1.0, 1.0], dtype=float)

        # Short integration window (dimensionless)
        t_span = (0.0, 0.1)

        from scipy.integrate import solve_ivp

        sol = solve_ivp(lambda tt, xx: kbaam_odes_nonisothermal_nd(tt, xx, params, 'ads'), t_span, y0, rtol=1e-6, atol=1e-9)

        print('Integration success:', sol.success)
        print('Message:', sol.message)
        print('t_final =', sol.t[-1])
        print('y_final =', sol.y[:, -1])

    except Exception as e:
        print('Smoke test failed with exception:')
        traceback.print_exc()


if __name__ == '__main__':
    main()
