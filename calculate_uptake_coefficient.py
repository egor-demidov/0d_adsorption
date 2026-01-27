from pathlib import Path

import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np
import json


DATASET_PATH = Path('uptake_curve_processing/NaCl-2')

FIXED_PARAMETERS_PATH = DATASET_PATH / 'drift_corrected.json'
FITTED_PARAMETERS_PATH = DATASET_PATH / 'fitted.json'

def model(t, y, p):
    X, X_gs, X_s, P = y
    V, R, F, X_in, k_diff, k_ads, k_des, k_rxn, S_tot, Y_tot = p
    return (
        F*(X_in - X) / V - k_diff * (X - X_gs),
        k_diff*(X - X_gs) - 2.0/R*k_ads*X_gs*(S_tot - X_s) + k_des*X_s,
        k_ads*X_gs*(S_tot - X_s) - R/2.0*k_des*X_s - k_rxn*X_s*(Y_tot - P),
        k_rxn*X_s*(Y_tot - P)
    )

def main():

    with open(FIXED_PARAMETERS_PATH, 'r') as f:
        fixed_parameters = json.load(f)

    with open(FITTED_PARAMETERS_PATH, 'r') as f:
        fitted_parameters = json.load(f)

    F = fixed_parameters['F'] * 760 / fixed_parameters['pressure']
    R = fixed_parameters['R']
    L = fixed_parameters['L']
    V = np.pi * R**2.0 * L
    X_in = fixed_parameters['X_in']
    Di = 40.0
    k_diff = 3.66 * Di / R**2.0
    k_ads = fitted_parameters['solution']['k_ads']
    k_des = fitted_parameters['solution']['k_des']
    k_rxn = fitted_parameters['solution']['k_rxn']
    S_tot = fitted_parameters['solution']['S_tot'] / 8  # TODO: remove the factor
    Y_tot = fitted_parameters['solution']['Y_tot']

    y0 = (X_in, X_in, 0.0, 0.0)
    p = (V, R, F, X_in, k_diff, k_ads, k_des, k_rxn, S_tot, Y_tot)

    ts = np.logspace(-2, 4, 50)

    sol = solve_ivp(lambda t, y: model(t,y,p),
                    t_span=(0.0, ts[-1]),
                    y0=y0,
                    method="BDF",   # or "Radau"
                    rtol=1e-8, atol=1e-12, t_eval=ts)

    uptake_rate = k_ads * sol.y[1] * (S_tot - sol.y[2]) - R/2.0 * k_des * sol.y[2]

    fig, ax = plt.subplots()

    l1, = ax.plot(sol.t, uptake_rate, label=R'$R_{\rm uptake}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Time, s')
    ax.set_ylabel(R'Rate of uptake, $1/(cm^2\cdot s)$')

    ax2 = ax.twinx()

    l2, = ax2.plot(sol.t, (S_tot - sol.y[2]), label='$S$', color='tab:orange')
    # ax2.set_yscale('log')
    ax2.set_ylabel('Site concentration, $1/cm^2$')

    ax.legend(handles=[l1, l2], loc="best")

    plt.show()

if __name__ == '__main__':
    main()
