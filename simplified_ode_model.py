import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
import json

def model(t, y, p):
    X, X_gs, X_s, P = y
    V, A, R, F, X0, k_diff, k_ads, k_des, k_rxn, S_tot, P_tot, t_ads_start, t_ads_end = p
    k_ads_eff = k_ads * ((0.5 + 0.5 * np.tanh((t - t_ads_start) / 2.646)) - (0.5 + 0.5 * np.tanh((t - t_ads_end) / 2.646)))
    return (
        (F*(X0 - X) - V * k_diff * (X - X_gs)) / V,
        (V * k_diff*(X - X_gs) - V * 2.0/R*k_ads_eff*X_gs*(S_tot - X_s) + V * k_des*X_s) / V,
        (A * k_ads_eff*X_gs*(S_tot - X_s) - A * R/2.0*k_des*X_s - A * k_rxn*X_s*(P_tot - P)) / A,
        (A * k_rxn*X_s*(P_tot - P)) / A
    )


def main():

    with open('uptake_curve_processing/NaCl-2/drift_corrected.json', 'r') as f:
        fixed_data = json.load(f)

    with open('uptake_curve_processing/NaCl-2/fitted.json', 'r') as f:
        fitted_data = json.load(f)

    F = fixed_data['F'] * 760 / fixed_data['pressure']
    R = 1.56 / 2.0
    L = 2.0
    V = np.pi * R**2.0 * L
    A = 2.0 * np.pi * R * L
    X0 = fixed_data['X_in']
    Di = 40.0
    k_diff = 3.66 * Di / R**2.0
    k_ads = fitted_data['solution']['k_ads']
    k_des = fitted_data['solution']['k_des']
    k_rxn = fitted_data['solution']['k_rxn']
    S0 = fitted_data['solution']['S_tot']
    Y0 = fitted_data['solution']['Y_tot']
    t_ads_start = fixed_data['t_ads_start']
    t_ads_end = fixed_data['t_ads_end']

    y0 = (X0, X0, 0.0, 0.0)
    p = (V, A, R, F, X0, k_diff, k_ads, k_des, k_rxn, S0, Y0, t_ads_start, t_ads_end)

    sol = solve_ivp(lambda t, y: model(t,y,p),
                    t_span=(0.0, fixed_data['experimental_data']['t_exp'][-1]),
                    y0=y0,
                    method="BDF",   # or "Radau"
                    rtol=1e-8, atol=1e-12,
                    # t_eval=data['solution_curves']['t_exp']
                    )

    # plt.plot(data['solution_curves']['t_exp'], data['solution_curves']['X_exp'], 'o', ms=1)
    # plt.plot(data['solution_curves']['t_exp'], data['solution_curves']['X_interp'], label='Resolved model')
    plt.plot(sol.t, sol.y[0], label='Bulk volume model')
    plt.legend()
    # plt.plot(data.)
    plt.show()


if __name__ == '__main__':
    main()
