import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
import json

def model(t, y, p):
    # X1, X_gs1, X_s1, P1, X2, X_gs2, X_s2, P2 = y
    V, A, R, F, X_feed, k_diff, k_ads, k_des, k_rxn, S_tot, P_tot, t_ads_start, t_ads_end, N_reactors = p
    k_ads_eff = k_ads * ((0.5 + 0.5 * np.tanh((t - t_ads_start) / 2.646)) - (0.5 + 0.5 * np.tanh((t - t_ads_end) / 2.646)))

    # Normalize volume and area by number of chained reactors
    V /= N_reactors
    A /= N_reactors

    dy = np.zeros(y.shape)
    for n in range(N_reactors):
        X_prev = X_feed if n == 0 else y[4*(n-1)]
        X, X_gs, X_s, P = y[4*n:4*n+4]
        dy[4*n] = (F*(X_prev - X) - V * k_diff * (X - X_gs)) / V
        dy[4*n+1] = (V * k_diff*(X - X_gs) - V * 2.0/R*k_ads_eff*X_gs*(S_tot - X_s) + V * k_des*X_s) / V
        dy[4*n+2] = (A * k_ads_eff*X_gs*(S_tot - X_s) - A * R/2.0*k_des*X_s - A * k_rxn*X_s*(P_tot - P)) / A
        dy[4*n+3] = (A * k_rxn*X_s*(P_tot - P)) / A

    return dy

    # return (
    #     (F*(X_feed - X1) - V * k_diff * (X1 - X_gs1)) / V,
    #     (V * k_diff*(X1 - X_gs1) - V * 2.0/R*k_ads_eff*X_gs1*(S_tot - X_s1) + V * k_des*X_s1) / V,
    #     (A * k_ads_eff*X_gs1*(S_tot - X_s1) - A * R/2.0*k_des*X_s1 - A * k_rxn*X_s1*(P_tot - P1)) / A,
    #     (A * k_rxn*X_s1*(P_tot - P1)) / A,
    #     (F*(X1 - X2) - V * k_diff * (X1 - X_gs2)) / V,
    #     (V * k_diff*(X2 - X_gs2) - V * 2.0/R*k_ads_eff*X_gs2*(S_tot - X_s2) + V * k_des*X_s2) / V,
    #     (A * k_ads_eff*X_gs2*(S_tot - X_s2) - A * R/2.0*k_des*X_s2 - A * k_rxn*X_s2*(P_tot - P2)) / A,
    #     (A * k_rxn*X_s2*(P_tot - P2)) / A
    # )

def solve_ode(p0, N_reactors, t_span):
    p = p0 + (N_reactors,)
    X0 = p[4]

    # Initialize the initial conditions
    y0 = np.zeros(N_reactors*4)
    for n in range(N_reactors):
        y0[n*4] = X0    # IC for X
        y0[n*4+1] = X0  # IC for Xgs

    sol = solve_ivp(lambda t, y: model(t,y,p),
                    t_span=(0.0, t_span),
                    y0=y0,
                    method="BDF",   # or "Radau"
                    rtol=1e-8, atol=1e-12,
                    # t_eval=data['solution_curves']['t_exp']
                    )

    return sol

def main():

    with open('uptake_curve_processing/levoglucosan/drift_corrected.json', 'r') as f:
        fixed_data = json.load(f)

    with open('uptake_curve_processing/levoglucosan/fitted.json', 'r') as f:
        fitted_data = json.load(f)

    F = fixed_data['F'] * 760 / fixed_data['pressure']
    R = 1.56 / 2.0
    L = fixed_data['L']
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

    p0 = (V, A, R, F, X0, k_diff, k_ads, k_des, k_rxn, S0, Y0, t_ads_start, t_ads_end)

    sols = (
               (solve_ode(p0, 1, fixed_data['experimental_data']['t_exp'][-1]), 1),
               (solve_ode(p0, 2, fixed_data['experimental_data']['t_exp'][-1]), 2),
               (solve_ode(p0, 3, fixed_data['experimental_data']['t_exp'][-1]), 3),
               (solve_ode(p0, 4, fixed_data['experimental_data']['t_exp'][-1]), 4),
               (solve_ode(p0, 5, fixed_data['experimental_data']['t_exp'][-1]), 5),
    )

    # plt.plot(data['solution_curves']['t_exp'], data['solution_curves']['X_exp'], 'o', ms=1)
    # plt.plot(data['solution_curves']['t_exp'], data['solution_curves']['X_interp'], label='Resolved model')
    for sol, N_reactors in sols:
        plt.plot(sol.t, sol.y[(N_reactors - 1)*4], label=f'{N_reactors} reactor(s)')
    plt.legend()
    # plt.plot(data.)
    plt.show()


if __name__ == '__main__':
    main()
