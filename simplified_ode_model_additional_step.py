import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
import json


k_ads_smooth = 2.5

def model(t, y, p):
    # X1, X_gs1, X_s1, P1, X2, X_gs2, X_s2, P2 = y
    V, A, R, F, X_feed, k_diff, k_ads, k_des, k_rxn, S_tot, Y_tot, t_ads_start, t_ads_end = p
    k_ads_eff = k_ads * ((0.5 + 0.5 * np.tanh((t - t_ads_start) / k_ads_smooth)) - (0.5 + 0.5 * np.tanh((t - t_ads_end) / k_ads_smooth)))

    # assert np.isclose(np.sum(Y_tot_dist), Y_tot)
    # assert np.isclose(np.dot(k_rxn_dist,Y_tot_dist) / Y_tot, k_rxn)

    dy = np.zeros(y.shape)
    I, X, X_gs, X_s, P = y

    dy[0] = (X - I) / 10.0
    dy[1] = F/V*(X_feed - X) - k_diff * (X - X_gs)
    dy[2] = k_diff * (X - X_gs) - 2.0 / R * k_ads_eff*X_gs*(S_tot - X_s) + 2.0 / R * k_des*X_s
    # dy[2] = k_ads_eff*X_gs*(S_tot - X_s) - k_des*X_s - k_rxn*X_s*(P_tot - P)
    dy[3] = k_ads_eff*X_gs*(S_tot - X_s) - k_des*X_s - k_rxn * (Y_tot - P) * X_s
    dy[4] = k_rxn * (Y_tot - P) * X_s

    return dy


def solve_ode(p0, t_span,):
    X0 = p0[4]

    y0 = np.zeros(5)
    y0[0] = X0    # IC for I
    y0[1] = X0    # IC for X
    y0[2] = X0      # IC for Xgs

    sol = solve_ivp(lambda t, y: model(t,y,p0),
                    t_span=(0.0, t_span),
                    y0=y0,
                    method="BDF",   # or "Radau"
                    rtol=1e-8, atol=1e-12,
                    t_eval=np.linspace(0.0, t_span, 300)
                    )

    return sol

def main():

    with open('paper_plots/figure_3/nacl_drift_corrected.json', 'r') as f:
        fixed_data = json.load(f)

    with open('paper_plots/figure_3/nacl_fitted.json', 'r') as f:
        fitted_data = json.load(f)

    F = fixed_data['F'] * 760 / fixed_data['pressure']
    R = 1.56 / 2.0
    L = fixed_data['L']
    V = np.pi * R**2.0 * L
    A = 2.0 * np.pi * R * L
    X0 = fixed_data['X_in']
    Di = fixed_data['Di'] * 760 / fixed_data['pressure']
    k_diff = 3.66 * Di / R**2.0
    k_ads = fitted_data['solution']['k_ads']
    k_des = fitted_data['solution']['k_des']
    k_rxn = fitted_data['solution']['k_rxn']
    S0 = fitted_data['solution']['S_tot']
    Y0 = fitted_data['solution']['Y_tot']
    t_ads_start = fixed_data['t_ads_start']
    t_ads_end = fixed_data['t_ads_end']

    p0 = (V, A, R, F, X0, k_diff, k_ads, k_des, k_rxn, S0, Y0, t_ads_start, t_ads_end)

    sol = solve_ode(p0, fixed_data['experimental_data']['t_exp'][-1])

    plt.plot(sol.t, sol.y[0])
    # plt.plot(sol.t, sol.y[3] / Y_tot_dist[3])
    # plt.plot(sol.t, sol.y[-1] / Y_tot_dist[-1])
    plt.show()


if __name__ == '__main__':
    main()
