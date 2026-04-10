import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
import json

def model(t, y, p):
    # X1, X_gs1, X_s1, P1, X2, X_gs2, X_s2, P2 = y
    V, A, R, F, X_feed, k_diff, k_ads, k_des, k_rxn, S_tot, Y_tot, t_ads_start, t_ads_end, Y_tot_dist, k_rxn_dist = p
    k_ads_eff = k_ads * ((0.5 + 0.5 * np.tanh((t - t_ads_start) / 2.646)) - (0.5 + 0.5 * np.tanh((t - t_ads_end) / 2.646)))

    # assert np.isclose(np.sum(Y_tot_dist), Y_tot)
    # assert np.isclose(np.dot(k_rxn_dist,Y_tot_dist) / Y_tot, k_rxn)

    dy = np.zeros(y.shape)
    X, X_gs, X_s, P_dist = (y[0], y[1], y[2], y[3:])

    dy[0] = F/V*(X_feed - X) - k_diff * (X - X_gs)
    dy[1] = k_diff*(X - X_gs) - 2.0/R*k_ads_eff*X_gs*(S_tot - X_s) + 2.0/R*k_des*X_s
    # dy[2] = k_ads_eff*X_gs*(S_tot - X_s) - k_des*X_s - k_rxn*X_s*(P_tot - P)
    dy[2] = k_ads_eff*X_gs*(S_tot - X_s) - k_des*X_s - np.dot(k_rxn_dist, (Y_tot_dist - P_dist)) * X_s
    dy[3:] = k_rxn_dist*X_s*(Y_tot_dist - P_dist)

    return dy


def solve_ode(p0, t_span, NY_types):
    X0 = p0[4]

    y0 = np.zeros(3 + NY_types)
    y0[0] = X0    # IC for X
    y0[1] = X0      # IC for Xgs

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

    # Number of sorptive site types
    NY_types = 5
    Y_fracs = np.array([0.3, 0.3, 0.2, 0.15, 0.05])          # Length: NS_types; adds up to 1
    k_rxn_multipliers = np.array([1.0e-3, 1.0e-2, 0.5, 1.0])     # Length: NS_types - 1
    # k_rxn_multipliers = np.ones(NY_types - 1)
    k_rxn_last_multiplier = (1.0 - np.dot(Y_fracs[:-1], k_rxn_multipliers)) / Y_fracs[-1]
    Y_tot_dist = Y0 * Y_fracs
    k_rxn_dist = k_rxn * np.concatenate([k_rxn_multipliers, [k_rxn_last_multiplier]])

    print(k_rxn_dist)
    print(np.dot(k_rxn_dist, Y_tot_dist) / Y0, k_rxn)
    print(Y_tot_dist)
    print(np.sum(Y_tot_dist), Y0)

    print('k_ads multipliers:', np.concatenate([k_rxn_multipliers, [k_rxn_last_multiplier]]))

    # Perform checks
    assert len(Y_fracs) == NY_types
    assert len(k_rxn_multipliers) == NY_types - 1
    assert np.isclose(np.sum(Y_fracs), 1.0, atol=0)
    assert k_rxn_dist[-1] > 0.0
    assert np.isclose(np.sum(Y_tot_dist), Y0, atol=0)
    assert np.isclose(np.dot(k_rxn_dist, Y_tot_dist) / Y0, k_rxn, atol=0)

    plt.bar(k_rxn_dist, Y_tot_dist, align='center', width=1e-16)
    plt.axvline(k_rxn, ls='--', color='red')
    plt.text(k_rxn * 1.1, np.max(Y_tot_dist) * 0.9, 'mean', color='red')
    plt.xlabel(R'$k_{\rm rxn,i}$')
    plt.ylabel(R'$Y_{\rm tot,i}$')
    plt.show()

    p0 = (V, A, R, F, X0, k_diff, k_ads, k_des, k_rxn, S0, Y0, t_ads_start, t_ads_end, Y_tot_dist, k_rxn_dist)

    sol = solve_ode(p0, fixed_data['experimental_data']['t_exp'][-1], NY_types)

    plt.plot(sol.t, sol.y[0])
    # plt.plot(sol.t, sol.y[3] / Y_tot_dist[3])
    # plt.plot(sol.t, sol.y[-1] / Y_tot_dist[-1])
    plt.show()


if __name__ == '__main__':
    main()
