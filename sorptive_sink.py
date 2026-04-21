import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
import json

def model(t, y, p):
    # X1, X_gs1, X_s1, P1, X2, X_gs2, X_s2, P2 = y
    V, A, R, F, X_feed, k_diff, k_ads, k_des, k_rxn, S_tot, P_tot, t_ads_start, t_ads_end, k_ads_sink, k_des_sink, S_tot_sink, V_sink, A_sink, N_reactors = p
    k_ads_eff = k_ads * ((0.5 + 0.5 * np.tanh((t - t_ads_start) / 2.646)) - (0.5 + 0.5 * np.tanh((t - t_ads_end) / 2.646)))

    # Normalize volume and area by number of chained reactors
    V /= N_reactors
    A /= N_reactors

    dy = np.zeros(y.shape)
    for n in range(N_reactors):
        X_prev = X_feed if n == 0 else y[4*(n-1)]
        X, X_gs, X_s, P = y[4*n:4*n+4]
        dy[4*n] = (F*(X_prev - X) - V * k_diff * (X - X_gs)) / V
        dy[4*n+1] = (V * k_diff*(X - X_gs) - V * 2.0/R*k_ads_eff*X_gs*(S_tot - X_s) + V * 2.0/R * k_des*X_s) / V
        dy[4*n+2] = (A * k_ads_eff*X_gs*(S_tot - X_s) - A*k_des*X_s - A * k_rxn*X_s*(P_tot - P)) / A
        dy[4*n+3] = (A * k_rxn*X_s*(P_tot - P)) / A

    X_prev = y[4*(N_reactors-1)]
    X, X_gs, X_s = y[4*N_reactors:]
    SV_ratio_sink = A_sink / V_sink


    dy[4*N_reactors] = F/V_sink*(X_prev - X) - k_diff * (X - X_gs)
    dy[4*N_reactors+1] = k_diff*(X - X_gs) - SV_ratio_sink*k_ads_sink*X_gs*(S_tot_sink - X_s) + SV_ratio_sink * k_des_sink*X_s
    dy[4*N_reactors+2] = k_ads_sink*X_gs*(S_tot_sink - X_s) - k_des_sink*X_s

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
    k_ads_sink, k_des_sink, S_tot_sink = (p[-6], p[-5], p[-4])
    K_ads_sink = k_ads_sink / k_des_sink

    print(S_tot_sink, k_ads_sink, k_des_sink)

    # Initialize the initial conditions
    y0 = np.zeros(N_reactors*4+3)
    for n in range(N_reactors):
        y0[n*4] = X0    # IC for X
        y0[n*4+1] = X0  # IC for Xgs

    # Initial condition for the sink and sink wall
    # TODO: uptade to langmuir later
    y0[N_reactors*4] = X0   # X
    y0[N_reactors*4+1] = X0    # Xgs
    y0[N_reactors*4+2] = S_tot_sink * K_ads_sink * X0 / (1.0 + K_ads_sink * X0)    # Xs

    print(y0)

    sol = solve_ivp(lambda t, y: model(t,y,p),
                    t_span=(0.0, t_span),
                    y0=y0,
                    method="BDF",   # or "Radau"
                    rtol=1e-8, atol=1e-12,
                    t_eval=np.linspace(0.0, t_span, 300)
                    )

    return sol

def main():

    with open('uptake_curve_processing/NaCl/drift_corrected.json', 'r') as f:
        fixed_data = json.load(f)

    with open('uptake_curve_processing/NaCl/fitted.json', 'r') as f:
        fitted_data = json.load(f)

    with open('uptake_curve_processing/NaCl/drift_corrected.json', 'r') as f:
        fixed_data_lo_conc = json.load(f)

    F = fixed_data_lo_conc['F'] * 760 / fixed_data_lo_conc['pressure']
    R = fixed_data_lo_conc['R'] / 2.0
    L = fixed_data_lo_conc['L']
    V = np.pi * R**2.0 * L
    A = 2.0 * np.pi * R * L
    X0 = fixed_data_lo_conc['X_feed']
    Di = fixed_data_lo_conc['Di'] * 760 / fixed_data_lo_conc['pressure']
    k_diff = 3.66 * Di / R**2.0
    k_ads = fitted_data['solution']['k_ads']
    k_des = fitted_data['solution']['k_des']
    k_rxn = fitted_data['solution']['k_rxn']
    S0 = fitted_data['solution']['S_tot']
    Y0 = fitted_data['solution']['Y_tot']
    t_ads_start = fixed_data_lo_conc['t_ads_start']
    t_ads_end = fixed_data_lo_conc['t_ads_end']

    # Sink parameters
    k_ads_sink = k_ads * 4.0
    k_des_sink = k_des
    S_tot_sink = S0 * 10.0
    V_sink = V
    A_sink = A

    p0 = (V, A, R, F, X0, k_diff, k_ads, k_des, k_rxn, S0, Y0, t_ads_start, t_ads_end, k_ads_sink, k_des_sink, S_tot_sink, V_sink, A_sink)

    N_reactors = 10

    sol = solve_ode(p0, N_reactors, fixed_data_lo_conc['experimental_data']['t_exp'][-1])

    # plt.plot(data['solution_curves']['t_exp'], data['solution_curves']['X_exp'], 'o', ms=1)
    # plt.plot(data['solution_curves']['t_exp'], data['solution_curves']['X_interp'], label='Resolved model')
    plt.plot(sol.t, sol.y[(N_reactors - 1)*4], label=f'No sink')
    plt.plot(sol.t, sol.y[N_reactors*4], label=f'Sink')
    plt.legend()
    # plt.plot(data.)
    plt.show()



if __name__ == '__main__':
    main()
