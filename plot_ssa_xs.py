import numpy as np
import matplotlib.pyplot as plt
import json


with open('uptake_curve_processing/NaCl-2/drift_corrected.json') as f:
    exp_data = json.load(f)

with open('uptake_curve_processing/NaCl-2/fitted.json') as f:
    fit_data = json.load(f)

ts = np.array(exp_data['experimental_data']['t_exp'])
X = np.array(fit_data['fitted_data']['X'])
Xs = np.array(fit_data['fitted_data']['Xs'])
P = np.array(fit_data['fitted_data']['P'])

k_ads = fit_data['solution']['k_ads']
k_des = fit_data['solution']['k_des']
k_rxn = fit_data['solution']['k_rxn']
S_tot = fit_data['solution']['S_tot']
Y_tot = fit_data['solution']['Y_tot']
R = exp_data['R']
Di = 40.0
k_diff = 3.66 * Di / R**2.0

K = k_rxn * (Y_tot - P)
E = R / 2.0 * k_des + K
T = 4.0 * K * S_tot * k_ads + R * k_diff * (2.0 * K + R * k_des + 2.0 * X * k_ads)

Xs_SSA = (T - np.sqrt(T**2.0 - 32.0*K*R*S_tot*X*k_ads**2.0*k_diff)) / (8.0 * k_ads * K)

K_Y0 = k_rxn * Y_tot
E_Y0 = R / 2.0 * k_des + K_Y0
T_Y0 = 4.0 * K_Y0 * S_tot * k_ads + R * k_diff * (2.0 * K_Y0 + R * k_des + 2.0 * X * k_ads)
Xs_SSA_Y0 = (T_Y0 - np.sqrt(T_Y0**2.0 - 32.0*K_Y0*R*S_tot*X*k_ads**2.0*k_diff)) / (8.0 * k_ads * K_Y0)

plt.plot(ts, Xs)
plt.plot(ts, Xs_SSA)
plt.plot(ts, Xs_SSA_Y0)
plt.show()
