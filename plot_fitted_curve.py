import matplotlib.pyplot as plt
import json
import numpy as np


DATASET = 'NaCl-2'

with open(f'uptake_curve_processing/{DATASET}/drift_corrected.json', 'r') as f:
    exp_data = json.load(f)

with open(f'uptake_curve_processing/{DATASET}/fitted.json', 'r') as f:
    fitted_data = json.load(f)

plt.plot(exp_data['experimental_data']['t_exp'], exp_data['experimental_data']['X_exp'], label='Experimental')
plt.plot(exp_data['experimental_data']['t_exp'], fitted_data['fitted_data']['X0'], label='Initial guess')
plt.plot(exp_data['experimental_data']['t_exp'], fitted_data['fitted_data']['X'], label='Fitted curve')
plt.xlabel('Time, s')
plt.ylabel(R'Concentration of X, $1/\rm cm^3$')
plt.legend()
# plt.savefig('levoglucosan-fit.png', dpi=300)
# plt.savefig('uptake-coeff-fig-1.png', dpi=300)


R_G = 8.314     # J/mol/K
T = 300.0       # K
M_W = 271.52e-3 # kg/mol

# Collision velocity
OMEGA = np.sqrt(R_G * T / 2.0 / np.pi / M_W) # m/s
OMEGA *= 1.0e2  # cm/s

# plt.figure()
# plt.plot(exp_data['experimental_data']['t_exp'], np.array(fitted_data['fitted_data']['uptake_rate']) / np.array(fitted_data['fitted_data']['Xgs']) / OMEGA)
# plt.xlabel('Time, s')
# plt.ylabel('Uptake coefficient')
# plt.savefig('uptake-coeff-fig-2.png', dpi=300)

plt.show()