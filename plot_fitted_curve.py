import matplotlib.pyplot as plt
import json


DATASET = 'levoglucosan'

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
plt.show()