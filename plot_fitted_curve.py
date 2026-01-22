import matplotlib.pyplot as plt
import json


DATASET = 'NaCl-5'

with open(f'uptake_curve_processing/{DATASET}/drift_corrected.json', 'r') as f:
    exp_data = json.load(f)

with open(f'uptake_curve_processing/{DATASET}/fitted.json', 'r') as f:
    fitted_data = json.load(f)

plt.plot(exp_data['experimental_data']['t_exp'], exp_data['experimental_data']['X_exp'], label='Experimental')
plt.plot(exp_data['experimental_data']['t_exp'], fitted_data['fitted_data']['X0'], label='Initial guess')
plt.plot(exp_data['experimental_data']['t_exp'], fitted_data['fitted_data']['X'], label='Fitted curve')
plt.legend()
plt.show()