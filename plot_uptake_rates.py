import matplotlib.pyplot as plt
import json


DATASET = 'NaCl-2'

with open(f'uptake_curve_processing/{DATASET}/drift_corrected.json', 'r') as f:
    exp_data = json.load(f)

with open(f'uptake_curve_processing/{DATASET}/fitted.json', 'r') as f:
    fitted_data = json.load(f)

plt.plot(exp_data['experimental_data']['t_exp'], fitted_data['uptake_rates']['exact'], label='Exact')
plt.plot(exp_data['experimental_data']['t_exp'], fitted_data['uptake_rates']['initial_approx'], label='Initial')
plt.plot(exp_data['experimental_data']['t_exp'], fitted_data['uptake_rates']['ss_approx'], label='SS')
plt.legend()
plt.show()