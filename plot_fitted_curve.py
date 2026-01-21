import matplotlib.pyplot as plt
import json


with open('experimental_data.json', 'r') as f:
    exp_data = json.load(f)

with open('build-release/fitted_curve.json', 'r') as f:
    fitted_data = json.load(f)

plt.plot(exp_data['t_exp'], exp_data['X_exp'], label='Experimental')
plt.plot(exp_data['t_exp'], fitted_data['X0'], label='Initial guess')
plt.plot(exp_data['t_exp'], fitted_data['X'], label='Fitted curve')
plt.legend()
plt.show()