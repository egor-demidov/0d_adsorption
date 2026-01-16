import numpy as np
import matplotlib.pyplot as plt
import json


data = np.loadtxt('cmake-build-debug/out_data.csv', delimiter=',', skiprows=1)

with open('/media/egor/Data/1d-adsorption/uptake_curve_processing/NaCl-2/fitted.json', 'r') as f:
    exp_data = json.load(f)

plt.plot(exp_data['solution_curves']['t_exp'], exp_data['solution_curves']['X_exp'])
plt.plot(data[:, 0], data[:, 1])
plt.show()
