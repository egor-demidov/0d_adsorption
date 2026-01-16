import numpy as np
import matplotlib.pyplot as plt
import json


data = np.loadtxt('build-release/out_data.csv', delimiter=',', skiprows=1)

with open('/media/egor/Data/1d-adsorption/uptake_curve_processing/NaCl-2/fitted.json', 'r') as f:
    exp_data = json.load(f)

fig, ax = plt.subplots()

# ax.plot(exp_data['solution_curves']['t_exp'], exp_data['solution_curves']['X_exp'])
ax.plot(data[:, 0], data[:, 1], 'k', label=R'$X$')
ax2 = ax.twinx()
ax2.plot(data[:, 0], data[:, 2], label=R'$k_{\rm ads}$')
ax2.plot(data[:, 0], data[:, 3], label=R'$k_{\rm des}$')
ax2.plot(data[:, 0], data[:, 4], label=R'$k_{\rm rxn}$')
ax2.plot(data[:, 0], data[:, 5], label=R'$S_0$')
ax2.plot(data[:, 0], data[:, 6], label=R'$Y_0$')

ax.set_xlabel(R'$t$, $\rm s$')
ax.set_ylabel(R'$X$, $\rm 1/cm^3$')
ax2.set_ylabel(R'$\theta\cdot(\partial X/\partial \theta$), $\rm 1/cm^3$')
plt.legend()
plt.tight_layout()
plt.savefig('sensitivities.png', dpi=300)
plt.show()
