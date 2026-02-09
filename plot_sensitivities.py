import matplotlib.pyplot as plt
import json
import numpy as np


DATASET = 'NaCl-2'

with open(f'uptake_curve_processing/{DATASET}/drift_corrected.json', 'r') as f:
    exp_data = json.load(f)

with open(f'uptake_curve_processing/{DATASET}/fitted.json', 'r') as f:
    fitted_data = json.load(f)


fig, ax = plt.subplots()

ax.plot(exp_data['experimental_data']['t_exp'], fitted_data['fitted_data']['X'], 'k', label=R'$X$')
ax2 = ax.twinx()
# ax2.plot(exp_data['experimental_data']['t_exp'], np.array(fitted_data['sensitivities']['dXdk_ads']) * fitted_data['solution']['k_ads'], label=R'$k_{\rm ads}$')
# ax2.plot(exp_data['experimental_data']['t_exp'], np.array(fitted_data['sensitivities']['dXdk_des']) * fitted_data['solution']['k_des'], label=R'$k_{\rm des}$')
# ax2.plot(exp_data['experimental_data']['t_exp'], np.array(fitted_data['sensitivities']['dXdk_rxn']) * fitted_data['solution']['k_rxn'], label=R'$k_{\rm rxn}$')
# ax2.plot(exp_data['experimental_data']['t_exp'], np.array(fitted_data['sensitivities']['dXdS_tot']) * fitted_data['solution']['S_tot'], label=R'$S_{\rm tot}$')
ax2.plot(exp_data['experimental_data']['t_exp'], np.array(fitted_data['sensitivities']['dXdY_tot']) * fitted_data['solution']['Y_tot'], label=R'$Y_{\rm tot}$')
ax2.set_ylim((-7.5e9, 5.5e9))

ax.set_xlabel(R'$t$, $\rm s$')
ax.set_ylabel(R'$X$, $\rm 1/cm^3$')
ax2.set_ylabel(R'$\theta\cdot(\partial X/\partial \theta$), $\rm 1/cm^3$')
plt.legend()
plt.tight_layout()
# plt.savefig('sensitivities_frame_6.png', dpi=300)
plt.show()
