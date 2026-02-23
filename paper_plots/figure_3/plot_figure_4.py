import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from matplotlib.ticker import MaxNLocator
import string


plt.rcParams.update({
    "font.size": 12,        # default text size
    # "axes.titlesize": 14,
    # "axes.labelsize": 13,
    # "xtick.labelsize": 11,
    # "ytick.labelsize": 11,
    # "legend.fontsize": 11
})



with open(f'nacl_drift_corrected.json', 'r') as f:
    exp_data = json.load(f)

with open(f'nacl_fitted.json', 'r') as f:
    fitted_data = json.load(f)

fig, axes = plt.subplot_mosaic(
    [
        ["A", "A" ,"B", "B"],
        ["C", "C", "D", "D"],
        [".", "E", "E", "."]
    ],
    figsize=(10.0, 4.3*2.4),
)

data_index = [
    (axes['A'], 'k_ads', R'$k_{\rm ads}\cdot \partial X/\partial k_{\rm ads}$', 'tab:blue'),
    (axes['B'], 'k_des', R'$k_{\rm des}\cdot \partial X/\partial k_{\rm des}$', 'tab:orange'),
    (axes['C'], 'k_rxn', R'$k_{\rm rxn}\cdot \partial X/\partial k_{\rm rxn}$', 'tab:green'),
    (axes['D'], 'S_tot', R'$S_{\rm tot}\cdot \partial X/\partial S_{\rm tot}$', 'tab:red'),
    (axes['E'], 'Y_tot', R'$Y_{\rm tot}\cdot \partial X/\partial Y_{\rm tot}$', 'tab:purple')
]

for ax, parameter, label, color in data_index:

    curve = f'dXd{parameter}'

    parameter_value = fitted_data['solution'][parameter]

    ax.plot(exp_data['experimental_data']['t_exp'], fitted_data['fitted_data']['X'], '-k')
    ax2 = ax.twinx()
    ax2.plot(exp_data['experimental_data']['t_exp'], np.array(fitted_data['sensitivities'][curve]) * parameter_value, color=color)
    ax2.set_ylim(np.array([-15.0, 12.0]) * 1.0e8)

    ax.set_ylabel(R'X concentration, $\rm cm^{-3}$')
    ax.set_xlabel(R'Time, $\rm s$')
    ax2.set_ylabel(f'{label}, {R'$\rm cm^{-3}$'}')

# for ax, dataset in zip((ax1, ax2), data):
#     ax.plot(dataset['experimental']['experimental_data']['t_exp'], dataset['experimental']['experimental_data']['X_exp'], label='Experimental')
#     ax.plot(dataset['experimental']['experimental_data']['t_exp'], dataset['fitted']['fitted_data']['X'], label='Fit')
#     ax.plot(dataset['experimental']['experimental_data']['t_exp'], dataset['fitted']['fitted_data']['X0'], '--', color='gray', label='Initial guess')
#     ax.set_xlabel(R'Time, $\rm s$')
#     ax.set_ylabel(R'X concentration, $\rm cm^{-3}$')



# Add a,b,c,d labels
for ax, label in zip(axes.values(), string.ascii_lowercase):
    ax.text(
        0.04, 0.96,                # position (x,y) in axes coords
        f'({label})',
        transform=ax.transAxes,    # use axes coordinates (0–1)
        fontsize=13,
        fontweight='bold',
        va='top',
        ha='left'
    )


plt.tight_layout()
plt.savefig('figure_4.pdf')
plt.show()
