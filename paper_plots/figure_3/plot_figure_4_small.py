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
    ],
    figsize=(10.0, 4.3),
)

twin_axes = [ax.twinx() for ax in axes.values()]

data_index = [
    (axes['A'], twin_axes[0], 'k_ads', R'$k_{\rm ads}$', 'tab:blue'),
    (axes['A'], twin_axes[0], 'k_des', R'$k_{\rm des}$', 'tab:orange'),
    (axes['B'], twin_axes[1], 'k_rxn', R'$k_{\rm rxn}$', 'tab:green'),
    (axes['A'], twin_axes[0], 'S_tot', R'$S_{\rm tot}$', 'tab:red'),
    (axes['B'], twin_axes[1], 'Y_tot', R'$Y_{\rm tot}$', 'tab:purple')
]

for ax in axes.values():

    ax.plot(exp_data['experimental_data']['t_exp'], fitted_data['fitted_data']['X'], '-k')
    ax.set_ylabel(R'X concentration, $\rm cm^{-3}$')
    ax.set_xlabel(R'Time, $\rm s$')
    ax.set_xlim((200.0, 800.0))

for ax, twinax, parameter, label, color in data_index:

    curve = f'dXd{parameter}'

    parameter_value = fitted_data['solution'][parameter]

    twinax.plot(exp_data['experimental_data']['t_exp'], np.array(fitted_data['sensitivities'][curve]) * parameter_value, color=color, label=label, lw='3', alpha=0.7)
    twinax.set_ylim(np.array([-75.0, 60.0]) * 1.0e8)


def align_y_values(ax1, y1_value, ax2, y2_value=0):
    """
    Align y1_value on ax1 with y2_value on ax2.
    Common use: align a custom baseline on the left axis with zero on the right axis.
    """

    y1min, y1max = ax1.get_ylim()
    y2min, y2max = ax2.get_ylim()

    # Position of y1_value within ax1, from 0 bottom to 1 top
    p = (y1_value - y1min) / (y1max - y1min)

    if not 0 < p < 1:
        raise ValueError("y1_value must be inside ax1 y-limits to align visibly.")

    # Choose new ax2 limits so y2_value appears at the same relative height p
    scale = max(
        (y2_value - y2min) / p,
        (y2max - y2_value) / (1 - p)
    )

    ax2.set_ylim(
        y2_value - p * scale,
        y2_value + (1 - p) * scale
    )

twin_axes[0].set_yticks([])
twin_axes[0].set_ylabel('')

axes['B'].set_yticks([])
axes['B'].set_ylabel('')

align_y_values(axes['A'], fitted_data['fitted_data']['X'][0], twin_axes[0])
align_y_values(axes['B'], fitted_data['fitted_data']['X'][0], twin_axes[1])

twin_axes[1].set_ylabel(f'{R'$\theta \cdot \partial X/\partial \theta$'}, {R'$\rm cm^{-3}$'}')

for ax in twin_axes:
    ax.legend(loc='lower right')

for ax in axes.values():
    ax.axhline(fitted_data['fitted_data']['X'][0], ls='--', color='k')

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
plt.savefig('figure_4_small.pdf')
plt.show()
