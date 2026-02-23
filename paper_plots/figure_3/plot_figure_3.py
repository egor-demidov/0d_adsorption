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



data = [
    {
        'name': 'nacl',
        'experimental': dict(),
        'fitted': dict()
    },
    {
        'name': 'levoglucosan',
        'experimental': dict(),
        'fitted': dict()
    }
]

# Load the datasets
for dataset in data:
    with open(f'{dataset['name']}_drift_corrected.json', 'r') as f:
        dataset['experimental'] = json.load(f)
    with open(f'{dataset['name']}_fitted.json', 'r') as f:
        dataset['fitted'] = json.load(f)


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.3))


for ax, dataset in zip((ax1, ax2), data):
    ax.plot(dataset['experimental']['experimental_data']['t_exp'], dataset['experimental']['experimental_data']['X_exp'], label='Experimental')
    ax.plot(dataset['experimental']['experimental_data']['t_exp'], dataset['fitted']['fitted_data']['X'], label='Fit')
    ax.plot(dataset['experimental']['experimental_data']['t_exp'], dataset['fitted']['fitted_data']['X0'], '--', color='gray', label='Initial guess')
    ax.set_xlabel(R'Time, $\rm s$')
    ax.set_ylabel(R'X concentration, $\rm cm^{-3}$')


ax1.legend()

# Add a,b,c,d labels
for ax, label in zip((ax1, ax2), string.ascii_lowercase):
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
plt.savefig('figure_3.pdf')
plt.show()
