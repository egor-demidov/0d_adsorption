import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from matplotlib.ticker import MaxNLocator
import string
import os

from generate_data_with_noise import DURATIONS, PREFIX, SEEDS


k_rxn_true = 2.5e-16
Y_tot_true = 9.0e13

plt.rcParams.update({
    "font.size": 12,        # default text size
    # "axes.titlesize": 14,
    # "axes.labelsize": 13,
    # "xtick.labelsize": 11,
    # "ytick.labelsize": 11,
    # "legend.fontsize": 11
})


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 4.3*1.8))

# prepare data for plot (a)
uptake_curves = []
run_durations = [100, 200, 300, 400, 500, 600]
for duration in run_durations:
    with open(f'combo_1/run_{duration}.json', 'r') as f:
        data = json.load(f)
        ts = data['solution']['t']
        uptake_curves.append(data['solution']['X'])

for curve in uptake_curves:
    ax1.plot(ts, curve)

ax1.plot([100.0, 100.0], [2.8e10, 2.94e10], '-k')
ax1.plot([700.0, 700.0], [2.8e10, 2.94e10], '-k')
ax1.plot([100.0, 700.0], [2.8e10, 2.8e10], '--k')
ax1.text(400.0, 2.73e10, 'Exposure', ha='center', va='top')

# Create inset axes
axins = inset_axes(ax1, width="60%", height="45%", loc="lower right")

# Read data file with noise
with open('combo_1/run_600_noise.json', 'r') as f:
    data_noise = json.load(f)

axins.plot(data_noise['solution']['t'], data_noise['solution']['X'], color='tab:brown')

# Remove ticks (optional)
axins.set_xticks([])
axins.set_yticks([])

ax1.set_xlabel(R'Time, $\rm s$')
ax1.set_ylabel(R'X concentration, $\rm cm^{-3}$')

# prepare data for plot (b)
# Load fitted parameters

k_rxn = [[] for _ in DURATIONS]
Y_tot = [[] for _ in DURATIONS]
k_rxn_comp_err = [[] for _ in DURATIONS]
Y_tot_comp_err = [[] for _ in DURATIONS]

for duration, j in zip(DURATIONS, range(len(DURATIONS))):
    parent_dir = PREFIX / f'run_{duration}'
    for i in range(len(SEEDS)):
        run_dir = parent_dir / f'run_{duration}_noisy_{i+1}'
        input_file = run_dir / 'fitted.json'

        if os.path.exists(input_file):
            with open(input_file, 'r') as f:
                data = json.load(f)

                k_rxn[j].append(data['solution']['k_rxn'])
                Y_tot[j].append(data['solution']['Y_tot'])

                k_rxn_comp_err[j].append(data['standard_error']['k_rxn'])
                Y_tot_comp_err[j].append(data['standard_error']['Y_tot'])

k_rxn_mean = np.array([np.mean(l) for l in k_rxn])
Y_tot_mean = np.array([np.mean(l) for l in Y_tot])
k_rxn_stdev = np.array([np.std(l) for l in k_rxn])
Y_tot_stdev = np.array([np.std(l) for l in Y_tot])
k_rxn_comp_err_mean = np.array([np.mean(l) for l in k_rxn_comp_err])
Y_tot_comp_err_mean = np.array([np.mean(l) for l in Y_tot_comp_err])
k_rxn_bias = k_rxn_mean - k_rxn_true
Y_tot_bias = Y_tot_mean - Y_tot_true
counts = np.array(len(l) for l in k_rxn)

k_rxn_rel_stdev = k_rxn_stdev / k_rxn_true * 100.0
Y_tot_rel_stdev = Y_tot_stdev / Y_tot_true * 100.0

k_rxn_rel_bias = k_rxn_bias / k_rxn_true * 100.0
Y_tot_rel_bias = Y_tot_bias / Y_tot_true * 100.0

k_rxn_rmse = np.sqrt(k_rxn_rel_stdev ** 2.0 + k_rxn_rel_bias ** 2.0)
Y_tot_rmse = np.sqrt(Y_tot_rel_stdev ** 2.0 + Y_tot_rel_bias ** 2.0)

ax2.plot(DURATIONS, k_rxn_rel_stdev, '-^', color='tab:green', label=R'$k_{\rm rxn}$')

# ax2_twin = ax2.twinx()
ax2.plot(DURATIONS, Y_tot_rel_stdev, '-*', color='tab:purple', label=R'$Y_{\rm tot}$')

ax2.set_ylabel('Relative standard deviation, %')
ax2.set_xlabel(R'Exposure time, $\rm s$')
ax2.legend()

ax3.plot(DURATIONS, k_rxn_rel_bias, '-^', color='tab:green', label=R'$k_{\rm rxn}$')
ax3.plot(DURATIONS, Y_tot_rel_bias, '-*', color='tab:purple', label=R'$Y_{\rm tot}$')

ax3.set_ylabel('Relative bias, %')
ax3.set_xlabel(R'Exposure time, $\rm s$')
ax3.legend()

ax4.plot(DURATIONS, k_rxn_rmse, '-^', color='tab:green', label=R'$k_{\rm rxn}$')
ax4.plot(DURATIONS, Y_tot_rmse, '-*', color='tab:purple', label=R'$Y_{\rm tot}$')

ax4.set_ylabel('R.M.S. error, %')
ax4.set_xlabel(R'Exposure time, $\rm s$')
ax4.legend()

# Add a,b,c,d labels
for ax, label in zip((ax1, ax2, ax3, ax4), string.ascii_lowercase):
    ax.text(
        0.04, 0.96,                # position (x,y) in axes coords
        f'({label})',
        transform=ax.transAxes,    # use axes coordinates (0–1)
        fontsize=13,
        fontweight='bold',
        va='top',
        ha='left'
    )

ax2.set_ylim(bottom=None, top=57.0)
ax3.set_ylim(bottom=None, top=43.0)
ax4.set_ylim(bottom=None, top=69.0)

plt.tight_layout()
plt.savefig('figure_5.pdf')
plt.show()
