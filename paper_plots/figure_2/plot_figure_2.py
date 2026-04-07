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


R = 0.78    # cm
L_LEVO = 5.0     # cm
L_NACL = 2.0     # cm
F_LEVO = 2.367965 * 760 / 1.924             # cm^3/s
F_NACL = 2.493333333333333 * 760 / 1.965    # cm^3/s

with open('levoglucosan_parameter_convergence.json', 'r') as f:
    data_plot_b = json.load(f)

with open('nacl_parameter_convergence.json', 'r') as f:
    data_plot_d = json.load(f)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 4.3*1.8))

# Prepare the data for plot (a)
uptake_curves = []

reactor_counts_plot_a = [1, 3, 10]
for count in reactor_counts_plot_a:
    with open(f'levoglucosan_runs/run_{count}.json', 'r') as f:
        data = json.load(f)
        ts = data['solution']['t']
        uptake_curves.append(data['solution']['X'])


for curve, count in zip(uptake_curves, reactor_counts_plot_a):
    ax1.plot(ts, curve, label=f'{count} reactor{'s' if count > 1 else ''}')

# Load the 1D model solution
with open('levoglucosan_runs/run_1d_model.json', 'r') as f:
    data_1d_model = json.load(f)

ax1.plot(data_1d_model['solution_curves']['t_exp'], data_1d_model['solution_curves']['X_interp'], label='1D model')

ax1.set_xlim((100, 400))
# ax1.set_ylim((0.15e10, 1.25e10))

ax1.legend(loc='lower right')
ax1.set_xlabel(R'Time, $\rm s$')
ax1.set_ylabel(R'X concentration, $\rm cm^{-3}$')


# Prepare the data for plot (c)
uptake_curves = []

reactor_counts_plot_c = [1, 3, 10]
for count in reactor_counts_plot_c:
    with open(f'nacl_runs/run_{count}.json', 'r') as f:
        data = json.load(f)
        ts = data['solution']['t']
        uptake_curves.append(data['solution']['X'])


for curve, count in zip(uptake_curves, reactor_counts_plot_c):
    ax3.plot(ts, curve, label=f'{count} reactor{'s' if count > 1 else ''}')

# Load the 1D model solution
with open('nacl_runs/run_1d_model.json', 'r') as f:
    data_1d_model = json.load(f)

ax3.plot(data_1d_model['solution_curves']['t_exp'], data_1d_model['solution_curves']['X_interp'], label='1D model')

ax3.set_xlim((200, 700))
ax3.set_ylim((ax3.get_ylim()[0]*0.95, None))

# Create inset axes
axins = inset_axes(ax3, width="60%", height="45%", loc="lower right")

# Plot same data in inset
for curve, count in zip(uptake_curves, reactor_counts_plot_c):
    axins.plot(ts, curve)
axins.plot(data_1d_model['solution_curves']['t_exp'], data_1d_model['solution_curves']['X_interp'])

# Define zoom region
x1, x2 = 260, 300
y1, y2 = 1.7e10, 2.4e10
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

# Remove ticks (optional)
axins.set_xticks([])
axins.set_yticks([])

# Draw box and connectors
mark_inset(ax3, axins, loc1=2, loc2=4, fc="none", ec="0.5")

ax3.set_xlabel(R'Time, $\rm s$')
ax3.set_ylabel(R'X concentration, $\rm cm^{-3}$')

# Prepare the data for plot (b)
plot_b = {
    "N_reactors": [],
    "k_ads": [],
    "k_des": [],
    "k_rxn": [],
    "S_tot": [],
    "Y_tot": []
}

# Load the data into lists
for point in data_plot_b:
    for key in plot_b.keys():
        plot_b[key].append(point[key])

plot_b_true_k_ads = plot_b["k_ads"][-1]
plot_b_true_S_tot = plot_b["S_tot"][-1]

# Convert lists into numpy arrays
for key in plot_b.keys():
    plot_b[key] = np.array(plot_b[key])

# Calculate relative error
for key in list(plot_b.keys())[1:]:
    plot_b[key] = (plot_b[key] - plot_b[key][-1]) / plot_b[key][-1] * 100

# Define transformation functions
def invert(x):
    # 1/x with special treatment of x == 0
    x = np.array(x).astype(float)
    near_zero = np.isclose(x, 0)
    x[near_zero] = np.inf
    x[~near_zero] = 2.0 * np.pi * R * L_LEVO * plot_b_true_k_ads * plot_b_true_S_tot / (x[~near_zero] * F_LEVO)
    return x

# for i in range(len(plot_b["N_reactors"])):
#     print(plot_b["N_reactors"][i], invert(plot_b["N_reactors"][i]), plot_b["k_ads"][i])

ax2.plot(plot_b["N_reactors"], plot_b["k_ads"], '-o', color='tab:blue', label=R'$k_{\rm ads}$')
ax2.plot(plot_b["N_reactors"], plot_b["k_des"], '-v', color='tab:orange', label=R'$k_{\rm des}$')
# ax2.plot(plot_b["N_reactors"], plot_b["k_rxn"], '-^', color='tab:green', label=R'$k_{\rm rxn}$')
ax2.plot(plot_b["N_reactors"], plot_b["S_tot"], '-s', color='tab:red', label=R'$S_{\rm tot}$')
# ax2.plot(plot_b["N_reactors"], plot_b["Y_tot"], '-*', color='tab:purple', label=R'$Y_{\rm tot}$')
# ax2.legend()
ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

ax2.set_ylim(bottom=None, top=130)
# ax2.set_xlim(left=1, right=18)

ax2.set_xlabel('Number of reactors')
ax2.set_ylabel('Relative parameter error, %')

secax = ax2.secondary_xaxis('top', functions=(invert, invert))
secax.set_xlabel(R'${\rm Da}/N_{\rm reactors}$')
secax.set_xticks((1.0, 0.3, 0.15, 0.1))

# Prepare the data for plot (d)
plot_d = {
    "N_reactors": [],
    "k_ads": [],
    "k_des": [],
    "k_rxn": [],
    "S_tot": [],
    "Y_tot": []
}

# Load the data into lists
for point in data_plot_d:
    for key in plot_d.keys():
        plot_d[key].append(point[key])

plot_d_true_k_ads = plot_d["k_ads"][-1]
plot_d_true_S_tot = plot_d["S_tot"][-1]

# Convert lists into numpy arrays
for key in plot_d.keys():
    plot_d[key] = np.array(plot_d[key])

# Calculate relative error
for key in list(plot_d.keys())[1:]:
    plot_d[key] = (plot_d[key] - plot_d[key][-1]) / plot_d[key][-1] * 100

# Define transformation functions
def invert_nacl(x):
    # 1/x with special treatment of x == 0
    x = np.array(x).astype(float)
    near_zero = np.isclose(x, 0)
    x[near_zero] = np.inf
    x[~near_zero] = 2.0 * np.pi * R * L_NACL * plot_d_true_k_ads * plot_d_true_S_tot / (x[~near_zero] * F_NACL)
    return x

for i in range(len(plot_d["N_reactors"])):
    print(plot_d["N_reactors"][i], invert_nacl(plot_d["N_reactors"][i]), plot_d["k_ads"][i])

ax4.plot(plot_d["N_reactors"], plot_d["k_ads"], '-o', label=R'$k_{\rm ads}$')
ax4.plot(plot_d["N_reactors"], plot_d["k_des"], '-v', label=R'$k_{\rm des}$')
ax4.plot(plot_d["N_reactors"], plot_d["k_rxn"], '-^', label=R'$k_{\rm rxn}$')
ax4.plot(plot_d["N_reactors"], plot_d["S_tot"], '-s', label=R'$S_{\rm tot}$')
ax4.plot(plot_d["N_reactors"], plot_d["Y_tot"], '-*', label=R'$Y_{\rm tot}$')
ax4.legend()
ax4.xaxis.set_major_locator(MaxNLocator(integer=True))

ax4.set_ylim(bottom=None, top=100)

ax4.set_xlabel('Number of reactors')
ax4.set_ylabel('Relative parameter error, %')

secax_nacl = ax4.secondary_xaxis('top', functions=(invert_nacl, invert_nacl))
secax_nacl.set_xlabel(R'${\rm Da}/N_{\rm reactors}$')
secax_nacl.set_xticks((0.5, 0.15, 0.07, 0.05))

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


plt.tight_layout()
# plt.savefig('figure_2.pdf')
plt.show()
