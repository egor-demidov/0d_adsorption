import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from matplotlib.ticker import MaxNLocator


R = 0.78    # cm
L = 5.0     # cm

with open('levoglucosan_parameter_convergence.json', 'r') as f:
    data_plot_b = json.load(f)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.3))

# Prepare the data for plot (a)
uptake_curves = []

reactor_counts_plot_a = [1, 3, 10, 18]
for count in reactor_counts_plot_a:
    with open(f'run_{count}.json', 'r') as f:
        data = json.load(f)
        ts = data['solution']['t']
        uptake_curves.append(data['solution']['X'])


for curve, count in zip(uptake_curves, reactor_counts_plot_a):
    ax1.plot(ts, curve, label=f'{count} reactor{'s' if count > 1 else ''}')

# Load the 1D model solution
with open('run_1d_model.json', 'r') as f:
    data_1d_model = json.load(f)

ax1.plot(data_1d_model['solution_curves']['t_exp'], data_1d_model['solution_curves']['X_interp'], label='1D model')

# Create inset axes
axins = inset_axes(ax1, width="40%", height="30%", loc="lower right")

# Plot same data in inset
for curve, count in zip(uptake_curves, reactor_counts_plot_a):
    axins.plot(ts, curve)
axins.plot(data_1d_model['solution_curves']['t_exp'], data_1d_model['solution_curves']['X_interp'])

# Define zoom region
x1, x2 = 150, 230
y1, y2 = 1.5e9, 8.5e9
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

# Remove ticks (optional)
axins.set_xticks([])
axins.set_yticks([])

# Draw box and connectors
mark_inset(ax1, axins, loc1=2, loc2=4, fc="none", ec="0.5")

ax1.legend(loc='upper left')
ax1.set_xlabel(R'Time, $\rm s$')
ax1.set_ylabel(R'X concentration, $\rm cm^{-3}$')

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

# Convert lists into numpy arrays
for key in plot_b.keys():
    plot_b[key] = np.array(plot_b[key])

# Calculate relative error
for key in list(plot_b.keys())[1:]:
    plot_b[key] = (plot_b[key] - plot_b[key][-1]) / plot_b[key][-1] * 100

# Define transformation functions
def forward(N):
    return R * N / L

def inverse(x):
    return x * L / R

ax2.plot(plot_b["N_reactors"], plot_b["k_ads"], '-o', label=R'$k_{\rm ads}$')
ax2.plot(plot_b["N_reactors"], plot_b["k_des"], '-v', label=R'$k_{\rm des}$')
ax2.plot(plot_b["N_reactors"], plot_b["k_rxn"], '-^', label=R'$k_{\rm rxn}$')
ax2.plot(plot_b["N_reactors"], plot_b["S_tot"], '-s', label=R'$S_{\rm tot}$')
ax2.plot(plot_b["N_reactors"], plot_b["Y_tot"], '-*', label=R'$Y_{\rm tot}$')
ax2.legend()
ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

ax2.set_xlabel('Number of reactors')
ax2.set_ylabel('Relative parameter error, %')

secax = ax2.secondary_xaxis('top', functions=(forward, inverse))
secax.set_xlabel(R'$RN_{\rm reactors}/L$')

plt.tight_layout()
plt.savefig('figure_2.pdf')
plt.show()
