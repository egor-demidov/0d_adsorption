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


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 4.3*1.8))

# prepare data for plot (a)
uptake_curves = []
run_durations = [100, 200, 300, 400, 500, 600]
for duration in run_durations:
    with open(f'combo_1/run_{duration}_noise.json', 'r') as f:
        data = json.load(f)
        ts = data['solution']['t']
        uptake_curves.append(data['solution']['X'])

for curve in uptake_curves:
    ax1.plot(ts, curve)

# prepare data for plot (b)
plot_b = {
    "duration": [],
    "k_ads": [],
    "k_des": [],
    "k_rxn": [],
    "S_tot": [],
    "Y_tot": []
}

with open('combo_1/parameter_convergence.json', 'r') as f:
    data_plot_b = json.load(f)

# Load the data into lists
for point in data_plot_b:
    for key in plot_b.keys():
        plot_b[key].append(point[key])

ax2.plot(plot_b["duration"], plot_b["k_rxn"])

plt.tight_layout()
plt.show()
