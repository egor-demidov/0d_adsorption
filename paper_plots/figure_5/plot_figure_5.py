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
    with open(f'combo_1/run_{duration}.json', 'r') as f:
        data = json.load(f)
        ts = data['solution']['t']
        uptake_curves.append(data['solution']['X'])

for curve in uptake_curves:
    ax1.plot(ts, curve)

plt.tight_layout()
plt.show()
