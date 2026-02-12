import matplotlib.pyplot as plt
import json
import numpy as np


# with open(f'run_only/run.json', 'r') as f:
#     run = json.load(f)

with open(f'paper_plots/figure_5/combo_1/run_600_noise.json', 'r') as f:
    run = json.load(f)

fig, ax = plt.subplots()
ax.plot(run['solution']['t'], run['solution']['X'], label='Experimental')
ax.set_xlabel('Time, s')
ax.set_ylabel(R'Concentration of X, $1/\rm cm^3$')
ax2 = ax.twinx()
# ax2.plot(run['solution']['t'], run['sensitivities']['dXdk_ads'], '-k', label=R'$k_{\rm ads}$')
# ax2.plot(run['solution']['t'], run['sensitivities']['dXfwd'], '-k', label=R'$k_{\rm fwd}$')
plt.legend()
# plt.savefig('levoglucosan-fit.png', dpi=300)
plt.show()