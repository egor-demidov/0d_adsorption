import numpy as np
import matplotlib.pyplot as plt
import json
from pathlib import Path


SOL_PATH_1D = Path(R'C:\Users\mail\CLionProjects\1d-adsorption\uptake_curve_processing\levoglucosan\fitted.json')
SOL_PATH_0D = Path('uptake_curve_processing/levoglucosan/fitted.json')
INP_PATH_0D = Path('uptake_curve_processing/levoglucosan/drift_corrected.json')

with open(SOL_PATH_1D, 'r') as f:
    sol_1d = json.load(f)

with open(SOL_PATH_0D, 'r') as f:
    sol_0d = json.load(f)

with open(INP_PATH_0D, 'r') as f:
    inp_0d = json.load(f)

plt.plot(sol_1d['solution_curves']['t_exp'], sol_1d['solution_curves']['X_interp_0'], label='1D')
plt.plot(inp_0d['experimental_data']['t_exp'], sol_0d['fitted_data']['X'], label='0D')
plt.xlabel('Time, s')
plt.ylabel(R'Concentration of X, $\rm 1/cm^3$')
plt.legend()
# plt.savefig('levoglucosan-performance.png', dpi=300)
plt.show()
