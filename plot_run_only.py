import matplotlib.pyplot as plt
import json


with open(f'run_only/NaCl/run.json', 'r') as f:
    run = json.load(f)

plt.plot(run['solution']['t'], run['solution']['X'], label='Experimental')
plt.xlabel('Time, s')
plt.ylabel(R'Concentration of X, $1/\rm cm^3$')
plt.legend()
# plt.savefig('levoglucosan-fit.png', dpi=300)
plt.show()