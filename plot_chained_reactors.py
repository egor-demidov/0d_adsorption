import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


RUN_DIR = Path('run_only/NaCl')

for i in range(1, 6):
    with open(RUN_DIR / f'run-{i}.json', 'r') as f:
        data = json.load(f)
        plt.plot(data['solution']['t'], data['solution']['X'])

plt.show()
