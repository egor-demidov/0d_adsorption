import subprocess
from generate_data_with_noise import DURATIONS, PREFIX, SEEDS

EXECUTABLE_PATH = '../../build-release/0d_adsorption_fit_chained.exe'

for duration in DURATIONS:
    parent_dir = PREFIX / f'run_{duration}'
    for i in range(len(SEEDS)):
        run_dir = parent_dir / f'run_{duration}_noisy_{i+1}'
        input_file = run_dir / 'drift_corrected.json'

        print(f'Fitting file {input_file}')

        subprocess.run([EXECUTABLE_PATH, input_file])

