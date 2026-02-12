import json
import numpy as np
import matplotlib.pyplot as plt
import sys

freq_lo = 0.3   # Number of noise points / experimental point (low frequency noise)
freq_hi = 1.0   # Number of noise points / experimental point (high frequency noise)

# Good values for nacl:
# cum. noise magnitude: 5.0e7
# rand. noise magnitude: 2.0e8
# freq_lo: 300
# high freq. rand. noise magnitude: 1.0e8
# freq_hi: 1000
# Seed used for testing: 486489
# Input used for testing: combo_1/run_600.json

def add_noise(cum_noise_magnitude, rand_noise_magnitude, hf_rand_noise_magnitude, seed, path_in, path_out):

    np.random.seed(seed=seed)

    with open(path_in, 'r') as in_f:
        in_data = json.load(in_f)

    t = np.array(in_data["solution"]["t"])
    X = np.array(in_data["solution"]["X"])

    M = int(freq_lo * X.shape[0])     # number of low frequency noise points
    L = int(freq_hi * X.shape[0])    # Number of high frequency noise points

    cum_noise = np.cumsum((np.random.rand(M) - 0.5) / 0.5 * cum_noise_magnitude)

    rand_noise = (np.random.normal(scale=1, size=M) - 0.5) / 0.5 * rand_noise_magnitude

    noise = cum_noise + rand_noise

    # Interpolate the noise filter
    x_n = np.linspace(0.0, 1.0, X.shape[0])
    x_m = np.linspace(0.0, 1.0, M)

    noise_interp = np.interp(x_n, x_m, noise)

    # Generate high frequency noise
    rand_noise_hf = (np.random.normal(scale=1, size=L) - 0.5) / 0.5 * hf_rand_noise_magnitude

    # Interpolate the high frequency noise filter
    x_l = np.linspace(0.0, 1.0, L)
    rand_noise_interp_hf = np.interp(x_n, x_l, rand_noise_hf)

    X += noise_interp + rand_noise_interp_hf

    with open(path_out, 'w') as out_f:
        out_data = {
            "solution": {
                "t": list(t),
                "X": list(X)
            }
        }
        json.dump(out_data, out_f, indent=4)

    # plt.plot(t, X)
    # plt.show()

if __name__ == '__main__':

    if len(sys.argv) < 7:
        print('Six arguments must be provided: cum. noise magnitude, rand. noise magnitude, HF rand. noise magnitude, seed, input path, and output path', file=sys.stderr)
        exit(-1)

    CUM_NOISE_MAGNITUDE = float(sys.argv[1])
    RAND_NOISE_MAGNITUDE = float(sys.argv[2])
    HF_RAND_NOISE_MAGNITUDE = float(sys.argv[3])
    SEED = int(sys.argv[4])
    PATH_IN = sys.argv[5]
    PATH_OUT = sys.argv[6]

    add_noise(CUM_NOISE_MAGNITUDE, RAND_NOISE_MAGNITUDE, HF_RAND_NOISE_MAGNITUDE, SEED, PATH_IN, PATH_OUT)
