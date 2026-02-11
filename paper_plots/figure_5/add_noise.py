import json
import numpy as np
import matplotlib.pyplot as plt
import sys

def add_noise(noise_magnitude, seed, path_in, path_out):

    np.random.seed(seed=seed)

    with open(path_in, 'r') as in_f:
        in_data = json.load(in_f)

    t = np.array(in_data["solution"]["t"])
    X = np.array(in_data["solution"]["X"])

    noise = (np.random.rand(X.shape[0]) - 0.5) / 0.5 * noise_magnitude

    X += noise

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

    if len(sys.argv) < 5:
        print('Three arguments must be provided: noise magnitude, seed, input path, and output path', file=sys.stderr)
        exit(-1)

    NOISE_MAGNITUDE = float(sys.argv[1])
    SEED = int(sys.argv[2])
    PATH_IN = sys.argv[3]
    PATH_OUT = sys.argv[4]

    add_noise(NOISE_MAGNITUDE, SEED, PATH_IN, PATH_OUT)
