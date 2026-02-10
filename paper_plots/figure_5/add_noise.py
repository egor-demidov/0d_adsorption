import json
import numpy as np
import matplotlib.pyplot as plt


DURATION = 100
NOISE_MAGNITUDE = 5.0e8
SUFFIX = 4

np.random.seed(seed=7899515)

PATH_IN = f'combo_1/run_{DURATION}.json'
PATH_OUT = f'combo_1/run_{DURATION}_noise_{SUFFIX}.json'

def main():

    with open(PATH_IN, 'r') as in_f:
        in_data = json.load(in_f)

    t = np.array(in_data["solution"]["t"])
    X = np.array(in_data["solution"]["X"])

    noise = (np.random.rand(X.shape[0]) - 0.5) / 0.5 * NOISE_MAGNITUDE

    X += noise

    with open(PATH_OUT, 'w') as out_f:
        out_data = {
            "solution": {
                "t": list(t),
                "X": list(X)
            }
        }
        json.dump(out_data, out_f, indent=4)

    plt.plot(t, X)
    plt.show()

if __name__ == '__main__':
    main()
