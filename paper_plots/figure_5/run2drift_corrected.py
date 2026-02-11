import json
import sys


def run2drit_corrected(duration, path_in, path_out):

    PARAMETERS = {
        "F": 2.493333333333333,
        "R": 0.78,
        "L": 2,
        "X_feed": 29424160260.695,
        "N_reactors": 5,
        "pressure": 1.965,
        "t_ads_start": 100.0,
        "t_ads_end": 100.0 + duration,
        "k_ads_smooth": 2.0,
        "initial_guess": {
            "k_ads": 2.5e-12,
            "k_des": 0.02,
            "k_rxn": 1.0e-17,
            "S_tot": 6.0e13,
            "Y_tot": 1.0e14
        }
    }

    with open(path_in, 'r') as in_f:
        in_data = json.load(in_f)

    PARAMETERS["experimental_data"] = {}
    PARAMETERS["experimental_data"]["t_exp"] = in_data["solution"]["t"]
    PARAMETERS["experimental_data"]["X_exp"] = in_data["solution"]["X"]

    with open(path_out, 'w') as out_f:
        json.dump(PARAMETERS, out_f, indent=4)

if __name__ == '__main__':

    if len(sys.argv) < 4:
        print('Three arguments must be provided: duration, input path, and output path', file=sys.stderr)
        exit(-1)

    DURATION = float(sys.argv[1])
    PATH_IN = sys.argv[2]
    PATH_OUT = sys.argv[3]

    run2drit_corrected(DURATION, PATH_IN, PATH_OUT)
