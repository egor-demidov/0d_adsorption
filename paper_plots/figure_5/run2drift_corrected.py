import json


DURATION = 100
SUFFIX = 4

PARAMETERS = {
    "F": 2.493333333333333,
    "R": 0.78,
    "L": 2,
    "X_feed": 29424160260.695,
    "N_reactors": 5,
    "pressure": 1.965,
    "t_ads_start": 100.0,
    "t_ads_end": 100.0 + DURATION,
    "k_ads_smooth": 2.0,
    "initial_guess": {
        "k_ads": 2.5e-12,
        "k_des": 0.02,
        "k_rxn": 1.0e-17,
        "S_tot": 6.0e13,
        "Y_tot": 1.0e14
    }
}

PATH_IN = f'combo_1/run_{DURATION}_noise_{SUFFIX}.json'
PATH_OUT = f'combo_1/run_{DURATION}_fit/drift_corrected_{SUFFIX}.json'

def main():

    with open(PATH_IN, 'r') as in_f:
        in_data = json.load(in_f)

    PARAMETERS["experimental_data"] = {}
    PARAMETERS["experimental_data"]["t_exp"] = in_data["solution"]["t"]
    PARAMETERS["experimental_data"]["X_exp"] = in_data["solution"]["X"]

    with open(PATH_OUT, 'w') as out_f:
        json.dump(PARAMETERS, out_f, indent=4)

if __name__ == '__main__':
    main()
