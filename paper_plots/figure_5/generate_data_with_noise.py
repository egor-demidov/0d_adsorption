from pathlib import Path
from add_noise import add_noise
from run2drift_corrected import run2drit_corrected

SEEDS = [478408758, 150371857, 796271257, 473363765, 369097843, 618748432, 133275616, 905182760, 817813893, 891406174, 525063725, 676640592, 631875005, 545785448, 904289640, 512087903, 801337777, 705770982, 189073802, 251646239, 695487012, 599308486, 507774040, 587610539, 845064615, 447940191, 999946091, 487862101, 912742488, 509315279, 435495860, 523568120, 399050268, 958714249, 672312241, 422252380, 362346973, 817565413, 937282650, 988053920, 642310551, 470556606, 450960559, 683462482, 828662228, 334123680, 197777833, 445193135, 334837196, 942338842, 250738185, 997557201, 289616562, 918386503, 762955910, 714801992, 335904720, 257811223, 180538123, 812847299, 315679312, 649019009, 848517079, 347683251, 351251291, 598313295, 662430543, 922027346, 701501527, 305423292, 957138254, 732214502, 895770776, 871916578, 487293820, 778542711, 861810377, 290356731, 601086845, 740630032, 563338822, 671704521, 497048347, 501233556, 919066435, 840052161, 407524024, 402850081, 843065430, 205868134, 588073611, 758004202, 584318109, 620829182, 648893371, 441145588, 764837740, 899840817, 919204067, 341874745]

NOISE_MAGNITUDE = 5.0e8

DURATIONS = [
    100, 200, 300, 400, 500, 600
]

PREFIX = Path('combo_1')

def main():
    # Create necessary directories
    for duration in DURATIONS:
        directory = PREFIX / f'run_{duration}'
        directory.mkdir(exist_ok=True)

    # Generate noisy curves
    for duration in DURATIONS:
        noisy_in_path = PREFIX / f'run_{duration}.json'
        for i, seed in zip(range(len(SEEDS)), SEEDS):
            noisy_out_path = PREFIX / f'run_{duration}/run_{duration}_noisy_{i+1}.json'
            add_noise(NOISE_MAGNITUDE, seed, noisy_in_path, noisy_out_path)

    # Convert to 0d_adsorption input
    for duration in DURATIONS:
        parent_dir = PREFIX / f'run_{duration}'
        for i in range(len(SEEDS)):

            # create a directory for a run
            run_dir = parent_dir / f'run_{duration}_noisy_{i+1}'
            run_dir.mkdir(exist_ok=True)

            path_in = parent_dir / f'run_{duration}_noisy_{i+1}.json'
            path_out = run_dir / 'drift_corrected.json'

            # Create an input
            run2drit_corrected(duration, path_in, path_out)

if __name__ == '__main__':
    main()

