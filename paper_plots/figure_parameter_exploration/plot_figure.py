from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import json
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import matplotlib.patheffects as pe

bkr = LinearSegmentedColormap.from_list(
    "bkr",
    ["blue", "black", "red"]
)

bkr_narrow = LinearSegmentedColormap.from_list(
    "bkr_narrow",
    [
        (0.00, "royalblue"),
        (0.45, "navy"),
        (0.495, "black"),
        (0.500, "black"),
        (0.505, "black"),
        (0.55, "darkred"),
        (1.00, "red"),
    ],
)

plt.rcParams.update({
    "font.size": 12,        # default text size
    # "axes.titlesize": 14,
    # "axes.labelsize": 13,
    # "xtick.labelsize": 11,
    # "ytick.labelsize": 11,
    # "legend.fontsize": 11
})


parameter_map = [
    {
        'name': 'k_ads',
        'dir': Path('k_ads'),
        'label': R'$k_{\rm ads}$',
        'datasets': []
    },
    {
        'name': 'k_des',
        'dir': Path('k_des'),
        'label': R'$k_{\rm des}$',
        'datasets': []
    },
    {
        'name': 'k_rxn',
        'dir': Path('k_rxn'),
        'label': R'$k_{\rm rxn}$',
        'datasets': []
    },
    {
        'name': 'S_tot',
        'dir': Path('S_tot'),
        'label': R'$S_{\rm tot}$',
        'datasets': []
    },
    {
        'name': 'Y_tot',
        'dir': Path('Y_tot'),
        'label': R'$Y_{\rm tot}$',
        'datasets': []
    }
]

for i, param in enumerate(parameter_map):
    for file_path in param['dir'].glob('*.json'):
        if file_path.is_file() and not file_path.name.startswith('run_'):
            print(f'Loading dataset {file_path.stem}...')

            run_file_path = file_path.with_name('run_' + file_path.name)

            with open(run_file_path, 'r') as in_file:
                data = json.load(in_file)

                dataset = {
                    't': data['solution']['t'],
                    'X': data['solution']['X'],
                    'label': file_path.stem
                }

                parameter_map[i]['datasets'].append(dataset)


fig, axs = plt.subplot_mosaic([
    ['A', 'A', 'B', 'B'],
    ['C', 'C', 'D', 'D'],
    ['.', 'E', 'E', '.'],
], figsize=(10, 4.3*2.8))

def plot_datasets(parameter_idx: int, ax: plt.Axes, inset_xs, inset_ys, inset_width, inset_height, inset_loc, hide_labels, reverse_order, label_at_min, show_inset, custom_label_pos):

    ax.set_xlabel('Time, s')
    ax.set_ylabel('X concentration, $\\rm cm^{-3}$')


    ax.text(
        0.05,
        0.95,
        parameter_map[parameter_idx]['label'],
        transform=ax.transAxes,
        ha="left",
        va="top",
        color="black",
        fontsize=12,
    )

    if show_inset:
        axins = inset_axes(ax, width=inset_width, height=inset_height, loc=inset_loc)

    factors = []
    for dataset in parameter_map[parameter_idx]['datasets']:
        if dataset['label'] != 'original':
            factors.append(float(dataset['label']))
        else:
            factors.append(1.0)

    factors = np.log(np.array(factors))
    normalized_factors = factors - np.min(factors)
    normalized_factors = normalized_factors / normalized_factors.max()

    order = np.argsort(normalized_factors)
    if reverse_order:
        order = order[::-1]

    print(normalized_factors[order])

    for n, i in enumerate(order):

        dataset = parameter_map[parameter_idx]['datasets'][i]

        if dataset['label'] == 'original':
            color = 'k'
            percent = '0 %'
        else:
            # cmap = plt.get_cmap('bwr')
            cmap = bkr_narrow
            color = cmap(normalized_factors[i])
            percent = (np.exp(factors[i]) - 1.0) * 100
            if percent > 0:
                percent = f'+{percent:.0f} %'
            else:
                percent = f'{percent:.0f} %'

        ax.fill_between(dataset['t'], dataset['X'], dataset['X'][0], alpha=1.0, color=color)
        if show_inset:
            axins.fill_between(dataset['t'], dataset['X'], dataset['X'][0], alpha=1.0, color=color)


        label_ax = axins if show_inset else ax

        if custom_label_pos is not None:
            if n not in hide_labels:
                next_X = parameter_map[parameter_idx]['datasets'][order[n+1]]['X'][custom_label_pos] if n < len(parameter_map[parameter_idx]['datasets']) - 1 else dataset['X'][0]

                label_ax.text(
                    dataset['t'][custom_label_pos],
                    (dataset['X'][custom_label_pos]+next_X) / 2.0,
                    percent,
                    color="white",
                    fontsize=10,
                    ha="center",
                    va="center",
                    path_effects=[
                        pe.withStroke(linewidth=3, foreground="black")
                    ],
                    )
        else:
            if label_at_min:
                min_loc = np.argmin(dataset['X'])
                next_min = np.min(parameter_map[parameter_idx]['datasets'][order[n+1]]['X']) if n < len(parameter_map[parameter_idx]['datasets']) - 1 else dataset['X'][0]

                if n not in hide_labels:
                    label_ax.text(
                        dataset['t'][min_loc]+1,
                        (dataset['X'][min_loc]+next_min) / 2.0,
                        percent,
                        color="white",
                        fontsize=10,
                        ha="center",
                        va="center",
                        path_effects=[
                            pe.withStroke(linewidth=3, foreground="black")
                        ],
                    )

            else:
                max_loc = np.argmax(dataset['X'])
                next_max = np.max(parameter_map[parameter_idx]['datasets'][order[n+1]]['X']) if n < len(parameter_map[parameter_idx]['datasets']) - 1 else dataset['X'][0]

                if n not in hide_labels:
                    label_ax.text(
                        dataset['t'][max_loc]+1,
                        (dataset['X'][max_loc]+next_max) / 2.0,
                        percent,
                        color="white",
                        fontsize=10,
                        ha="center",
                        va="center",
                        path_effects=[
                            pe.withStroke(linewidth=3, foreground="black")
                        ],
                        )


    if show_inset:
        # Define zoom region
        x1, x2 = inset_xs
        y1, y2 = inset_ys
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)

        # Remove ticks (optional)
        axins.set_xticks([])
        axins.set_yticks([])

        # Draw box and connectors
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")


        # ax.plot(
        #     dataset['t'],
        #     dataset['X'],
        #     color=color,
        #     linewidth=2,
        #     alpha=0.6
        # )

plot_datasets(
    0,
              axs['A'],
              (95, 120),
              (1.1e10, 3e10),
              '65%',
              '60%',
              'lower right',
              (2, 4),
    True,
    True,
    True,
    None
)

plot_datasets(
    1,
    axs['B'],
    (395, 465),
    (2.9e10, 3.6e10),
    '55%',
    '50%',
    'lower right',
    (),
    True,
    False,
    True,
    None
)

plot_datasets(
    2,
    axs['C'],
    (95, 120),
    (1.1e10, 3e10),
    '65%',
    '60%',
    'lower right',
    (),
    True,
    True,
    False,
    250
)

plot_datasets(
    3,
    axs['D'],
    (95, 120),
    (1.1e10, 3e10),
    '65%',
    '60%',
    'lower right',
    (),
    True,
    True,
    False,
    250
)

plot_datasets(
    4,
    axs['E'],
    (95, 120),
    (1.1e10, 3e10),
    '65%',
    '60%',
    'lower right',
    (),
    True,
    True,
    False,
    350
)

fig.tight_layout()
fig.savefig('figure_parameter_exploration.png', dpi=300)
plt.show()
