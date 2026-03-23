#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt, butter, filtfilt
from scipy.signal import savgol_filter, find_peaks
import pandas as pd
import openpyxl as xl
import simplejson as json
from decimal import Decimal
import argparse
from pathlib import Path
import string

parser = argparse.ArgumentParser(description="Preprocessing of uptake curves")

parser.add_argument("input", help="Path to excel workbook")
parser.add_argument("--worksheet", required=True, help="Name of the worksheet containing the uptake curve")
parser.add_argument("--output", default="drift_corrected.json", help="Name of output file")

args = parser.parse_args()

file_name = Path(args.input)
sheet_name = args.worksheet
# median_filter_kernel = 51
out_name = Path(args.output)

default_k_ads_smooth = 2.0
default_n_reactors = 10
default_initial_guess = {
    'k_ads': 2.628413090171913e-12,
    'k_des': 0.009199064,
    'k_rxn': 1.3703757029773324e-17,
    'S_tot': 58264568110461.61,
    'Y_tot': 96507325673713.16
}

# Formats can be prescribed for select fields
formats = {
    "k_ads": "{:.6e}",
    "k_des": "{:.6e}",
    "k_rxn": "{:.6e}",
    "S_tot": "{:.6e}",
    "Y_tot": "{:.6e}",
}

def format_value(key, value):
    if key in formats:
        return Decimal(formats[key].format(value))
    return value

def to_jsonable(obj, key=None):
    """Convert pandas/numpy objects to JSON-serializable Python objects,
    and apply per-key numeric formatting as Decimal (stays a JSON number).
    """
    # --- pandas containers ---
    if isinstance(obj, pd.Series):
        # choose ONE:
        return obj.to_dict()          # label->value mapping
        # return obj.tolist()         # list of values (no index)
    if isinstance(obj, pd.DataFrame):
        return obj.to_dict(orient="records")  # list[dict] rows

    # --- numpy scalars/arrays ---
    if isinstance(obj, np.generic):
        obj = obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()

    # --- dict/list recursion ---
    if isinstance(obj, dict):
        return {k: to_jsonable(v, key=k) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v, key=key) for v in obj]

    # per-field numeric formatting (keeps JSON numbers)
    if isinstance(obj, (int, float)) and key in formats:
        return Decimal(formats[key].format(obj))

    return obj

with open(file_name, 'rb') as file:
    data = pd.read_excel(file, sheet_name=sheet_name, usecols='A,B')

xl_wb = xl.load_workbook(file_name)
xl_ws = xl_wb[sheet_name]

reactor_diameter = xl_ws['L2'].internal_value # cm
reactor_cross_section = np.pi * (reactor_diameter / 2.0) ** 2.0

concentration_correction_factor = xl_ws['M2'].internal_value

# Normalize data
data.signal /= data.signal[0]
# Convert minutes to seconds
data.time *= 60.0

def lowpass_filtfilt(y, fs, fc, order=3):
    # fs: sampling frequency (Hz), fc: cutoff frequency (Hz)
    # b, a = butter(order, fc/(0.5*fs), btype='low')
    b, a = butter(order, fc, btype="low", fs=fs)
    return filtfilt(b, a, y)

dt = np.median(np.diff(data.time))   # seconds/sample
fs = 1.0 / dt                # samples/second (Hz)
filtered_signal = lowpass_filtfilt(medfilt(data.signal, kernel_size=9), fs=fs, fc=fs/10.0, order=3)

# Determine indices of extrema
max_idx, _ = find_peaks(filtered_signal, plateau_size=True)
min_idx, _ = find_peaks(-filtered_signal, plateau_size=True)

imin_global = np.flatnonzero(filtered_signal == filtered_signal.min())[0]
imax_global = np.flatnonzero(filtered_signal == filtered_signal.max())[0]

imax_left = max_idx[np.searchsorted(max_idx, imin_global, side="right") - 1]
imin_left = min_idx[np.searchsorted(min_idx, imax_global, side="right") - 1]

t_ads_start = (0.3 * data.time[imax_left] + 0.7 * data.time[imin_global])
t_ads_end = (0.3 * data.time[imin_left] + 0.7 * data.time[imax_global])
print(f't_ads_start: {t_ads_start}')
print(f't_ads_end: {t_ads_end}')

### PLOT FOR PEAK DETECTION DESCRIPTION ###
# plt.rcParams.update({
#     "font.size": 12
# })
# fig, axs = plt.subplots(2, 2, figsize=(10, 4.3*1.8))
# (ax1, ax2), (ax3, ax4) = axs
# ax1.plot(data.time, data.signal)
# ax2.plot(data.time, filtered_signal)
# ax2.plot(data.time[imax_left], filtered_signal[imax_left], 'o', color='tab:orange', markersize=8)
# ax2.annotate("Local maximum", (data.time[imax_left], filtered_signal[imax_left]), xytext=(8, 9), textcoords="offset points", ha='right')
# ax2.plot(data.time[imin_global], filtered_signal[imin_global], 's', color='tab:orange', markersize=8)
# ax2.annotate("Global minimum", (data.time[imin_global], filtered_signal[imin_global]), xytext=(5, 5), textcoords="offset points")
# ax2.plot(data.time[imin_left], filtered_signal[imin_left], 'o', color='tab:green', markersize=8)
# ax2.annotate("Local minimum", (data.time[imin_left], filtered_signal[imin_left]), xytext=(5, -5), textcoords="offset points", va='top')
# ax2.plot(data.time[imax_global], filtered_signal[imax_global], 's', color='tab:green', markersize=8)
# ax2.annotate("Global maximum", (data.time[imax_global], filtered_signal[imax_global]), xytext=(-5, -5), textcoords="offset points", ha='right', va='top')
#
# baseline_pts = [data.time[0], data.time[len(data.time)-121]]
# baseline_pts_y = [data.signal[0], data.signal[len(data.time)-121]]
# ax3.plot(data.time, data.signal)
# ax3.plot([data.time[0], data.time[len(data.time) - 1]], [1.0, 1.0], '--k')
# ax3.plot(baseline_pts, baseline_pts_y, '--r')
# ax3.plot([baseline_pts[-1], data.time[len(data.time)-1]], [baseline_pts_y[-1], baseline_pts_y[-1] ], '--r')
# ax3.plot(baseline_pts, baseline_pts_y, 'or')
#
# max_error = data.signal[len(data.time)-121] - 1
# x1 = baseline_pts[1]
# x0 = baseline_pts[0]
# max_distance = x1 - x0
# def correction_function(x):
#     return (np.heaviside(x - x0, 0.5) - np.heaviside(x - x1, 0.5)) * (x - x0) * max_error / (x1 - x0) + np.heaviside(x - x1, 0.5) * max_error
# data.signal -= correction_function(data.time)
# ax4.plot(data.time, data.signal)
# ax4.plot([data.time[0], data.time[len(data.time) - 1]], [1.0, 1.0], '--k')
# ax4.plot(baseline_pts, np.ones(len(baseline_pts)) * data.signal[0], 'or')
#
# for ax in axs.flatten():
#     ax.set_xlabel(R'Time, $\rm s$')
#     ax.set_ylabel(R'Normalized concentration')
#
# for ax, label in zip(axs.flatten(), string.ascii_lowercase):
#     ax.text(
#         0.04, 0.96,                # position (x,y) in axes coords
#         f'({label})',
#         transform=ax.transAxes,    # use axes coordinates (0–1)
#         fontsize=13,
#         fontweight='bold',
#         va='top',
#         ha='left'
#     )
#
# plt.tight_layout()
# plt.savefig('figure_peak_detection.pdf')
# plt.show()
# exit()
### END OF PLOT FOR PEAK DETECTION DESCRIPTION ###

flow_rate = xl_ws['F2'].internal_value + xl_ws['G2'].internal_value # flow rate, sccm
flow_rate /= 60  # flow rate, scc/s
pressure = xl_ws['H2'].internal_value   # pressure, Torr
concentration = xl_ws['I2'].internal_value * concentration_correction_factor # concentration, 1/cm^3 (including the correction factor)
length = xl_ws['J2'].internal_value # reactor length, cm
diffusion_coeff = xl_ws['K2'].internal_value # reactor length, cm^2/s

class DriftCorrection:

    def __init__(self):

        self.first_point_selected = False

        # First baseline point
        self.baseline_pts = [
            data.time[0]
        ]

        self.fig, self.ax = plt.subplots(1, 1)

        self.fig.suptitle('Step 1: drift correction')

        self.ax.plot([data.time[0], data.time[len(data.time) - 1]], [1.0, 1.0], '--k')
        self.ax.plot(data.time, data.signal)
        self.ax.plot(self.baseline_pts, np.ones(len(self.baseline_pts)), 'or')
        self.ax.set_xlabel('time, s')

        self.fig.canvas.mpl_connect('button_press_event', lambda event : self.onclick(event))
        self.fig.canvas.mpl_connect('close_event', lambda event: self.onclose(event))

        plt.show()

    def onclick(self, event):
        max_error = event.ydata - 1.0
        x0 = self.baseline_pts[-1]
        x1 = event.xdata
        max_distance = x1 - x0

        if max_distance < 0:
            return

        def correction_function(x):
            return (np.heaviside(x - x0, 0.5) - np.heaviside(x - x1, 0.5)) * (x - x0) * max_error / (x1 - x0) + np.heaviside(x - x1, 0.5) * max_error

        self.baseline_pts.append(x1)

        # Apply the correction
        data.signal -= correction_function(data.time)

        self.ax.clear()
        self.ax.plot([data.time[0], data.time[len(data.time) - 1]], [1.0, 1.0], '--k')
        self.ax.plot(data.time, data.signal)
        self.ax.plot(self.baseline_pts, np.ones(len(self.baseline_pts)), 'or')
        self.fig.canvas.draw()

    def onclose(self, event):
        data_to_write = {
            'F': flow_rate,
            'R': reactor_diameter / 2.0,
            'L': length,
            'N_reactors': default_n_reactors,
            'X_feed': concentration,
            'pressure': pressure,
            'Di': diffusion_coeff,
            't_ads_start': t_ads_start,
            't_ads_end': t_ads_end,
            'k_ads_smooth': default_k_ads_smooth,
            'initial_guess': default_initial_guess,
            'experimental_data': {
                't_exp': data.time.tolist(),
                'X_exp': (data.signal * concentration).tolist()
            }
        }

        # Apply the formats
        formatted = to_jsonable(data_to_write)

        with open(out_name, 'w') as f:
            json.dump(formatted, f, use_decimal=True, indent=4)
            # json.dump(data_to_write, f, indent=4, sort_keys=False, ensure_ascii=False)


if __name__ == '__main__':
    DriftCorrection()
