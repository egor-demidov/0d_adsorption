import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import medfilt, butter, filtfilt
from scipy.signal import savgol_filter, find_peaks
import pandas as pd
import openpyxl as xl
import json


file_name = 'uptake_curve_processing/uptake_curves.xlsx'
sheet_name = 'levoglucosan'
# median_filter_kernel = 51
reactor_diameter = 1.56 # cm
reactor_cross_section = np.pi * (reactor_diameter / 2.0) ** 2.0
concentration_correction_factor = 1.0 / 8.0
out_name = f'uptake_curve_processing/{sheet_name}/drift_corrected.json'

default_initial_guess = {
    'k_ads': 2.628413090171913e-12,
    'k_des': 0.02358734401037852,
    'k_rxn': 1.3703757029773324e-17,
    'S_tot': 58264568110461.61,
    'Y_tot': 96507325673713.16
}

with open(file_name, 'rb') as file:
    data = pd.read_excel(file, sheet_name=sheet_name, usecols='A,B')

xl_wb = xl.load_workbook(file_name)
xl_ws = xl_wb[sheet_name]

# Normalize data
data.signal /= data.signal[0]
# Convert minutes to seconds
data.time *= 60.0

def lowpass_filtfilt(y, fs, fc, order=3):
    # fs: sampling frequency (Hz), fc: cutoff frequency (Hz)
    # b, a = butter(order, fc/(0.5*fs), btype='low')
    b, a = butter(order, fc, btype="low", fs=fs)
    return filtfilt(b, a, y)

# Apply median filter to smooth the signal
dt = np.median(np.diff(data.time))   # seconds/sample
fs = 1.0 / dt                # samples/second (Hz)
filtered_signal = lowpass_filtfilt(medfilt(data.signal, kernel_size=9), fs=fs, fc=fs/10.0, order=3)
# data.signal = medfilt(data.signal, kernel_size=median_filter_kernel)

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
            'X_in': concentration,
            'pressure': pressure,
            'Di': diffusion_coeff,
            't_ads_start': t_ads_start,
            't_ads_end': t_ads_end,
            'initial_guess': default_initial_guess,
            'experimental_data': {
                't_exp': data.time.tolist(),
                'X_exp': (data.signal * concentration).tolist()
            }
        }

        with open(out_name, 'w') as f:
            json.dump(data_to_write, f, indent=4, sort_keys=False, ensure_ascii=False)


if __name__ == '__main__':
    DriftCorrection()
