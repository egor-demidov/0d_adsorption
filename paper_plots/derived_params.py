import json
import string
import matplotlib.pyplot as plt
from uncertainties import ufloat
import uncertainties.unumpy as unp
import uncertainties.umath as um
import math
import numpy as np
from pathlib import Path
from scipy.optimize import fsolve

plt.rcParams.update({
    "font.size": 12,        # default text size
})

R =  0.78       # cm
T = 300         # K
M = 271.52e-3   # kg/mol
Rg = 8.314      # J/mol/K
omega = np.sqrt(Rg * T / 2.0 / np.pi / M) * 1.0e2   # cm/s

nacl = {
    "name": "Nacl",
    "input_path": Path('figure_3/nacl_drift_corrected.json'),
    "fitted_path": Path('figure_3/nacl_fitted.json'),
    "gamma_qss_t": 539.0,
    "C_ml": 1.0 / (5.59e-8) ** 2.0      # surface concentration of units cells
}

levoglucosan = {
    "name": "levoglucosan",
    "input_path": Path('figure_3/levoglucosan_drift_corrected.json'),
    "fitted_path": Path('figure_3/levoglucosan_fitted.json'),
    "gamma_qss_t": 280.0,
    "C_ml": 1.0/0.4e-14                   # surface concentration of molecules, 0.4 nm^2 is projected area of a molecule
}

compounds = [nacl, levoglucosan]

fig, axes = plt.subplots(1, 2, figsize=(10, 4.3))

for compound, ax in zip(compounds, axes):

    # Load the fitted parameters and solution
    with open(compound["fitted_path"], "r") as f:
        data = json.load(f)

    # Load input parameters
    with open(compound["input_path"], "r") as f:
        input_data = json.load(f)

    k_ads = ufloat(data["solution"]["k_ads"], data["standard_error"]["k_ads"])
    k_des = ufloat(data["solution"]["k_des"], data["standard_error"]["k_des"])
    k_rxn = ufloat(data["solution"]["k_rxn"], data["standard_error"]["k_rxn"])
    S_tot = ufloat(data["solution"]["S_tot"], data["standard_error"]["S_tot"])
    Y_tot = ufloat(data["solution"]["Y_tot"], data["standard_error"]["Y_tot"])

    ts = np.array(input_data["experimental_data"]["t_exp"])
    t_ads_start = np.array(input_data["t_ads_start"])
    t_ads_end = np.array(input_data["t_ads_end"])
    k_ads_smooth = np.array(input_data["k_ads_smooth"])
    Xgs = np.array(data["fitted_data"]["Xgs"])
    Xs = np.array(data["fitted_data"]["Xs"])
    P = np.array(data["fitted_data"]["P"])
    r_uptake = np.array(data["fitted_data"]["uptake_rate"])
    gamma = r_uptake / (Xgs * omega)
    r_uptake_qss = k_ads*Xgs*S_tot * (1.0 - 1.0 / (1.0 + k_rxn * (Y_tot - P) / (k_ads*Xgs + R/2.0*k_des)))
    gamma_qss = r_uptake_qss / (Xgs * omega)
    gamma_qss_nominal = unp.nominal_values(gamma_qss)
    dXs_dt = r_uptake - k_rxn.n * Xs * (Y_tot.n - P)

    qss_idx = np.abs(ts - compound["gamma_qss_t"]).argmin()

    print(f'{compound["name"]} FITTED PARAMETERS')

    print(f"k_ads\t\t{k_ads:.2uE}")
    print(f"k_des\t\t{k_des:.2uE}")
    print(f"k_rxn\t\t{k_rxn:.2uE}")
    print(f"S_tot\t\t{S_tot:.2uE}")
    print(f"Y_tot\t\t{Y_tot:.2uE}")

    print(f'{compound["name"]} DERIVED PARAMETERS')

    R_2_k_des = k_des * R / 2.0
    K_ads = k_ads / R_2_k_des
    K_sa = K_ads * S_tot
    gamma_0 = k_ads * S_tot / omega
    chi_S = S_tot / compound["C_ml"]
    chi_Y = Y_tot / compound["C_ml"]

    print(f"R/2*k_des\t{R_2_k_des:.2uE}")
    print(f"K_ads\t\t{K_ads:.2uE}")
    print(f"K_sa\t\t{K_sa:.2uE}")
    print(f"gamma_0\t\t{gamma_0:.2uE}")
    print(f"gamma_qss\t{gamma_qss[qss_idx]:.2uE}")
    print(f"C_ml\t\t{compound["C_ml"]:.2e}")
    print(f"chi_S\t\t{chi_S:.2uE}")
    print(f"chi_Y\t\t{chi_Y:.2uE}")

    # ax2 = ax.twinx()
    ax.plot(ts, gamma, '-k')
    ax.plot(ts, gamma_qss_nominal, '--', color='tab:blue')
    ax.axhline(y=gamma_0.n, linestyle=':', color='tab:red', linewidth=2)
    ax.text(50.0, gamma_0.n - 0.001, f"$\\gamma_0={gamma_0.n:.03f}$", color="tab:red", va="top", ha="left")
    ax.plot(ts[qss_idx], gamma_qss_nominal[qss_idx], 'o', color='tab:blue')
    exponent = math.floor(math.log10(abs(gamma_qss_nominal[qss_idx])))
    mantissa = gamma_qss_nominal[qss_idx] / 10**exponent
    ax.text(ts[qss_idx]+10, gamma_qss_nominal[qss_idx] + 0.001, "$\\gamma_{\\rm qss}=" + f"{mantissa:.01f}\\times 10^" + "{" + f"{exponent}" + "}$", color="tab:blue", va="bottom", ha="left")
    # ax2.plot(range(len(dXs_dt)), dXs_dt)
    # ax.set_title(compound["name"])

    ylim_old = ax.get_ylim()
    xlim_old = ax.get_xlim()

    dXs_dt_abs = np.abs(dXs_dt)
    shift = 0.8
    gradient_fun = 1.0 - dXs_dt_abs / np.max(dXs_dt_abs) - shift
    gradient_fun = np.clip(gradient_fun, 0, None)

    # Create a mask for gradient
    mask = (ts > t_ads_start - k_ads_smooth) & (ts < t_ads_end + k_ads_smooth)
    gradient_fun *= mask

    # make 2D gradient image
    gradient = np.tile(gradient_fun, (10, 1))

    ax.imshow(
        gradient,
        extent=[ts[0], ts[-1], ylim_old[0], ylim_old[1]*1.1],
        aspect='auto',
        cmap='Purples',
        alpha=0.3,
        origin='lower',
        zorder=0
    )

    ax.set_ylim((None, ylim_old[1]*1.1))

    # get primary axis limits
    # y1min, y1max = ax.get_ylim()

    # compute zero position on primary axis
    # zero_frac = (0 - y1min) / (y1max - y1min)

    # current limits on secondary
    # y2min, y2max = ax2.get_ylim()

    # adjust secondary limits so zero aligns
    # new_span = max(y2max, -y2min)
    # ax2.set_ylim(-new_span * zero_frac/(1-zero_frac), new_span)

    ax.set_xlabel(R'Time, $\rm s$')
    ax.set_ylabel(R'Uptake coefficient')


for ax, label in zip(axes, string.ascii_lowercase):
    ax.text(
        0.04, 0.96,                # position (x,y) in axes coords
        f'({label})',
        transform=ax.transAxes,    # use axes coordinates (0–1)
        fontsize=13,
        fontweight='bold',
        va='top',
        ha='left'
    )

plt.tight_layout()
# plt.savefig('uptake_coefficients.pdf')
plt.show()
