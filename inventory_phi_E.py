from c_surf_phi_E import X, Y, T, c_max_instant, phi_, E_
from inventory_T_c import interp_extrap_2D, compute_inv_T_c
from extract_JET_data import E_vs_phi
from extract_ITER_data import x_inner, E_inner, phi_inner,\
    x_outer, E_outer, phi_outer

from scipy.interpolate import interp2d

import math

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.patches import Rectangle

try:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)
except:
    pass

print("Compute inventory as function of phi_inc and E_inc")


def calculate_inv(T, c, inv_T_c):
    values = np.zeros(T.shape)
    for i in range(len(T)):
        for j in range(len(T[i])):
            if T[i][j] < 1200:
                val = inv_T_c(T[i][j], c[i][j])
                values[i][j] = val
    return values


def find_min(z):
    ''' finds min of a 2D array'''
    minimum = np.max(z)
    for l in z:
        for e in l:
            if not math.isnan(e) and e > 0:
                if e < minimum:
                    minimum = e
    return minimum


def compute_inv_phi_E(t):
    inv = calculate_inv(T, c_max_instant,
                        compute_inv_T_c(T, c_max_instant, t=t))
    inv_phi_E = interp2d(phi_, E_, inv, kind='cubic')
    return inv_phi_E


if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    fontsize = 20
    labelsize = 14
    inv = calculate_inv(T, c_max_instant,
                        compute_inv_T_c(T, c_max_instant, t=1e7))

    minimum = find_min(inv)
    levels = np.logspace(
        np.log10(1e18),
        np.log10(np.max(inv)),
        100)
    levels2 = np.logspace(
        np.log10(minimum),
        np.log10(np.max(inv)),
        10)
    locator = ticker.LogLocator()
    CS = plt.contourf(X, Y, inv, levels=levels, locator=locator)
    for c in CS.collections:  # for avoiding white lines in pdf
        c.set_edgecolor("face")

    formatter = ticker.LogFormatter(10, labelOnlyBase=False, minor_thresholds=(np.inf, np.inf))
    cbar = plt.colorbar(CS, ticks=locator)
    cbar.formatter = formatter
    CS2 = plt.contour(
        X, Y, inv, levels=levels2, locator=locator,
        colors="white", linewidths=1)
    plt.scatter(
        phi_inner + phi_outer,
        E_inner + E_outer,
        facecolors='none',
        edgecolors='white'
        )
    ax.tick_params(axis='both', which='major', labelsize=labelsize)
    plt.xlabel(r"$\varphi_\mathrm{inc}$ (m$^{-2}$ s$^{-1}$)", fontsize=fontsize)
    plt.ylabel(r"$E$ (eV)", fontsize=fontsize)
    plt.xscale("log")
    plt.yscale("log")

    cbar.ax.tick_params(labelsize=labelsize)
    cbar.set_label(label="Inventory per monoblock (H)", size=fontsize)

    xmin, xmax, ymin, ymax = plt.axis()
    currentAxis = plt.gca()
    currentAxis.add_patch(
        Rectangle(
            (xmin, ymin), xmax-xmin, ymax-ymin, color="grey",
            alpha=0.4, fill=True, zorder=0))
    plt.tight_layout()
    plt.savefig("Figures/inventory_phi_E.pdf")
    plt.show()
