from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import csv

# Read E and phi from Denis et al (NME 2020)
x_inner, E_inner, phi_inner = [], [], []
x_outer, E_outer, phi_outer = [], [], []

with open('exposure_conditions_divertor/ITER/inner_target.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=';')
    next(readCSV)
    for row in readCSV:
        x_inner.append(float(row[0]))
        E_inner.append(float(row[1]))
        phi_inner.append(float(row[2]))
with open('exposure_conditions_divertor/ITER/outer_target.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=';')
    next(readCSV)
    for row in readCSV:
        x_outer.append(float(row[0]))
        E_outer.append(float(row[1]))
        phi_outer.append(float(row[2]))


if __name__ == "__main__":
    # plt.scatter(E_vs_phi[:, 0], E_vs_phi[:, 1])
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row')
    (ax1, ax2), (ax3, ax4) = axs
    ax1.scatter(x_inner, E_inner, zorder=3,
                facecolors=(*to_rgb('tab:blue'), 0.3),
                edgecolors=(*to_rgb('tab:blue'), 1))
    ax1.set_ylabel(r"E (eV)")
    ax2.scatter(x_outer, E_outer,
                facecolors=(*to_rgb('tab:blue'), 0.3),
                edgecolors=(*to_rgb('tab:blue'), 1))
    ax3.scatter(x_inner, phi_inner,
                facecolors=(*to_rgb('tab:blue'), 0.3),
                edgecolors=(*to_rgb('tab:blue'), 1))
    ax3.set_ylabel(r"$\varphi_\mathrm{inc}$ (m$^{-2}$ s$^{-1}$)")
    ax3.set_xlabel(r"Distance along inner target (m)")
    ax4.scatter(x_outer, phi_outer,
                facecolors=(*to_rgb('tab:blue'), 0.3),
                edgecolors=(*to_rgb('tab:blue'), 1))
    ax4.set_xlabel(r"Distance along outer target (m)")
    for tup in axs:
        for ax in tup:
            ax.minorticks_on()
            ax.grid(which='minor', alpha=0.3)
            ax.grid(which='major', alpha=0.7)
            ax.set_yscale("log")
    # plt.xscale("log")
    plt.show()
