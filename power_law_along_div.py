from extract_JET_data import E_vs_phi
from extract_ITER_data import x_inner, E_inner, phi_inner,\
    x_outer, E_outer, phi_outer
from inventory_phi_E import compute_inv_phi_E
from c_surf_phi_E import phi_H, T_phi_H, r, R_p, D

from scipy.interpolate import interp1d
from scipy.stats import linregress
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import numpy as np

try:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)
except:
    pass

times = [1e5, 4e5, 1.2e6, 3e6, 1e7]
# times = np.logspace(5, 7, num=10, endpoint=True)
distances = [2.54e-2, 0.2, 0.4, 0.8]


def scientificNotation(value):
    if value == 0:
        return '0'
    else:
        e = np.log10(np.abs(value))
        m = np.sign(value) * 10 ** (e - int(e))
        return r'${:.0f} \times 10^{{{:d}}}$'.format(m, int(e))


if __name__ == "__main__":
    invs_inner, invs_outer = [], []
    for t in times:
        inv_phi_E = compute_inv_phi_E(t)
        for x, phi, E, invs in zip(
                [x_inner, x_outer],
                [phi_inner, phi_outer],
                [E_inner, E_outer],
                [invs_inner, invs_outer]):
            inv_along_div = []
            for i in range(len(E)):
                inv_local = inv_phi_E(phi[i], E[i])
                inv_along_div.append(float(inv_local))
            invs.append(interp1d(x, inv_along_div))

    # plotting inventory vs time
    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for i in range(len(distances)):
        inv_d = []
        for j in range(0, len(times)):
            inv_d.append(invs_inner[j](distances[i]))

        slope, intercept, r_value, p_value, std_err = \
            linregress(np.log10(times), np.log10(inv_d))
        plt.plot(times, 10**intercept*np.array(times)**slope, linestyle="--")
        plt.scatter(
            times, inv_d,
            label=r"$x = {:.2f}$ m".format(distances[i], 10**intercept, slope),
            zorder=3)
        txt = scientificNotation(10**intercept) + r"$t^{" + r"{:.2f}".format(slope) + r"}$"
        plt.annotate(txt, xy=(1.1e7, 10**intercept*1.1e7**slope), color=cycle[i], fontsize=15)
    plt.xlabel(r"t (s)", fontsize=20)
    plt.ylabel(r"Inventory per monoblock (H)", fontsize=20)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(right=5e7)
    plt.ylim(top=2e20)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.minorticks_on()
    plt.grid(which='minor', alpha=0.3)
    plt.grid(which='major', alpha=0.7)
    plt.legend(fontsize=18)
    plt.show()

    # plotting b along divertor
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row')
    b_inner, b_outer = [], []
    a_inner, a_outer = [], []
    for x, a, b, invs in zip(
            [x_inner, x_outer], [a_inner, a_outer],
            [b_inner, b_outer], [invs_inner, invs_outer]):
        for d in x:
            inv_d = []
            for j in range(0, len(times)):
                inv_d.append(invs[j](d))

            slope, intercept, r_value, p_value, std_err = \
                linregress(np.log10(times), np.log10(inv_d))
            b.append(slope)
            a.append(10**intercept)
    for x, a, b, i in zip(
            [x_inner, x_outer], [a_inner, a_outer],
            [b_inner, b_outer], [0, 1]):
        coeffs = [a, b]

        for ax, coeff in zip(axs, coeffs):
            ax[i].scatter(
                x, coeff, zorder=3,
                facecolors=(*to_rgb('tab:blue'), 0.1),
                edgecolors=(*to_rgb('tab:blue'), 1))
            ax[i].minorticks_on()
            ax[i].grid(which='minor', alpha=0.3)
            ax[i].grid(which='major', alpha=0.7)

    axs[0][0].set_yscale("log")
    axs[0][1].set_yscale("log")
    axs[1][0].set_xlabel(r"Distance along inner target (m)")
    axs[1][1].set_xlabel(r"Distance along outer target (m)")
    axs[0][0].set_ylabel(r"$a$ (H s $^{-b}$)")
    axs[1][0].set_ylabel(r"$b$")

    plt.show()

    # plot b vs a
    plt.scatter(a_inner, b_inner)
    plt.scatter(a_outer, b_outer)
    plt.xscale("log")
    plt.show()

    # plotting b vs T_surf
    T_inner, T_outer = [], []
    for E, phi, b, T_along_div in zip([E_inner, E_outer], [phi_inner, phi_outer], [b_inner, b_outer], [T_inner, T_outer]):
        for i in range(len(E)):
            T_local = T_phi_H(phi_H(phi[i], E[i]))
            T_along_div.append(T_local)
        plt.scatter(
            T_along_div, b, zorder=3,
            facecolors=(*to_rgb('tab:blue'), 0.1),
            edgecolors=(*to_rgb('tab:blue'), 1))
    plt.minorticks_on()
    plt.grid(which='minor', alpha=0.3)
    plt.grid(which='major', alpha=0.7)
    plt.xlabel(r"$T_\mathrm{surface}$ (K)")
    plt.ylabel(r"$b$ coefficient")
    plt.show()

    # plot b vs T, c
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    c_max_inner = R_p(np.array(E_inner))*phi_inner*(1-r(np.array(E_inner)))/D(np.array(T_inner))
    ax.scatter(T_inner, np.log10(c_max_inner), b_inner)
    c_max_outer = R_p(np.array(E_outer))*phi_outer*(1-r(np.array(E_outer)))/D(np.array(T_outer))
    ax.scatter(T_outer, np.log10(c_max_outer), b_outer)
    ax.set_ylabel(r"$c_\mathrm{surface}$")
    ax.set_xlabel(r"$T_\mathrm{surface}$")
    ax.set_zlabel(r"$b$ coefficient")
    plt.show()

    # plot b vs c
    plt.scatter(
        c_max_inner, b_inner, zorder=3,
        facecolors=(*to_rgb('tab:blue'), 0.1),
        edgecolors=(*to_rgb('tab:blue'), 1))
    plt.scatter(
        c_max_outer, b_outer, zorder=3,
        facecolors=(*to_rgb('tab:blue'), 0.1),
        edgecolors=(*to_rgb('tab:blue'), 1))
    plt.xscale("log")
    plt.minorticks_on()
    plt.grid(which='minor', alpha=0.3)
    plt.grid(which='major', alpha=0.7)
    plt.xlabel(r"$c_\mathrm{surface}$")
    plt.ylabel(r"$b$ coefficient")
    plt.show()
