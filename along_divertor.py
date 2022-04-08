from inventory_phi_E import compute_inv_phi_E
from c_surf_phi_E import phi_H, T_phi_H
from extract_JET_data import E_vs_phi
from extract_ITER_data import x_inner, E_inner, phi_inner,\
    x_outer, E_outer, phi_outer
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.ticker import FuncFormatter
import numpy as np
from scipy.stats import linregress

try:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)
except:
    pass


def scientificNotation(value):
    if value == 0:
        return '0'
    else:
        e = np.log10(np.abs(value))
        m = np.sign(value) * 10 ** (e - int(e))
        return r'${:.1f} \times 10^{{{:d}}}$'.format(m, int(e))


color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

if __name__ == "__main__":
    fontsize = 20
    labelsize = 14
    print("Compute along divertor")
    # flux along divertor
    plt.figure()
    phi_H_along_div_outer, phi_H_along_div_inner = [], []
    for i in range(len(E_outer)):
        phi_H_along_div_outer.append(phi_H(phi_outer[i], E_outer[i]))
    for i in range(len(E_inner)):
        phi_H_along_div_inner.append(phi_H(phi_inner[i], E_inner[i]))
    plt.scatter(
        x_outer,
        phi_H_along_div_outer,
        zorder=3
        )
    plt.yscale("log")
    plt.ylabel(r"$\varphi_{H}$ (W m$^{-2}$)")
    plt.xlabel(r"Distance along divertor (m)")
    plt.minorticks_on()
    plt.grid(which='minor', alpha=0.3)
    plt.grid(which='major', alpha=0.7)

    # inventory along divertor
    plt.figure()
    invs_inner, invs_outer = [], []
    times = [1e5, 1e6, 1e7]
    times = np.logspace(5, 7, num=8)
    for t in times:
        inv_phi_E = compute_inv_phi_E(t)

        inv_along_div = []
        for i in range(len(E_inner)):
            inv_local = inv_phi_E(phi_inner[i], E_inner[i])
            inv_along_div.append(float(inv_local))
        invs_inner.append(inv_along_div)

        inv_along_div = []
        for i in range(len(E_outer)):
            inv_local = inv_phi_E(phi_outer[i], E_outer[i])
            inv_along_div.append(float(inv_local))
        invs_outer.append(inv_along_div)

    plt.scatter(
        x_outer,
        inv_along_div,
        zorder=3
        )
    plt.yscale("log")
    plt.ylabel(r"Inventory (H m$^{-1}$)")
    plt.xlabel(r"Distance along divertor (m)")
    plt.minorticks_on()
    plt.grid(which='minor', alpha=0.3)
    plt.grid(which='major', alpha=0.7)

    # temperature along divertor
    plt.figure()
    T_along_div_inner, T_along_div_outer = [], []
    for i in range(len(E_inner)):
        T_local = T_phi_H(phi_H(phi_inner[i], E_inner[i]))
        T_along_div_inner.append(T_local)
    for i in range(len(E_outer)):
        T_local = T_phi_H(phi_H(phi_outer[i], E_outer[i]))
        T_along_div_outer.append(T_local)
    plt.scatter(
        x_outer,
        T_along_div_outer,
        zorder=3
        )
    # plt.yscale("log")
    plt.ylabel(r"$T_{surf}$ (K)", fontsize=fontsize)
    plt.xlabel(r"Distance along divertor (m)", fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(which='minor', alpha=0.3)
    plt.grid(which='major', alpha=0.7)

    # plot everything
    fig, axs = plt.subplots(3, 2, sharex='col', sharey='row', figsize=(10, 8))
    (ax1, ax2), (ax3, ax4), (ax5, ax6) = axs
    ax1.scatter(
            x_inner, phi_H_along_div_inner, zorder=3,
            facecolors=(*to_rgb('tab:blue'), 0.3),
            edgecolors=(*to_rgb('tab:blue'), 1))
    ax1.set_yscale("log")
    ax1.set_ylabel(r"$\varphi_{H}$ (W m$^{-2}$)", fontsize=fontsize)
    ax2.scatter(
            x_outer, phi_H_along_div_outer, zorder=3,
            facecolors=(*to_rgb('tab:blue'), 0.3),
            edgecolors=(*to_rgb('tab:blue'), 1))
    ax2.set_yscale("log")
    ax3.scatter(
        x_inner, T_along_div_inner, zorder=3,
        facecolors=(*to_rgb('tab:blue'), 0.3),
        edgecolors=(*to_rgb('tab:blue'), 1))
    ax3.set_ylabel(r"$T_\mathrm{surface}$ (K)", fontsize=fontsize)
    ax4.scatter(
        x_outer, T_along_div_outer, zorder=3,
        facecolors=(*to_rgb('tab:blue'), 0.3),
        edgecolors=(*to_rgb('tab:blue'), 1))

    for inv, color, t in zip(invs_inner, color_cycle, times):
        ax5.scatter(
            x_inner, inv, zorder=3, label=scientificNotation(t) + ' s',
            facecolors=(*to_rgb(color), 0.3),
            edgecolors=(*to_rgb(color), 1))
        ax5.set_yscale("log")
    ax5.set_ylabel(r"Inventory per \\ monoblock (H)", fontsize=fontsize)
    ax5.set_xlabel(r"Distance along inner target (m)", fontsize=fontsize)
    ax5.legend()
    for inv, color, t in zip(invs_outer, color_cycle, times):
        ax6.scatter(
            x_outer, inv, zorder=3, label=scientificNotation(t) + ' s',
            facecolors=(*to_rgb(color), 0.3),
            edgecolors=(*to_rgb(color), 1))
        ax6.set_yscale("log")
    ax6.set_xlabel(r"Distance along outer target (m)", fontsize=fontsize)
    for tup in axs:
        for ax in tup:
            ax.minorticks_on()
            ax.grid(which='minor', alpha=0.3)
            ax.grid(which='major', alpha=0.7)
            ax.tick_params(axis='both', which='major', labelsize=labelsize)

    # integration
    nb_cassettes = 54
    new_inv = []
    # for x, val in zip(x_inner, invs_inner[-1]):
    #     if 0 < x < 0.2:
    #         new_inv.append(val*1.2)
    #     else:
    #         new_inv.append(val)
    total_invs = []
    for i in range(len(times)):
        val = 16*np.trapz(np.array(invs_inner[i])/12e-3, x_inner) + \
              22*np.trapz(np.array(invs_outer[i])/12e-3, x_outer)
        val *= nb_cassettes
        val /= 6.02214076e23
        total_invs.append(val)
    # print(np.trapz(((np.array(invs_inner[-1]) - np.array(new_inv))/np.array(invs_inner[-1])), x_inner))
    slope, intercept, r_value, p_value, std_err = \
        linregress(np.log10(times), np.log10(total_invs))
    plt.figure()
    plt.scatter(times, total_invs, zorder=3)
    txt = scientificNotation(10**intercept) + \
        r"$t^{" + r"{:.2f}".format(slope) + r"}$"
    plt.plot(np.logspace(5, 7, num=100), 10**intercept*np.logspace(5, 7, num=100)**slope, linestyle="--", label=txt)

    plt.xlabel("t (s)", fontsize=20)
    plt.ylabel("Inventory in divertor (g)", fontsize=20)
    plt.xscale("log")
    plt.legend(fontsize=18)
    plt.minorticks_on()
    plt.grid(which='minor', alpha=0.3)
    plt.grid(which='major', alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.show()
