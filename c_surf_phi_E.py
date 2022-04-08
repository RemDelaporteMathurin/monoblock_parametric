import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import ticker, cm, colors

from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from scipy.stats import linregress

from extract_JET_data import E_vs_phi
from extract_ITER_data import x_inner, E_inner, phi_inner,\
    x_outer, E_outer, phi_outer
import properties

try:
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)
except:
    pass


print("Compute c_surface as function of phi_inc and E_inc")

phi_ = np.logspace(19, 24+np.log10(2), endpoint=True, num=300)
E_ = np.logspace(-1+np.log10(5), 2, endpoint=True, num=300)

phi_H_ = np.logspace(4, 7, endpoint=True)
T_ = [
    324, 324, 324, 324, 325, 325.1, 325.6, 325.7, 326, 326.5, 327, 327.6, 328.3, 329.1, 330, 331,
    332, 333, 335.28, 337, 339, 341, 344, 347, 351, 355, 360, 366, 373, 380, 389, 400, 412,
    425, 441, 460, 482, 507, 536, 570, 610, 657, 712, 777, 854, 944, 1051, 1178, 1327, 1503
    ]
T_phi_H = interp1d(phi_H_, T_, fill_value='extrapolate')


def K(T):
    return 3e-25/(T**0.5)*np.exp(-(-2.06)/8.6e-5/T)
    # return 1.3e-17*np.exp(-0.84/8.6e-5/T)
    # return 2.9e-18*np.exp(-1.16/8.6e-5/T)


def D(T):
    return properties.D_0_W*np.exp(-properties.E_D_W/8.6e-5/T)
    # return 2.06e-7*np.exp(-0.28/8.6e-5/T)
    # return 2.9e-7*np.exp(-0.39/8.6e-5/T)


def phi_H(phi_inc, E):
    return phi_inc*1.6e-19*(E + 13.6)*2.2


def R_p(E):
    return 1.3636e-10*E**(6.372e-1)


def r(E):
    return 2e-8*E**2 - 6e-5*E + 0.8096


def c_0(phi, T, K):
    def fun(x):
        # return K(T, x)*x**2 - phi
        return K(T)*x**2 - phi
    # print(fun(0), fun(1e24))
    sol = root_scalar(fun, bracket=(0, 1e30))
    return sol.root


def calculate_c_0(phi, T):
    values = np.zeros(phi.shape)
    for i in range(len(phi)):
        for j in range(len(phi[i])):
            # print(x2, y2)
            val = c_0(phi[i][j], T[i][j], K)
            # print(val)
            values[i][j] = val
    return values


X, Y = np.meshgrid(phi_, E_)
T = T_phi_H(phi_H(X, Y))

c_max_recomb = R_p(Y)*X*(1-r(Y))/D(T) + (X*(1-r(Y))/K(T))**0.5
# c_max_recomb = R_p(Y)*(1-r(E))X/D(T_phi_H(phi_H(X, Y))) + calculate_c_0(X, T_phi_H(phi_H(X, Y)))
c_max_instant = R_p(Y)*X*(1-r(Y))/D(T)

if __name__ == "__main__":
    colorbar = "viridis"
    # contour c_max non instant
    fig, ax = plt.subplots()
    plt.title("Non-instantaneous recombination")
    levels = np.logspace(
        np.log10(np.min(c_max_recomb)),
        np.log10(np.max(c_max_recomb)),
        100)
    levels2 = np.logspace(
        np.log10(np.min(c_max_recomb)),
        np.log10(np.max(c_max_recomb)),
        12)
    locator = ticker.LogLocator(base=10)
    CS = ax.contourf(X, Y, c_max_recomb, levels=levels, locator=locator, cmap=colorbar)
    fig.colorbar(CS, label=r"$c_{max}$ (m$^{-3}$)", ticks=locator)
    CS2 = ax.contour(X, Y, c_max_recomb, levels=levels2, locator=locator, colors="white")
    plt.scatter(
        phi_inner + phi_outer,
        E_inner + E_outer,
        facecolors='none',
        edgecolors='white'
        )
    for c in CS.collections:  # for avoiding white lines in pdf
        c.set_edgecolor("face")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\varphi_\mathrm{inc}}$ (m$^{-2}$ s$^{-1}$)")
    plt.ylabel(r"$E$ (eV)")

    # contour c_max instant
    fig, ax = plt.subplots()
    # plt.title("Instantaneous recombination")
    levels = np.logspace(
        np.log10(np.min(c_max_instant)),
        np.log10(np.max(c_max_instant)),
        100)
    levels2 = np.logspace(
        np.log10(np.min(c_max_instant)),
        np.log10(np.max(c_max_instant)),
        12)
    locator = ticker.LogLocator(base=10)
    CS = ax.contourf(X, Y, c_max_instant, levels=levels, locator=locator, cmap=colorbar)

    fig.colorbar(CS, label=r"$c_\mathrm{max}$ (m$^{-3}$)", ticks=locator)
    CS2 = ax.contour(X, Y, c_max_instant, levels=levels2, locator=locator, colors="white")
    plt.scatter(
        phi_inner + phi_outer,
        E_inner + E_outer,
        facecolors='none',
        edgecolors='white'
        )
    for c in CS.collections:  # for avoiding white lines in pdf
        c.set_edgecolor("face")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\varphi_\mathrm{inc}}$ (m$^{-2}$ s$^{-1}$)")
    plt.ylabel(r"$E$ (eV)")

    # temperature
    fig, ax = plt.subplots()
    levels = np.linspace(323, 3695, 100, endpoint=True)
    levels2 = np.linspace(323, 3695, 10, endpoint=True)
    CS = ax.contourf(X, Y, T, levels=levels, extend="max", cmap=colorbar)

    plt.scatter(
        phi_inner + phi_outer,
        E_inner + E_outer,
        facecolors='none',
        edgecolors='white'
        )

    fig.colorbar(CS,  label=r"$T_\mathrm{surface}$ (K)", ticks=[*levels[0:len(levels):10], levels[-1]])
    CS2 = ax.contour(X, Y, T, levels=levels2, locator=locator, colors="white")
    plt.xlabel(r"$\varphi_\mathrm{inc}}$ (m$^{-2}$ s$^{-1}$)")
    plt.ylabel(r"$E$ (eV)")
    for c in CS.collections:  # for avoiding white lines in pdf
        c.set_edgecolor("face")
    plt.xscale("log")
    plt.yscale("log")

    # heat flux
    fig, ax = plt.subplots()
    flux = phi_H(X, Y)
    levels = np.logspace(
        np.log10(np.min(flux)),
        np.log10(np.max(flux)),
        100)
    levels2 = np.logspace(
        np.log10(np.min(flux)),
        np.log10(np.max(flux)),
        12)
    CS = ax.contourf(X, Y, flux, levels=levels, locator=ticker.LogLocator(), cmap=colorbar)
    CS2 = ax.contour(X, Y, flux, levels=levels2, locator=locator, colors="white")
    plt.scatter(
        phi_inner + phi_outer,
        E_inner + E_outer,
        facecolors='none',
        edgecolors='white'
        )

    fig.colorbar(CS,  label=r"$\varphi_{H}$ (W m$^{-2}$)", ticks=locator)
    plt.xlabel(r"$\varphi_\mathrm{inc}$ (m$^{-2}$ s$^{-1}$)")
    plt.ylabel(r"$E$ (eV)")

    plt.xscale("log")
    plt.yscale("log")
    for c in CS.collections:  # for avoiding white lines in pdf
        c.set_edgecolor("face")
    plt.figure()
    slope, intercept, r_value, p_value, std_err = linregress(phi_H_, T_)
    # plt.plot(phi_H_, T_phi_H(phi_H_))
    plt.plot(phi_H_, slope*phi_H_ + intercept)
    plt.scatter(phi_H_, T_, zorder=3, facecolors=(*colors.to_rgb('tab:blue'),0.3), edgecolors=(*colors.to_rgb('tab:blue'),1))
    plt.ylabel(r"$T_{surface}$ (K)")
    plt.xlabel(r"$\varphi_{H}$ (W m$^{-2}$)")
    plt.xscale("log")
    plt.minorticks_on()
    plt.grid(which='minor', alpha=0.3)
    plt.grid(which='major', alpha=0.7)
    plt.show()
