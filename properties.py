import matplotlib.pyplot as plt
import numpy as np

# try:
#     plt.rc('text', usetex=True)
#     plt.rc('font', family='serif', size=12)
# except:
#     pass

k_B = 8.6e-5


def polynomial(coeffs, x, main):
    val = coeffs[0]
    for i in range(1, 4):
        if main:
            val += coeffs[i]*np.float_power(x, i)
        else:
            val += coeffs[i]*x**i
    return val


def rhoCp_W(T, main=False):
    coeffs = [2.48160e6, 5.98312e2, -8.30703e-2, 5.15356e-6]
    return polynomial(coeffs, T, main=main)


def thermal_cond_W(T, main=False):
    coeffs = [1.75214e2, -1.07335e-1, 5.03006e-5, -7.84154e-9]
    return polynomial(coeffs, T, main=main)


def rhoCp_Cu(T, main=False):
    coeffs = [3.45899e6, 4.67353e2, 6.14079e-2, 1.68402e-4]
    return polynomial(coeffs, T, main=main)


def thermal_cond_Cu(T, main=False):
    coeffs = [4.02301e+02, -7.88669e-02, 3.76147e-05, -3.93153e-08]
    return polynomial(coeffs, T, main=main)


def rhoCp_CuCrZr(T, main=False):
    coeffs = [3.46007e6, 6.22091e2, 1.51383e-1, -1.79134e-4]
    return polynomial(coeffs, T, main=main)


def thermal_cond_CuCrZr(T, main=False):
    coeffs = [3.12969e2, 2.57678e-01, -6.45110e-4, 5.25780e-7]
    return polynomial(coeffs, T, main=main)


D_0_W = 1.9e-7
E_D_W = 0.2


def D_W(T, main=False):
    return D_0_W*np.exp(-E_D_W/k_B/T)


D_0_Cu = 6.6e-7
E_D_Cu = 0.387


def D_Cu(T, main=False):
    return D_0_Cu*np.exp(-E_D_Cu/k_B/T)


D_0_CuCrZr = 3.92e-7
E_D_CuCrZr = 0.418


def D_CuCrZr(T, main=False):
    return D_0_CuCrZr*np.exp(-E_D_CuCrZr/k_B/T)


if __name__ == "__main__":
    def tick_function(X):
        V = 1000/(X)
        return ["%.0f" % z for z in V]

    T = np.arange(300, 1300, step=1)
    T_Cu = np.arange(470, 1200, step=1)
    T_CuCrZr = np.arange(561, 769, step=1)
    # T_W = np.arange()
    grey_W = "tab:grey"
    orange_Cu = (228/255, 146/255, 64/255)
    yellow_CuCrZr = (180/255, 95/255, 6/255)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.plot(1000/T, D_W(T, main=True), label="W", color=grey_W)
    ax1.plot(1000/T_Cu, D_Cu(T_Cu, main=True), label="Cu", color=orange_Cu)
    ax1.plot(1000/T_CuCrZr, D_CuCrZr(T_CuCrZr, main=True), label="CuCrZr", color=yellow_CuCrZr)

    ax1.annotate(
        "W", (1000/T[-1], D_W(T[-1], main=True)),
        fontsize=12,
        color=grey_W)
    ax1.annotate(
        "Cu", (1000/T_Cu[-1], D_Cu(T_Cu[-1], main=True)),
        fontsize=12,
        color=orange_Cu)
    ax1.annotate(
        "CuCrZr", (1000/T_CuCrZr[-1], D_CuCrZr(T_CuCrZr[-1], main=True)),
        fontsize=12,
        color=yellow_CuCrZr)

    new_tick_locations = np.array([1, 1.5, 2, 2.5, 3])
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(tick_function(new_tick_locations))
    plt.yscale("log")
    ax1.minorticks_on()
    ax1.grid(which='minor', alpha=0.3)
    ax1.grid(which='major', alpha=0.7)
    ax1.set_ylabel(r"$D$ (m$^{2}$ s$^{-1}$)", fontsize=15)
    ax1.set_xlabel(r"1000/T (K$^{-1}$)", fontsize=15)
    ax2.set_xlabel(r"T (K)", fontsize=15)
    plt.tight_layout()

    plt.figure()
    ax1 = plt.subplot(211)
    # T_Cu = np.arange(300, 1085+273.15, step=1)
    plt.plot(T, rhoCp_W(T, main=True), label="W", color=grey_W)
    plt.plot(T_Cu, rhoCp_Cu(T_Cu, main=True), label="Cu", color=orange_Cu)
    plt.plot(T_CuCrZr, rhoCp_CuCrZr(T_CuCrZr, main=True), label="CuCrZr", color=yellow_CuCrZr)

    plt.ylabel(r"$\rho \cdot C_p$ (J m$^{-3}$ K$^{-1}$)", fontsize=15)
    # plt.yscale("log")
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.minorticks_on()
    plt.grid(which='minor', alpha=0.3)
    plt.grid(which='major', alpha=0.7)

    # share x and y
    ax2 = plt.subplot(212, sharex=ax1)
    plt.xlabel(r"T (K)", fontsize=15)
    plt.ylabel(r"$\lambda$ (W m$^{-1}$ K$^{-1}$)", fontsize=15)
    plt.plot(T, thermal_cond_W(T, main=True), label="W", color=grey_W)
    plt.plot(T_Cu, thermal_cond_Cu(T_Cu, main=True), label="Cu", color=orange_Cu)
    plt.plot(T_CuCrZr, thermal_cond_CuCrZr(T_CuCrZr, main=True), label="CuCrZr", color=yellow_CuCrZr)

    ax1.annotate("W", (T[-1], rhoCp_W(T[-1], main=True)), color=grey_W)
    ax1.annotate("Cu", (T_Cu[-1], rhoCp_Cu(T_Cu[-1], main=True)), color=orange_Cu)
    ax1.annotate("CuCrZr", (T_CuCrZr[-1], rhoCp_CuCrZr(T_CuCrZr[-1], main=True)), color=yellow_CuCrZr)
    ax2.annotate("W", (T[-1], thermal_cond_W(T[-1], main=True)), color=grey_W)
    ax2.annotate("Cu", (T_Cu[-1], thermal_cond_Cu(T_Cu[-1], main=True)), color=orange_Cu)
    ax2.annotate("CuCrZr", (T_CuCrZr[-1], thermal_cond_CuCrZr(T_CuCrZr[-1], main=True)), color=yellow_CuCrZr)

    ax1.set_ylim(bottom=0, top=5e6)
    ax2.set_ylim(bottom=0, top=400)
    # plt.yscale("log")
    plt.minorticks_on()
    plt.grid(which='minor', alpha=0.3)
    plt.grid(which='major', alpha=0.7)
    plt.show()
