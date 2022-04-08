import numpy as np
import scipy as sp
from scipy.interpolate import griddata, interp2d, interp1d
import matplotlib.pyplot as plt
from matplotlib import ticker
from process_T_c_data import points, data
import globalspine


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
        return r'${:.0f} \times 10^{{{:d}}}$'.format(m, int(e))


def inv(points, time=1e7):
    ''' returns a 1D array'''
    values = []
    for p in points:
        for d in data:
            if d["T"] == p[0] and d["c"] == p[1]:
                values.append(interp1d(d["t"], d["inventory"])(time))
                break
    return np.asarray(values)


def compute_inv(X, Y):
    ''' returns a 2D array (meshgrid)'''
    values = np.zeros(X.shape)
    for i in range(len(X)):
        for j in range(len(X[i])):
            # print(x2, y2)
            val = interp_extrap_2D(X[i][j], Y[i][j])
            # val = interp(X[i][j], Y[i][j])
            # print(val)
            values[i][j] = val
    return values


def compute_inv_T_c(T, c, t):
    # grid the data
    grid_z = griddata(
        (points[:, 0], points[:, 1]/1e20),
        inv(points, time=t)/1e20,
        (X, Y), method='cubic', fill_value=1)

    # interpolate the gridded data
    # nearest neighbor extrapolation
    e = 12e-3  # monoblock thickness
    # interp = interp2d(x, y*1e20, e*grid_z*1e20, kind='cubic')

    # better extrapolation
    interp_extrap_2D = \
        globalspine.GlobalSpline2D(x, y*1e20, e*grid_z*1e20, kind='cubic')
    return interp_extrap_2D


print("Compute inventory as function of T and c_surface")


# create the meshgrid
scaling_factor = 1e20
x, y = np.linspace(330, 1200, num=9), \
    np.logspace(20, 22+np.log10(10), num=9, endpoint=True)/scaling_factor
X, Y = np.meshgrid(x, y)

# grid the data
grid_z = griddata(
    (points[:, 0], points[:, 1]/1e20),
    inv(points, time=1e7)/1e20,
    (X, Y), method='cubic', fill_value=1)

# interpolate the gridded data
# nearest neighbor extrapolation
e = 12e-3  # monoblock thickness
# interp = interp2d(x, y*1e20, e*grid_z*1e20, kind='cubic')

# better extrapolation
interp_extrap_2D = \
    globalspine.GlobalSpline2D(x, y*scaling_factor, e*grid_z*scaling_factor, kind='cubic')
if __name__ == "__main__":
    # with inference-tools
    points2 = []
    z = []
    for p in points:
        if 320 <= p[0] <= 1100 and 1e20 <= p[1] <= 1e23:
            points2.append([p[0], np.log10(p[1])])
            z.append(np.log10(e*inv([p], time=1e7)))
    points2 = np.array(points2)

    # scatter the dots
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter([i[0] for i in points2], [i[1] for i in points2], z)
    plt.tight_layout()
    plt.show()


    # Train the GP on the data
    from inference.gp_tools import GpRegressor
    from inference.gp_tools import RationalQuadratic, SquaredExponential
    GP = GpRegressor(points2, z, kernel=RationalQuadratic)

    # evaluate the estimate
    N = 100
    gp_x = np.linspace(320, 1100, N)
    gp_y = np.log10(np.logspace(20, 23, N))
    gp_coords = [(i, j) for i in gp_x for j in gp_y]
    mu, sig = GP(gp_coords)

    fig = plt.figure()
    locator = ticker.LogLocator(base=10)
    levels = np.logspace(
        min(mu),
        max(mu),
        1000)
    levels2 = np.logspace(
        min(mu),
        max(mu),
        8)

    CS = plt.contourf(*np.meshgrid(gp_x, 10**gp_y), 10**mu.reshape([N, N]).T, locator=locator, levels=levels)
    plt.colorbar(CS, label=r"Inventory per monoblock (H)", ticks=locator)
    CS2 = plt.contour(
        *np.meshgrid(gp_x, 10**gp_y), 10**mu.reshape([N, N]).T, levels=levels2,
        locator=locator, colors="white")
    CLS = plt.clabel(CS2, inline=True, fontsize=10, fmt=scientificNotation)

    # overplot the data points
    plt.scatter(points2[:, 0], 10**points2[:, 1], color='black')
    plt.yscale("log")
    plt.show()


    fig = plt.figure()
    levels = np.linspace(
        min(sig),
        max(sig),
        1000)

    X2, Y2 = np.meshgrid(
        np.linspace(320, 1100, num=100),
        np.logspace(20, 22+np.log10(7), num=100, endpoint=True))
    CS = plt.contourf(*np.meshgrid(gp_x, 10**gp_y), 10**sig.reshape([N, N]).T, levels=1000)
    plt.colorbar(CS, label=r"$\sigma$")

    # overplot the data points
    plt.scatter(points2[:, 0], 10**points2[:, 1], color='black')
    plt.yscale("log")
    plt.show()

    # plot profile along T axis
    n = 0
    for c in gp_y:
        n += 1
        if  3e21 < 10**c < 4e21:
            break
    plt.figure()
    mu_reshaped = mu.reshape([N, N]).T
    sig_reshaped = sig.reshape([N, N]).T
    plt.plot(gp_x, 10**mu_reshaped[n], label='posterior mean c_surf = {:.1e} m-3'.format(10**gp_y[n]))
    plt.fill_between(
        gp_x, 10**(mu_reshaped[n]-sig_reshaped[n]), 10**(mu_reshaped[n]-sig_reshaped[n]*2),
        facecolor='blue', alpha=0.15, label=r'$\pm 2 \sigma$ interval')
    plt.fill_between(
        gp_x, 10**(mu_reshaped[n]+sig_reshaped[n]), 10**(mu_reshaped[n]+sig_reshaped[n]*2),
        facecolor='blue', alpha=0.15)
    plt.fill_between(
        gp_x, 10**(mu_reshaped[n]-sig_reshaped[n]), 10**(mu_reshaped[n]+sig_reshaped[n]),
        facecolor='blue', alpha=0.3, label=r'$\pm 1 \sigma$ interval')
    plt.legend()
    plt.yscale("log")
    plt.show()


    # plot profile along c axis
    n = 0
    for T in gp_x:
        n += 1
        if  840 < T < 850:
            break
    plt.figure()
    mu_reshaped = mu.reshape([N, N]).T
    sig_reshaped = sig.reshape([N, N]).T
    plt.plot(gp_y, 10**mu_reshaped[n], label='posterior mean T_surf = {:.0f} K'.format(gp_x[n]))
    plt.fill_between(
        gp_y, 10**(mu_reshaped[:, n]-sig_reshaped[:, n]),
        10**(mu_reshaped[:, n]-sig_reshaped[:, n]*2),
        facecolor='blue', alpha=0.15, label=r'$\pm 2 \sigma$ interval')
    plt.fill_between(
        gp_y, 10**(mu_reshaped[:, n]+sig_reshaped[:, n]),
        10**(mu_reshaped[:, n]+sig_reshaped[:, n]*2),
        facecolor='blue', alpha=0.15)
    plt.fill_between(
        gp_y, 10**(mu_reshaped[:, n]-sig_reshaped[:, n]),
        10**(mu_reshaped[:, n]+sig_reshaped[:, n]),
        facecolor='blue', alpha=0.3, label=r'$\pm 1 \sigma$ interval')
    plt.legend()
    plt.yscale("log")
    plt.show()
    # ploting the T, c space
    fig, ax = plt.subplots()
    X2, Y2 = np.meshgrid(
        np.linspace(320, 1100, num=100),
        np.logspace(20, 22+np.log10(7), num=100, endpoint=True))
    inventory_per_monoblock = compute_inv(X2, Y2)
    levels = np.logspace(
        np.log10(np.min(inventory_per_monoblock)),
        np.log10(np.max(inventory_per_monoblock)),
        # np.log10(1e17),
        # np.log10(2e20),
        1000)
    levels2 = np.logspace(
        # np.log10(np.min(grid_z*1e20)),
        np.log10(np.min(inventory_per_monoblock)),
        np.log10(np.max(inventory_per_monoblock)),
        10)
    locator = ticker.LogLocator(base=10)

    CS = ax.contourf(
        X2, Y2, inventory_per_monoblock,
        levels=levels, locator=locator)
    # plt.xlim(left=320, right=900)
    # plt.ylim(bottom=1e20, top=6e22)
    # Scatter simulation points
    points_2D = []
    points_1D = []
    for p in points:
        if p[0] >= 500:
            points_2D.append(p)
        else:
            points_1D.append(p)
    # all points are grey
    # plt.scatter(points[:, 0], points[:, 1], marker="+", color="grey")
    # red and grey (1D, 2D)
    # plt.scatter(
    #     np.asarray(points_2D)[:, 0], np.asarray(points_2D)[:, 1],
    #     marker="+", color="grey")
    # plt.scatter(
    #     np.asarray(points_1D)[:, 0], np.asarray(points_1D)[:, 1],
    #     marker="+", color="tab:red")

    # add contour lines
    CS2 = ax.contour(
        X2, Y2, inventory_per_monoblock, levels=levels2,
        locator=locator, colors="white")
    CLS = plt.clabel(CS2, inline=True, fontsize=10, fmt=scientificNotation)

    # to remove the ones outside axis go to:
    # https://stackoverflow.com/questions/25873681/matplotlib-contour-plot-labels-overlap-axes
    for c in CS.collections:  # for avoiding white lines in pdf
        c.set_edgecolor("face")
    plt.colorbar(CS, label=r"Inventory per monoblock (H)", ticks=locator)
    plt.yscale("log")
    plt.xlabel(r"$T_{surface}$ (K)")
    plt.ylabel(r"$c_{surface}$ (m$^{-3}$)")
    plt.minorticks_on()
    plt.show()

    # plot abaqus
    T = np.linspace(330, 950, num=100)
    concentrations = np.logspace(20, 22+np.log10(3), num=10)
    for c_surf in concentrations:
        plt.plot(
            T, interp_extrap_2D(T, c_surf),
            label="{:.1e}".format(c_surf),
            color="tab:grey")
        plt.annotate(
            scientificNotation(c_surf) + " m$^{-3}$",
            xy=(T[-1]+5, 0.90*interp_extrap_2D(T[-1], c_surf)),
            color="tab:grey")
    # try:
    #     from labellines import labelLine, labelLines
    #     labelLines(plt.gca().get_lines(), zorder=2.5)
    # except:
    #     pass
    plt.xlim(left=330, right=1110)
    plt.minorticks_on()
    plt.grid(which='minor', alpha=0.3)
    plt.grid(which='major', alpha=0.7)
    plt.yscale("log")
    plt.xlabel(r"$T$ (K)")
    plt.ylabel("Inventory per monoblock (H)")
    plt.show()
