from fenics import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


meshfile = "Mesh_ITER/mesh_domains_44722.xdmf"
meshfile = "Mesh_ITER/Mesh_ITER_58734/mesh_domains_58734.xdmf"
field = "T"
field = "retention"
folder = "Solution_T/10000000.0/"
folder = "Fields/T=7.000e+02;c=1.00e+20/"
folder = "Fields/T=1.00e+03;c=1.00e+21/"
fieldfile = folder + "T.xdmf"
fieldfile = folder + "retention.xdmf"

outfile = "retentiongif/"


def scientificNotation(value, pos=0):
    if value == 0:
        return '0'
    else:
        e = np.log10(np.abs(value))
        m = np.sign(value) * 10 ** (e - int(e))
        return r'${:.1f} \times 10^{{{:d}}}$'.format(m, int(e))

fmt_labels = scientificNotation
fmt_colorbar = ticker.FuncFormatter(scientificNotation)

fontsize = 9
inline_spacing = 50
linewidths = 0.5

ext = ["pdf", "svg"]
ext = "pdf"
transparent = True


def load_mesh(meshfile):
    """Loads mesh

    Args:
        meshfile (str): path of the meshfile

    Returns:
        fenics.mesh: the mesh of the solution
    """
    mesh = Mesh()
    XDMFFile(meshfile).read(mesh)
    return mesh


def load_field(mesh, fieldfile, field, counter=-1):
    V = FunctionSpace(mesh, "DG", 1)
    u = Function(V)

    XDMFFile(fieldfile).read_checkpoint(u, field, counter)
    return u


def save(filename, ext, transparent=True):
    """Saves active figure

    Args:
        filename (str): path of the file (without extension)
        ext (str, list): extension(s) of the output file(s)
        (ex: "pdf", "svg", ["png", "pdf"])
        transparent (bool, optional): Sets the background transparent.
        Defaults to True.
    """
    if not isinstance(ext, list):
        ext = [ext]
    for e in ext:
        plt.savefig(filename + "." + e, transparent=transparent)
    return


def create_figure(u, size=(4.8, 4.8), levels_fill="lin",
                  extend="min", vmin=0, isovalues=True,
                  levels_iso="lin",
                  nb_iso=4, inline_labels=False, thresh=0.10,
                  levels_colorbar="lin", nb_colorbar_ticks=8):
    """[summary]

    Args:
        u (fenics.Function): function to plot
        size (tuple, optional): size of the figure. Defaults to (4.8, 4.8).
        levels_fill (str, optional): levels for contourf
        {"lin", "log", array, float}. if "lin" or "log",
        creates a lin or logspace between min and max value. if array, this
        array will be used. Defaults to "lin".
        extend (str, optional): {"neither", "min", "max", "both"}
        Defaults to "neither". Used for extending the plotted data.
        vmin (float, optional): If not None, values below this value will be
        plotted with the same color as vmin. Defaults to None.
        isovalues (bool, optional): If False, isovalues won't be displayed.
        Defaults to True.
        levels_iso (str, optional): levels for contour
        {"lin", "log", array, float}. if "lin" or "log",
        creates a lin or logspace between min and max value. if array, this
        array will be used. Defaults to "lin".
        nb_iso (int, optional): Number of isovalues. Defaults to 10.
        inline_labels (bool, optional): If False, inline labels for isovalues
        won't be displayed. Defaults to True.
        thresh (float, optional): threshold for overlap detection (inline
        labels outside of box). Defaults to 0.10.
        levels_colorbar (str, optional): levels for colorbar
        {"lin", "log", array, float}. if "lin" or "log",
        creates a lin or logspace between min and max value. if array, this
        array will be used. Defaults to "lin".
        nb_colorbar_ticks (int, optional): number of colorbar ticks. Defaults
        to 10.

    Returns:
        plt.fig(): fig object
    """
    fig = plt.figure(figsize=size)

    if levels_fill == "lin":
        levels = np.linspace(
            -1e22, u.vector().max(),
            num=400, endpoint=True)
    elif levels_fill == "log":
        levels = np.logspace(
            np.log10(u.vector().min()), np.log10(u.vector().max()),
            num=400, endpoint=True)
    else:
        levels = levels_fill
    CS = plot(u, levels=levels, extend=extend, vmin=vmin)

    if isovalues:
        if levels_iso == "lin":
            levels_contours = np.linspace(
                2e24, 9.3e24,
                num=nb_iso)
        elif levels_iso == "log":
            levels_contours = np.logspace(
                np.log10(u.vector().min()), np.log10(u.vector().max()),
                num=nb_iso)
        else:
            levels_contours = levels_iso
        CS2 = plot(
            u, mode="contour", levels=levels_contours,
            colors="white", linewidths=linewidths)
        if inline_labels:
            CL = plt.clabel(
                CS2, inline=True, fontsize=fontsize,
                inline_spacing=inline_spacing,
                fmt=fmt_labels)
            CS2, CL = check_inline_labels(
                CS2, CL, thresh=thresh,
                inline_spacing=inline_spacing, fontsize=fontsize,
                fmt=fmt_labels)

    if levels_colorbar == "lin":
        levels = np.linspace(
            0, u.vector().max(),
            num=nb_colorbar_ticks, endpoint=True)

    elif levels_colorbar == "log":
        levels = np.logspace(
            0, u.vector().max(),
            num=nb_colorbar_ticks, endpoint=True)
    else:
        levels = levels_colorbar
    cb = fig.colorbar(CS, ticks=levels, format=fmt_colorbar, extendfrac=0)
    cb.ax.set_title(r'm$^{-3}$')

    for c in CS.collections:  # for avoiding white lines in pdf
        c.set_edgecolor("face")
    fig.patch.set_visible(False)
    plt.axis('off')
    plt.tight_layout()
    return fig


def check_inline_labels(CS, CL, thresh=0.10,
                        inline_spacing=10, fontsize=10, fmt='%.2f'):
    """Checks that the contour labels don't overlapp with axes and remove them
    if so.

    Args:
        CS ([type]): [description]
        CL ([type]): [description]
        thresh (float, optional): [description]. Defaults to 0.10.
        inline_spacing (int, optional): [description]. Defaults to 10.
        fontsize (int, optional): [description]. Defaults to 10.
        fmt (str, optional): [description]. Defaults to '%.2f'.

    Returns:
        [type]: [description]
    """
    # check that the labels don't overlap axes
    # get limits if they're automatic
    xmin, xmax, ymin, ymax = plt.axis()
    Dx = xmax-xmin
    Dy = ymax-ymin

    # check which labels are near a border
    keep_labels = []
    for label in CL:
        lx, ly = label.get_position()
        if xmin+thresh*Dx < lx < xmax-thresh*Dx and \
           ymin+thresh*Dy < ly < ymax-thresh*Dy:
            # inlier, redraw it later
            pass
        elif lx < xmin+thresh*Dx:
            lx = xmin+2*thresh*Dx
        keep_labels.append((lx, ly))

    colors = []
    linewidths = []
    for cline in CS.collections:
        cline.remove()
        # colors.append(cline.get_colors())
        linewidths.append(cline.get_linewidths())
    for label in CL:
        label.remove()
    CS = plot(u, mode="contour", levels=CS.levels, colors="white", linewidths=linewidths)
    CL = plt.clabel(CS, inline=True, fontsize=fontsize,
                    inline_spacing=inline_spacing, manual=keep_labels, fmt=fmt)
    return CS, CL


if __name__ == "__main__":
    mesh = load_mesh(meshfile)
    u = load_field(mesh, fieldfile, field, -1)
    fig = create_figure(u)
    save("out", "pdf", transparent)
