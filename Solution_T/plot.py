from fenics import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

mesh = Mesh()
XDMFFile("../Mesh_ITER/mesh_domains_44722.xdmf").read(mesh)

V = FunctionSpace(mesh, "CG", 1)
u = Function(V)

field = "Solutions/10000000"
folder = field + "/"
filename = "T.xdmf"
output_filename = "T_1e7"

XDMFFile(folder + filename).read_checkpoint(u, "T")
fig, ax = plt.subplots(figsize=(4.8, 4.8))
CS = plot(u, levels=400, cmap="plasma")

# levels2 = np.linspace(u.vector().min(), u.vector().max(), 10, endpoint=True)
# levels2 = np.linspace(u.vector().min(), u.vector().max(), num=10)
levels2 = np.linspace(400, 1400, num=10)
CS2 = plot(u, mode="contour", levels=levels2, colors="white", linewidths=0.5)
CLS = plt.clabel(CS2, inline=True, fontsize=6, inline_spacing=25, fmt='%.0f')

# -------------------------------------------------------
# check that the labels don't overlap axes
# get limits if they're automatic
xmin, xmax, ymin, ymax = plt.axis()
Dx = xmax-xmin
Dy = ymax-ymin

# check which labels are near a border
keep_labels = []
thresh = 0.10
for label in CLS:
    lx, ly = label.get_position()
    if xmin+thresh*Dx < lx < xmax-thresh*Dx and \
       ymin+thresh*Dy < ly < ymax-thresh*Dy:
        # inlier, redraw it later
        pass
    elif lx < xmin+thresh*Dx:
        lx = xmin+2*thresh*Dx
    keep_labels.append((lx, ly))

for cline in CS2.collections:
    cline.remove()
for label in CLS:
    label.remove()
CS2 = plot(u, mode="contour", levels=levels2, colors="white", linewidths=0.5)
CLS = plt.clabel(CS2, inline=True, fontsize=10,
                 inline_spacing=25, manual=keep_labels, fmt='%.0f')
# -------------------------------------------------------
levels = np.linspace(u.vector().min(), u.vector().max(), num=10)
cb = fig.colorbar(CS, ticks=levels)
cb.ax.set_title(r'K')
for c in CS.collections:  # for avoiding white lines in pdf
    c.set_edgecolor("face")
fig.patch.set_visible(False)
ax.axis('off')
plt.tight_layout()
plt.savefig(folder + output_filename + ".pdf", transparent=True)
plt.savefig(folder + output_filename + ".svg", transparent=True)
