from process_baking_data import process_data
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import numpy as np

data_1 = process_data(folder="solution_baking/noninstant_recomb_on_coolant")
data_2 = process_data(folder="solution_baking/instant_recomb_on_coolant")
data_3 = process_data(folder="solution_baking/instant_recomb_on_coolant_and_left")
data_4 = process_data(folder="solution_baking/heat_flux")

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
plt.figure()
handles = []
j = 0
for i in [0, 2, len(data_1)-1]:
    add_handle = True
    for d, style in zip([data_1, data_2, data_3, data_4], ["solid", "dashed", "dashdot", "dotted"]):
        series = d[i]
        if d == data_4:
            t = np.asarray(series["t"])/3600/24 + 6
        else:
            t = np.asarray(series["t"])/3600/24
        plt.plot(t, np.asarray(series["inventory"])/series["inventory"][0]*100, color=cycle[j], linestyle=style)
        if add_handle:
            handles.append(mlines.Line2D([], [], color=cycle[j], label=r"$T_\mathrm{baking} = $" + "{:.2e} K".format(series["T_baking"])))
        add_handle = False
    j += 1
plt.legend(handles=handles)

plt.xlabel("Time (days)")
plt.ylabel("Inventory (%)")
plt.show()
