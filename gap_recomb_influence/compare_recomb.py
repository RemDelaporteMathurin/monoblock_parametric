import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy.interpolate import interp1d


def read_inv(filename, indexes):
    t, inv = [], []
    if not isinstance(indexes, list):
        indexes = [indexes]
    with open(filename) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        next(readCSV)
        for row in readCSV:
            t.append(float(row[0]))
            val = 0
            for index in indexes:
                val += float(row[index])
            inv.append(val * 12e-3)
    return t, inv


t = np.logspace(np.log10(6e4), 7, num=100)
temperatures = [700, 1000]
for T in temperatures:
    t_recomb, inv_recomb = read_inv("T={:.3e};c={:.2e}_no_recomb/derived_quantities.csv".format(T, 1e20), [9, 10, 11])
    t_no_recomb, inv_no_recomb = read_inv("T={:.3e};c={:.2e}_recomb/derived_quantities.csv".format(T, 1e20), [9, 10, 11])
    plt.figure(1)
    line, = plt.plot(t_recomb, inv_recomb, label="{:.0f} K".format(T))
    plt.plot(t_no_recomb, inv_no_recomb, linestyle="--", color=line.get_color())

    plt.figure(2)
    interp_recomb = interp1d(t_recomb, inv_recomb)(t)
    interp_norecomb = interp1d(t_no_recomb, inv_no_recomb)(t)
    plt.plot(t, 100*abs(interp_norecomb-interp_recomb)/interp_recomb, label="{:.0f} K".format(T))

plt.figure(1)
plt.yscale("log")
plt.xscale("log")
plt.minorticks_on()
plt.grid(which='minor', alpha=0.3)
plt.grid(which='major', alpha=0.7)
plt.xlabel("Time (s)")
plt.ylabel("Inventory (H)")
plt.legend()


plt.figure(2)
plt.xscale("log")
plt.minorticks_on()
plt.grid(which='minor', alpha=0.3)
plt.grid(which='major', alpha=0.7)
# plt.yscale("log")
plt.xlabel("Time (s)")
plt.ylabel("Relative difference (%)")
plt.show()