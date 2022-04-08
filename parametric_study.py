from sim_monoblock import simulate_monoblock, simulate_monoblock_1D
from scipy.interpolate import interp1d
import numpy as np
from os import path
import matplotlib.pyplot as plt


# for T_surf in [400, 500, 650, 800, 1000, 1200, 1320, 1500, 1600, 2100, 2500, 3000]:
# for T_surf in [500, 650, 800, 1000, 1200, 1320, 1500, 1600, 2100, 2500, 3000]:
    # for c_surf in np.logspace(20, 23, num=7, endpoint=True):
# create (50, 50) points

points = np.c_[np.random.uniform(low=400, high=405, size=(3,)), np.random.uniform(low=20+np.log10(4), high=20+np.log10(6), size=(3,))]
print(points)
for point in points:
    T_surf = point[0]
    c_surf = 10**point[1]
    if T_surf < 500:
        sim = simulate_monoblock_1D
        folder = "Solution_instantaneous_recomb_low_temp"
    else:
        sim = simulate_monoblock
        folder = "Solution_instantaneous_recomb_rand"
    if not path.exists(folder + "/T={:.2e};c={:.2e}".format(T_surf, c_surf)): # quite likely...
        print(T_surf, c_surf)
        sim(
            c_surf=c_surf,
            T_surf=T_surf,
            folder=folder + "/" + "T={:.2e}".format(T_surf) +
            ";c={:.2e}".format(c_surf))
