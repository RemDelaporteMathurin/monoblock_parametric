from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import csv

# Read E and phi from Denis et al (NME 2020)
E = []
with open('exposure_conditions_divertor/JET/E.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=';')
    for row in readCSV:
        E.append([float(row[0]), float(row[1])])
phi = []
with open('exposure_conditions_divertor/JET/phi.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=';')
    for row in readCSV:
        phi.append([float(row[0]), float(row[1])])

E = np.asarray(E)
phi = np.asarray(phi)


# create interpolated object E_vs_s and phi_vs_s
E_vs_s = interp1d(E[:, 0] - min(E[:, 0]), E[:, 1], fill_value="extrapolate")
phi_vs_s = interp1d(phi[:, 0] - min(phi[:, 0]), phi[:, 1], fill_value="extrapolate")

# create s array
s = np.arange(1.75 - min(E[:, 0]), 2.2 - min(E[:, 0]), step=0.01)
s = np.append(s, np.arange(2.7 - min(E[:, 0]), 4 - min(E[:, 0]), step=0.01))

# create E_vs_phi : [distance_along_div, E, phi]
E_vs_phi = []
for d in s:
    E_vs_phi.append([d, E_vs_s(d), phi_vs_s(d)])
E_vs_phi = np.asarray(E_vs_phi)

if __name__ == "__main__":
    # plt.scatter(E_vs_phi[:, 0], E_vs_phi[:, 1])
    plt.scatter(E_vs_phi[:, 0], E_vs_phi[:, 2])
    # plt.scatter(phi[:, 0], phi_vs_s(phi[:, 0]))
    # plt.scatter(E[:, 0], E[:, 1])
    # plt.scatter(s, E_vs_s(s))
    # plt.scatter(s, phi_vs_s(s))
    plt.yscale("log")
    # plt.xscale("log")
    plt.show()
