from context import FESTIM
import numpy as np
from FESTIM.generic_simulation import run
import properties
# atom_density  =  density(g/m3)*Na(/mol)/M(g/mol)
atom_density_W = 6.3222e28  # atomic density m^-3
atom_density_Cu = 8.4912e28  # atomic density m^-3
atom_density_CuCrZr = 2.6096e28  # atomic density m^-3

# IDs for edges and surfaces (must be the same as in xdmf files)
id_W = 8
id_Cu = 7
id_CuCrZr = 6

id_top_surf = 9
id_coolant_surf = 10
id_left_surf = 11


# for heat_flux in np.arange(1e5, 1e6, step=2e5):
# for heat_flux in np.logspace(4, 7, endpoint=True):


def compute_min_max_T(heat_fluxes):
    max_T = []
    min_T = []
    for heat_flux in heat_fluxes:
        parameters = {
            "mesh_parameters": {
                "mesh_file": "Mesh_ITER/mesh_domains_44722.xdmf",
                "cells_file": "Mesh_ITER/mesh_domains_44722.xdmf",
                "facets_file": "Mesh_ITER/mesh_lines_44722.xdmf",
            },
            "materials": [
                {
                    # Tungsten
                    "D_0": 2.9e-7,
                    "E_D": 0.39,
                    # "S_0": atom_density_W*1.3e-4,  # at/m3.Pa0.5 (from Grislia 2015)
                    # "E_S": 0.34,  # eV
                    #"S_0": 1.0,  # case without solubility
                    #"E_S": 0.0,  # case without solubility
                    "thermal_cond": properties.thermal_cond_W,
                    "heat_capacity": 1,
                    "rho": properties.rhoCp_W,
                    "id": id_W,
                },
                {
                    # Cu
                    "D_0": 6.6e-7,
                    "E_D": 0.387,
                    # "S_0": 3.12e28,  # at/m3.Pa0.5 (from ITER)
                    # "E_S": 0.572,  # eV
                    #"S_0": 1.0,  # case without solubility
                    #"E_S": 0.0,  # case without solubility
                    "thermal_cond": properties.thermal_cond_Cu,
                    "heat_capacity": 1,
                    "rho": properties.rhoCp_Cu,
                    "id": id_Cu,
                },
                {
                    # CuCrZr
                    "D_0": 3.92e-7,
                    "E_D": 0.418,
                    # "S_0": 4.28e23,  # at/m3.Pa0.5 (from ITER)
                    # "E_S": 0.387,  # eV
                    #"S_0": 1.0,  # case without solubility
                    #"E_S": 0.0,  # case without solubility
                    "thermal_cond": properties.thermal_cond_CuCrZr,
                    "heat_capacity": 1,
                    "rho": properties.rhoCp_CuCrZr,
                    "id": id_CuCrZr,
                },
                ],
                "traps": [],
            "boundary_conditions": [],
            "temperature": {
                "type": "solve_stationary",
                "boundary_conditions": [
                    {
                        "type": "flux",
                        "value": heat_flux,
                        "surfaces": id_top_surf
                    },
                    {
                        "type": "convective_flux",
                        "h_coeff": 70000,
                        "T_ext": 273.15+50,
                        "surfaces": id_coolant_surf
                    }
                    ],
                "source_term": [
                ],
                "initial_condition": 273.15+200
                },
            "solving_parameters": {
                "type": "solve_stationary",
                "newton_solver": {
                    "absolute_tolerance": 1e10,
                    "relative_tolerance": 1e-9,
                    "maximum_iterations": 10,
                }
                },
            "exports": {}
        }
        output = run(parameters, log_level=20)
        T = output["T"]
        max_T.append(T.vector().max())
        min_T.append(T(0, 0.0145))
        print(min_T)
        print(max_T)
    return min_T, max_T

import matplotlib.pyplot as plt

min_T = [323.97886484164712, 324.93614596272505, 327.00004238582727, 331.45352717681266, 341.0806229747626, 361.97244918530663, 407.6893796694178, 509.51752665157051, 744.55934268980957, 1315.0831844141353]
max_T = [324.0850178919163, 325.164979333802, 327.4936693137756, 332.5198982332075, 343.3914814535244, 367.0138243442832, 418.84616837181767, 534.9582565010357, 805.9666988099833, 1471.31936244766]

heat_fluxes = np.logspace(4, 7, endpoint=True, num=30)

min_T, max_T = compute_min_max_T(heat_fluxes)

plt.plot(heat_fluxes, max_T, color="tab:blue")
plt.plot(heat_fluxes, min_T, color="tab:blue")
plt.fill_between(heat_fluxes, min_T, max_T, alpha=0.6)

plt.annotate("Middle", (1.2e7, min_T[-1]), color="tab:blue", fontsize=12)
plt.annotate("Top left corner", (1.2e7, max_T[-1]), color="tab:blue", fontsize=12)

plt.xscale("log")
plt.xlim(right=8e7)
plt.xlabel(r"Heat flux (W m$^{-2}$)", fontsize=15)
plt.ylabel(r"Surface temperature (K)", fontsize=15)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.tight_layout()
plt.savefig("out.pdf")
