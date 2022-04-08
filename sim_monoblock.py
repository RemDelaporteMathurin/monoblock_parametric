from context import FESTIM
from FESTIM.generic_simulation import run
import properties


def simulate_monoblock(c_surf=5e22, T_surf=1200, folder="Solution_temp"):

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

    parameters = {
        "mesh_parameters": {
            # "mesh_file": "Mesh_ITER/mesh_domains_44722.xdmf",
            # "cells_file": "Mesh_ITER/mesh_domains_44722.xdmf",
            # "facets_file": "Mesh_ITER/mesh_lines_44722.xdmf",

            "mesh_file": "Mesh_ITER_58734/mesh_domains_58734.xdmf",
            "cells_file": "Mesh_ITER_58734/mesh_domains_58734.xdmf",
            "facets_file": "Mesh_ITER_58734/mesh_lines_58734.xdmf",

            # "mesh_file": "Mesh_ITER/Mesh_ITER_99950/mesh_domains_99950.xdmf",
            # "cells_file": "Mesh_ITER/Mesh_ITER_99950/mesh_domains_99950.xdmf",
            # "facets_file": "Mesh_ITER/Mesh_ITER_99950/mesh_lines_99950.xdmf",
        },
        "materials": [
            {
                # Tungsten
                "D_0": properties.D_0_W,
                "E_D": properties.E_D_W,
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
                "D_0": properties.D_0_Cu,
                "E_D": properties.E_D_Cu,
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
                "D_0": properties.D_0_CuCrZr,
                "E_D": properties.E_D_CuCrZr,
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
        "traps": [
            {
                "E_k": properties.E_D_W,
                "k_0": properties.D_0_W/(1.1e-10**2)*0.8165/atom_density_W,
                "E_p": 1,
                "p_0": 8.4e12,
                # "density": 5e-4*atom_density_W,
                "density": 1.1e-3*atom_density_W,
                "materials": [id_W]
            },
            # {
            #     "E_k": 0.39,
            #     "k_0": 2.9e-7/(1.1e-10**2)*0.8165/atom_density_W,
            #     "E_p": 1.4,
            #     "p_0": 8.4e12,
            #     "density": 5e-3*atom_density_W*(FESTIM.y > 0.014499),
            #     "materials": [id_W]
            # },
            {
                "E_k": properties.E_D_Cu,
                "k_0": properties.D_0_Cu/(3.61e-10**2)/atom_density_Cu,
                "E_p": 0.5,
                "p_0": 7.98e13,
                "density": 5e-5*atom_density_Cu,
                "materials": [id_Cu]
            },
            {
                "E_k": properties.E_D_CuCrZr,
                "k_0": properties.D_0_CuCrZr/(3.61e-10**2)/atom_density_CuCrZr,
                "E_p": 0.85,
                "p_0": 7.98e13,
                "density": 5e-5*atom_density_CuCrZr,
                "materials": [id_CuCrZr]
            },
            ],

        "boundary_conditions": [
            {
                "type": "dc",
                "surfaces": id_top_surf,
                "value": c_surf,#*(FESTIM.t>=1e5) + (c_surf*(FESTIM.t/1e5))*(FESTIM.t<1e5)
            },
            {
                "type": "recomb",
                "surfaces": id_coolant_surf,
                "Kr_0": 2.9e-14,
                "E_Kr": 1.92,
                "order": 2,
            },
            ],
        "temperature": {
            "type": "solve_stationary",
            "boundary_conditions": [
                # {
                #     "type": "flux",
                #     "value": phi_h,
                #     "surfaces": id_top_surf
                # },
                {
                    "type": "dc",
                    "value": T_surf,
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
            # "type": "solve_stationary",
            "traps_finite_element": 'DG',
            "final_time": 1e7,
            "initial_stepsize": 0.5e5,
            "adaptive_stepsize": {
                "stepsize_change_ratio": 1.1,
                "t_stop": 1e618,
                "stepsize_stop_max": 1/10,
                "dt_min": 1e4,
                },
            "newton_solver": {
                "absolute_tolerance": 1e10,
                "relative_tolerance": 1e-9,
                "maximum_iterations": 10,
            }
            },
        "exports": {
            "xdmf": {
                "functions": ['T', '0', '1', '2', '3', 'retention'],
                "labels": ['T', 'solute', '1', '2', '3', 'retention'],
                "folder": folder,
                "all_timesteps": True,
            },
            "derived_quantities": {
                "total_volume": [
                    {
                        "volumes": [id_W, id_Cu, id_CuCrZr],
                        "field": "solute"
                    },
                    {
                        "volumes": [id_W],
                        "field": "1"
                    },
                    {
                        "volumes": [id_Cu],
                        "field": "2"
                    },
                    {
                        "volumes": [id_CuCrZr],
                        "field": "3"
                    },
                    {
                        "volumes": [id_W, id_Cu, id_CuCrZr],
                        "field": "retention"
                    },
                ],
                "surface_flux": [
                    {
                        "surfaces": [id_coolant_surf, id_left_surf],
                        "field": "solute"
                    }
                ],
                "file": "derived_quantities.csv",
                "folder": folder
            }
        }

    }
    run(parameters, log_level=30)


def simulate_monoblock_1D(c_surf=5e22, T_surf=1200, folder="Solution_temp"):
    print(c_surf, T_surf)
    # atom_density  =  density(g/m3)*Na(/mol)/M(g/mol)
    atom_density_W = 6.3222e28  # atomic density m^-3
    atom_density_Cu = 8.4912e28  # atomic density m^-3
    atom_density_CuCrZr = 2.6096e28  # atomic density m^-3

    # IDs for edges and surfaces (must be the same as in xdmf files)
    id_W = 8
    id_Cu = 7
    id_CuCrZr = 6

    id_top_surf = 1
    id_coolant_surf = 2
    id_left_surf = 11

    size = 8.5e-3
    parameters = {
        "mesh_parameters": {
            "size": size,
            "initial_number_of_cells": 5000,
            "refinements": [
                {
                    "x": 8e-4,
                    "cells": 1000
                },
                {
                    "x": 8e-5,
                    "cells": 1000
                },
                {
                    "x": 1e-6,
                    "cells": 2000
                },
                # {
                #     "x": 1e-8,
                #     "cells": 1000
                # },
                # {
                #     "x": 1e-8,
                #     "cells": 1000
                # },
            ]
        },
        "materials": [
            {
                # Tungsten
                "D_0": 1.9e-7,
                "E_D": 0.2,
                # "S_0": atom_density_W*1.3e-4,  # at/m3.Pa0.5 (from Grislia 2015)
                # "E_S": 0.34,  # eV
                #"S_0": 1.0,  # case without solubility
                #"E_S": 0.0,  # case without solubility
                "thermal_cond": properties.thermal_cond_W,
                "heat_capacity": 1,
                "rho": properties.rhoCp_W,
                "id": id_W,
                "borders": [0, 6e-3],
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
                "borders": [6e-3, 7e-3],
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
                "borders": [7e-3, size],
            },
            ],
        "traps": [
            {
                "E_k": 0.2,
                "k_0": 1.9e-7/(1.1e-10**2)*0.8165/atom_density_W,
                "E_p": 1,
                "p_0": 8.4e12,
                # "density": 5e-4*atom_density_W,
                "density": 1.1e-3*atom_density_W,
                "materials": [id_W]
            },
            # {
            #     "E_k": 0.39,
            #     "k_0": 2.9e-7/(1.1e-10**2)*0.8165/atom_density_W,
            #     "E_p": 1.4,
            #     "p_0": 8.4e12,
            #     "density": 5e-3*atom_density_W*(FESTIM.y > 0.014499),
            #     "materials": [id_W]
            # },
            # {
            #     "E_k": 0.387,
            #     "k_0": 6.6e-7/(3.61e-10**2)/atom_density_Cu,
            #     "E_p": 0.5,
            #     "p_0": 7.98e13,
            #     "density": 5e-5*atom_density_Cu,
            #     "materials": [id_Cu]
            # },
            # {
            #     "E_k": 0.418,
            #     "k_0": 3.92e-7/(3.61e-10**2)/atom_density_CuCrZr,
            #     "E_p": 0.85,
            #     "p_0": 7.98e13,
            #     "density": 5e-5*atom_density_CuCrZr,
            #     "materials": [id_CuCrZr]
            # },
            ],

        "boundary_conditions": [
            {
                "type": "dc",
                "surfaces": id_top_surf,
                # "value": c_surf*(FESTIM.t>=1e5) + (c_surf*(FESTIM.t/1e5))*(FESTIM.t<1e5)
                "value": c_surf
            },
            {
                "type": "recomb",
                "surfaces": id_coolant_surf,
                "Kr_0": 2.9e-14,
                "E_Kr": 1.92,
                "order": 2,
            },
            ],
        "temperature": {
            "type": "solve_stationary",
            "boundary_conditions": [
                # {
                #     "type": "flux",
                #     "value": phi_h,
                #     "surfaces": id_top_surf
                # },
                {
                    "type": "dc",
                    "value": T_surf,
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
            "final_time": 1e6,
            "initial_stepsize": 50,
            "adaptive_stepsize": {
                "stepsize_change_ratio": 1.05,
                "t_stop": 1e618,
                "stepsize_stop_max": 1/10,
                "dt_min": 1e-4,
                },
            "newton_solver": {
                "absolute_tolerance": 1e10,
                "relative_tolerance": 1e-9,
                "maximum_iterations": 20,
            }
            },
        "exports": {
            # "xdmf": {
            #     "checkpoint": False,
            #     "functions": ['0', '1', 'retention'],
            #     # "functions": ['0'],
            #     "labels": ['solute', '1', 'retention'],
            #     # "labels": ['solute'],
            #     "folder": folder,
            #     "all_timesteps": True,
            # },
            "derived_quantities": {
                "total_volume": [
                    {
                        "volumes": [id_W, id_Cu, id_CuCrZr],
                        "field": "solute"
                    },
                    {
                        "volumes": [id_W],
                        "field": "1"
                    },
                    # {
                    #     "volumes": [id_Cu],
                    #     "field": "2"
                    # },
                    # {
                    #     "volumes": [id_CuCrZr],
                    #     "field": "3"
                    # },
                    {
                        "volumes": [id_W, id_Cu, id_CuCrZr],
                        "field": "retention"
                    },
                ],
                "surface_flux": [
                    {
                        "surfaces": [id_coolant_surf, id_left_surf],
                        "field": "solute"
                    }
                ],
                "file": "derived_quantities.csv",
                "folder": folder
            }
        }

    }
    run(parameters, log_level=30)


if __name__ == "__main__":
    c_surf = 1e20
    T_surf = 700
    # simulate_monoblock_1D(c_surf=c_surf, T_surf=T_surf, folder="Solution_instantaneous_recomb_low_temp/T={:.2e};c={:.2e}".format(T_surf, c_surf))
    simulate_monoblock(c_surf=c_surf, T_surf=T_surf, folder="Fields_3/T={:.3e};c={:.2e}".format(T_surf, c_surf))
