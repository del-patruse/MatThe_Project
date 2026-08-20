# Teo Tassa 03.08.2026
import numpy as np
import matplotlib.pyplot as plt

from mat_params import MatParams
from Model import Model
from Solver import hill_params_from_yield_stresses, elastic_params_from_E_nu, build_von_mises_params, Driver
from plots_and_verifs import (verify_von_mises_recovery, plot_load_history, plot_stress_response,
                               plot_plastic_strains, plot_stress_strain, plot_yield_locus)



#--------------------------------------
#define what the code should do:
#1: run the same load cases through the von Mises special case and check that Hill recovers von Mises
#2: plot load history, stress response, plastic strains, stress-strain and yield locu
#--------------------------------------
task = 2
#--------------------------------------


#--------------------------------------
#deifnig the load type
#--------------------------------------
load_type = "ramp" # "ramp" or "cyclic"#
#--------------------------------------


#--------------------------------------
#defining task-specific values
#--------------------------------------
E = 200000
nu = .3
sy11 = 580
sy22 = 510
sy33 = 460
sy12 = 590
sy23 = 530
sy13 = 550
#--------------------------------------
#defining task-specific values
#--------------------------------------


#--------------------------------------
#defining time points
#--------------------------------------
deltaT = .5
match load_type:
    case "ramp":
        time_points = np.arange(0, 100 + deltaT, deltaT)
    case "cyclic":
        time_points = np.arange(0, 400 + deltaT, deltaT)
    case _:
        raise ValueError(f"Unknown load type: {load_type}")



#--------------------------------------
#generating Material parammeter objects
#--------------------------------------
if task == 1:
    hillMatParams =  MatParams(elastic_params=elastic_params_from_E_nu(E, nu), plastic_params=hill_params_from_yield_stresses(sy11, sy11, sy11, sy11/np.sqrt(3), sy11/np.sqrt(3), sy11/np.sqrt(3)))
else:
    hillMatParams =  MatParams(elastic_params=elastic_params_from_E_nu(E, nu), plastic_params=hill_params_from_yield_stresses(sy11, sy22, sy33, sy12, sy23, sy13))

misesMatParams = MatParams(elastic_params=elastic_params_from_E_nu(E, nu), plastic_params=build_von_mises_params(sigma_y=sy11))
#--------------------------------------
#generating Material parammeter objects
#--------------------------------------



#--------------------------------------
#defining loads
#--------------------------------------
def ramp_load(t)-> float:
    # 0 -> 0.05 over t in [0, 100]
    m = 0.05 / 100
    return m * t


def cyclic_load(t)-> float:
    #values from plot:
    f    = .01 #frequency [1/s]
    ampl = 0.05 #strain amplitude
    # triangular wave, starting at 0 and rising: closed form via arcsin(sin(.))
    return ampl * (2 / np.pi) * np.arcsin(np.sin(2 * np.pi * f * t))
#--------------------------------------
#defining loads
#--------------------------------------







def main():
    params = hillMatParams
    model = Model(params)
    driver = Driver(model, prescribed_components={0, 1})

    #run the sim
    match load_type:
                    case "ramp":
                        history = driver.run(ramp_load, time_points=time_points)
                    case "cyclic":
                        history = driver.run(cyclic_load, time_points=time_points)

    match task:         
        case 1:
            vm_model = Model(misesMatParams)
            match load_type:
                case "ramp":
                    vm_history = Driver(vm_model, prescribed_components={0, 1}).run(ramp_load, time_points=time_points)
                case "cyclic":
                    vm_history = Driver(vm_model, prescribed_components={0, 1}).run(cyclic_load, time_points=time_points)
            verify_von_mises_recovery(history, vm_history)

        case 2:
            plot_load_history(history)
            plot_stress_response(history)
            plot_plastic_strains(history)
            plot_stress_strain(history)
            plot_yield_locus(history, model.plastic)



if __name__ == "__main__":
    main()
