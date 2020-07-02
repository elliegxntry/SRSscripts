# import Python modules
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('GRvis/scripts/modules/')
from raw_data_utils import read_athdf

# paths to load data
datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
configA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2"
configB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2"
datapath_baseA = datapath + configA
datapath_baseB = datapath + configB

times_to_look_at = np.arange(500, 501)
quantity_to_load = "vel2"

# dictionary for quantities
quantities = ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
quantity_names = {"rho":r"$\rho(t)/\rho_{0,max}$", "press":"$P$", "vel1":"$v^1$", "vel2":"$v^2$",
                  "vel3":"$v^3$", "Bcc1":"$B^1$", "Bcc2":"$B^2$",
                  "Bcc3":"$B^3$"}

for timestep in times_to_look_at:
    timestep = "{:05d}".format(int(timestep))
    filepathA = datapath_baseA + "/" + configA + ".prim." + timestep + ".athdf"
    filepathB = datapath_baseB + "/" + configB + ".prim." + timestep + ".athdf"
    print("Loading time step {}".format(timestep))
    dataA = read_athdf(filepathA, quantities=[quantity_to_load])
    dataB = read_athdf(filepathB, quantities=[quantity_to_load])

    # Get variables
    simulation_timeA = dataA["Time"]
    simulation_timeB = dataB["Time"]
    quantity_dataA = dataA[quantity_to_load]
    quantity_dataB = dataB[quantity_to_load]
    radial_valuesA = dataA['x1v']
    radial_valuesB = dataB['x1v']
    #print(radial_valuesA[int(radial_valuesA.size/2)])

    radius_in_codeunits = 5
    radiusind = (np.abs(dataA['x1v'] - radius_in_codeunits)).argmin()
   # print(radiusind)

    # get line out at constant radius
    phi_index = 0; radius_index = 0; theta_indexA = int(quantity_dataA.shape[1]/2)
    phi_index = 0; radius_index = 0; theta_indexB = int(quantity_dataB.shape[1]/2)
    radial_dataA = quantity_dataA[phi_index, theta_indexA, :]
    radial_dataB = quantity_dataB[phi_index, theta_indexB, :]
    plt.plot(radial_valuesA, radial_dataA, label = "Constant Beta", linestyle = "--")
    plt.plot(radial_valuesB, radial_dataB, label = "Constant B")
    plt.xlabel("Radius")
    plt.ylabel(quantity_names[quantity_to_load])
    plt.title("Time: {}".format(simulation_timeA))
    plt.legend()

    # save figure
    figname = quantity_to_load + "_at_timestep_{}".format(timestep)
    filedir ="C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/" + quantity_to_load + "_rprofile/"
    if not os.path.isdir(filedir):
        os.mkdir(filedir)
    print(filedir)
    print("Saving figure " + filedir + figname)
    plt.savefig(filedir + figname)
    # plt.show()
    plt.gca().clear()

