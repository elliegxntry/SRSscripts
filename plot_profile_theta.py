# import Python modules
import numpy as np
import matplotlib.pyplot as plt
import os

# Athena++ modules
from scripts import athena_read

# paths to load data
datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
configA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2"
configB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2"
datapath_baseA = datapath + configA
datapath_baseB = datapath + configB

times_to_look_at = np.arange(0, 5)
quantity_to_load = "Bcc2"

# dictionary for quantities
quantities = ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
quantity_names = {"rho":"Density", "press":"Pressure", "vel1":"Radial velocity", "vel2":"Theta velocity",
                  "vel3":"Azimuthal velocity", "Bcc1":"Radial magnetic field", "Bcc2":"Theta magnetic field",
                  "Bcc3":"Azimuthal magnetic field"}

for timestep in times_to_look_at:
    timestep = "{:05d}".format(int(timestep))
    filepathA = datapath_baseA + "/" + configA + ".prim." + timestep + ".athdf"
    filepathB = datapath_baseB + "/" + configB + ".prim." + timestep + ".athdf"
    #print("Loading time step {}".format(timestep))
    dataA = athena_read.athdf(filepathA, quantities=[quantity_to_load])
    dataB = athena_read.athdf(filepathB, quantities=[quantity_to_load])

    # define variables
    simulation_timeA = dataA["Time"]
    simulation_timeB = dataB["Time"]
    quantity_dataA = dataA[quantity_to_load]
    quantity_dataB = dataB[quantity_to_load]
    theta_valuesA = dataA['x2v']
    theta_valuesB = dataB['x2v']
    #print(theta_valuesA[int(theta_valuesA.size/2)])

    radius_in_codeunits = 5
    radiusind = (np.abs(dataA['x1v'] - radius_in_codeunits)).argmin()

    # get line out at variable theta
    phi_index = 0; radius_index = radiusind; theta_indexA = int(quantity_dataA.shape[1]/2)
    phi_index = 0; radius_index = radiusind; theta_indexB = int(quantity_dataB.shape[1]/2)
    theta_dataA = quantity_dataA[phi_index, :, radius_index]
    theta_dataB = quantity_dataB[phi_index, :, radius_index]
    plt.plot(theta_valuesA, theta_dataA, label = "Constant Beta", linestyle = "--")
    plt.plot(theta_valuesB, theta_dataB, label = "Constant B")
    plt.xlabel("Theta")
    plt.xlim([np.pi/2-0.1, np.pi/2+0.1])
    plt.ylabel(quantity_names[quantity_to_load])
    plt.title("Time: {}".format(simulation_timeA))
    plt.legend()
    plt.tight_layout()

    # Save figure
    figname = quantity_to_load + "_at_timestep_{}".format(timestep)
    filedir ="C:/Users/Ellie/Downloads/nerd/SRSProfiles/" + quantity_to_load + "_thetaprofile_c01/"
    if not os.path.isdir(filedir):
        os.mkdir(filedir)
    #print(filedir)
    #print("Saving figure " + filedir + figname)
    plt.savefig(filedir + figname)
    plt.show()
    plt.gca().clear()