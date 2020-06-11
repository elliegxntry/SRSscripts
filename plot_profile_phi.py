# import Python modules
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('GRvis-master/scripts/modules/')
from raw_data_utils import read_athdf

# path to load data
datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
configA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2"
configB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2"
datapath_baseA = datapath + configA
datapath_baseB = datapath + configB

# Specifications
times_to_look_at = np.arange(0, 10)
quantity_to_load = "rho"

# dictionary for quantities
quantities = ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
quantity_names = {"rho":"Density", "press":"Pressure", "vel1":"Radial velocity", "vel2":"Theta velocity",
                  "vel3":"Azimuthal velocity", "Bcc1":"Radial magnetic field", "Bcc2":"Theta magnetic field",
                  "Bcc3":"Azimuthal magnetic field"}

initial_datapath_baseA = "C:/Users/Ellie/Downloads/nerd/SRSData/1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2/"
initial_datapath_baseB = "C:/Users/Ellie/Downloads/nerd/SRSData/1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2/"
initial_filenameA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2.prim.00000.athdf"
initial_filenameB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2.prim.00000.athdf"

# load initial data
initial_dataA = read_athdf(initial_datapath_baseA + initial_filenameA, quantities=[quantity_to_load])
quantity_maxA = np.max(initial_dataA[quantity_to_load])
initial_dataB = read_athdf(initial_datapath_baseA + initial_filenameA, quantities=[quantity_to_load])
quantity_maxB = np.max(initial_dataB[quantity_to_load])


for timestep in times_to_look_at:
    timestep = "{:05d}".format(int(timestep))
    filepathA = datapath_baseA + "/" + configA + ".prim." + timestep + ".athdf"
    filepathB = datapath_baseB + "/" + configB + ".prim." + timestep + ".athdf"
    #print("Loading time step {}".format(timestep))
    dataA = read_athdf(filepathA, quantities=[quantity_to_load])
    dataB = read_athdf(filepathB, quantities=[quantity_to_load])

    #define variables
    simulation_timeA = dataA["Time"]
    simulation_timeB = dataB["Time"]
    quantity_dataA = dataA[quantity_to_load]
    quantity_dataB = dataB[quantity_to_load]
    phi_valuesA = dataA['x3v']
    phi_valuesB = dataB['x3v']
    #print(phi_valuesA[int(phi_valuesA.size)])

    radius_in_codeunits = 2
    radiusind = (np.abs(dataA['x1v'] - radius_in_codeunits)).argmin()

    # get line out at constant radius
    phi_index = 0; radius_index = radiusind; theta_indexA = int(quantity_dataA.shape[1]/2)
    phi_index = 0; radius_index = radiusind; theta_indexB = int(quantity_dataB.shape[1]/2)
    phi_dataA = quantity_dataA[:, theta_indexA, radius_index]
    phi_dataB = quantity_dataB[:, theta_indexB, radius_index]

    # plot according to phi
    plt.plot(phi_valuesA, phi_dataA/quantity_maxA, label = "Constant Beta", linestyle = "--")
    #plt.plot(phi_valuesB, phi_dataB/quantity_maxB, label = "Constant B")
    plt.xlabel("Phi")
    plt.ylabel(quantity_names[quantity_to_load])
    plt.title("Time: {}".format(simulation_timeA))
    plt.legend()
    plt.tight_layout()

    # Save figure
    figname = quantity_to_load + "_at_timestep_{}".format(timestep)
    filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/" + quantity_to_load + "_r2phiprofile/normalized/"
    if not os.path.isdir(filedir):
        os.mkdir(filedir)
    #print(filedir)
    #print("Saving figure " + filedir + figname)
    plt.savefig(filedir + figname)
    plt.show()
    plt.gca().clear()