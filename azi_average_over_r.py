"""
This script takes an azimuthal average of a quantity and plots it over the radius for each timestep
Azimuthal averages as a function of radius give a more complete picture of the direction of the
magnetic field throughout the torus.
INPUTS:
    - Timesteps to take data from
    - Quantity to average (most useful for B field)
OUTPUT:
    - One plot for each timestep with both initial simulations plotted over radius
"""

# import Python modules
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append("C:/Users/Ellie/Downloads/nerd/scripts/modules/")
import new_athena_read

# path to load data
datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
configA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2"
configB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2"
datapath_baseA = datapath + configA
datapath_baseB = datapath + configB

# specifications
times_to_look_at = np.arange(414, 671)
quantity_to_load = "Bcc3"

# dictionary for quantities
quantities = ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
quantity_names = {"rho":"Density", "press":"Pressure", "vel1":"Radial velocity", "vel2":"Theta velocity",
                  "vel3":"Azimuthal velocity", "Bcc1":"Radial magnetic field", "Bcc2":"Theta magnetic field",
                  "Bcc3":"Azimuthal magnetic field"}

for timestep in times_to_look_at:
    # load data
    timestep = "{:05d}".format(int(timestep))
    filepathA = datapath_baseA + "/" + configA + ".prim." + timestep + ".athdf"
    filepathB = datapath_baseB + "/" + configB + ".prim." + timestep + ".athdf"
    print("Loading timestep {}".format(timestep))
    dataA = new_athena_read.athdf(filepathA, quantities=[quantity_to_load])
    dataB = new_athena_read.athdf(filepathB, quantities=[quantity_to_load])

    # get variables from file
    simulation_timeA = dataA["Time"]
    simulation_timeB = dataB["Time"]
    quantity_dataA = dataA[quantity_to_load]
    quantity_dataB = dataB[quantity_to_load]
    radial_valuesA = dataA['x1v']
    radial_valuesB = dataB['x1v']

    # azimuthally average system
    azi_avgA = np.mean(quantity_dataA, axis=0)
    azi_avgB = np.mean(quantity_dataB, axis=0)

    # get line out at constant radius
    radius_index = 0; theta_indexA = int(quantity_dataA.shape[1]/2)
    radius_index = 0; theta_indexB = int(quantity_dataB.shape[1]/2)
    radial_dataA = azi_avgA[theta_indexA, :]
    radial_dataB = azi_avgB[theta_indexB, :]

    # plot azimuthal average
    plt.plot(radial_valuesA, radial_dataA, label="Constant Beta", linestyle="--")
    plt.plot(radial_valuesB, radial_dataB, label="Constant B")
    plt.xlabel("Radius [GM/c^2]")
    plt.ylabel(quantity_names[quantity_to_load])
    plt.title("Time: {}".format(simulation_timeA) + " GM/c^3")
    plt.legend()

    # save figure
    figname = quantity_to_load + "_at_timestep_{}".format(timestep)
    filedir ="C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/" + quantity_to_load + "_rprofile_aziAverage/"
    if not os.path.isdir(filedir):
        os.mkdir(filedir)
    #print(filedir)
    #print("Saving figure " + filedir + figname)
    plt.savefig(filedir + figname)
    #plt.show()
    plt.gca().clear()

