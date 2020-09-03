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
sys.path.append('GRvis/scripts/modules/')
from raw_data_utils import read_athdf


# path to load data
datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
configA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2"
configB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2"
datapath_baseA = datapath + configA
datapath_baseB = datapath + configB

# specifications
times_to_look_at = np.arange(0, 2)
quantity_to_load = "Bcc1"
normalized = True

# dictionary for quantities
quantities = ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
if normalized:
    quantity_names = {"rho": r"$\rho(t)/\rho_{0,max}$", "press": "$P(t)/P_{0,max}$",
                      "vel1": r"$v^1(t)/v^1_{0,max}$", "vel2": r"$v^2(t)/v^2_{0,max}$",
                      "vel3": r"$v^3(t)/v^3_{0,max}$", "Bcc1": r"$B^1(t)/B^1_{0,max}$",
                      "Bcc2": r"$B^2(t)/B^2_{0,max}$", "Bcc3": r"$B^3(t)/B^3_{0,max}$"}
    # load initial data
    initial_datapath_baseA = "C:/Users/Ellie/Downloads/nerd/SRSData/1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2/"
    initial_datapath_baseB = "C:/Users/Ellie/Downloads/nerd/SRSData/1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2/"
    initial_filenameA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2.prim.00000.athdf"
    initial_filenameB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2.prim.00000.athdf"
    initial_dataA = read_athdf(initial_datapath_baseA + initial_filenameA, quantities=[quantity_to_load])
    quantity_maxA = np.max(initial_dataA[quantity_to_load])
    initial_dataB = read_athdf(initial_datapath_baseA + initial_filenameA, quantities=[quantity_to_load])
    quantity_maxB = np.max(initial_dataB[quantity_to_load])
else:
    quantity_names = {"rho": "$rho$", "press": "$P$", "vel1": "$v^1$", "vel2": "$v^2$", "vel3": "$v^3$",
                      "Bcc1": "$B^1$", "Bcc2": "$B^2$", "Bcc3": "$B^3$"}

for timestep in times_to_look_at:
    # load data
    timestep = "{:05d}".format(int(timestep))
    filepathA = datapath_baseA + "/" + configA + ".prim." + timestep + ".athdf"
    filepathB = datapath_baseB + "/" + configB + ".prim." + timestep + ".athdf"
    print("Loading timestep {}".format(timestep))
    dataA = read_athdf(filepathA, quantities=[quantity_to_load])
    dataB = read_athdf(filepathB, quantities=[quantity_to_load])

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
    plt.figure()
    plt.plot(radial_valuesA, radial_dataA, label="Constant Beta", linestyle="--", color='red')
    plt.plot(radial_valuesB, radial_dataB, label="Constant B", color='purple')
    plt.xlabel("Radius $[GM/c^2]$")
    plt.ylabel(quantity_names[quantity_to_load])
    titstr = "Time: {}".format(simulation_timeA) + " $GM/c^3$"
    if normalized:
        titstr += "\nNormalized to maximuminitial value"
    plt.title(titstr)
    plt.legend()

    # save figure
    figname = quantity_to_load + "_at_timestep_{}".format(timestep)
    filedir ="C:/Users/Ellie/Downloads/nerd/SRSPlots/both_sims/direction_quantity_profiles/" \
             + quantity_to_load + "_rprofile_aziAverage/"
    if normalized:
        filedir += "/Normalized/"
    if not os.path.isdir(filedir):
        os.mkdir(filedir)
    plt.savefig(filedir + figname)
    #plt.show()
    plt.gca().clear()