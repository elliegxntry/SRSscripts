"""
This script plots a variable over each direction at each time step\
Can't do normalized for vel2, Bcc1, Bcc2, Bcc3
INPUTS:
    - Variable
    - Timesteps to plot
    - To normalize or not
    -
OUTPUTS:
    - Plot of each variable at each timestep
TO DO:
    - Add max normalization for initial zero components (probably won't do this cuz it's not worth
     it for the amount I would learn from it)
"""

# import Python modules
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('GRvis/scripts/modules/')
from raw_data_utils import read_athdf

# specifications
times_to_look_at = np.arange(420, 1189)
quantity_to_load = "Bcc3"
normalized = False
radius = True
theta = False
phi = False

# paths to load data
datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
configA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2"
configB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2"
datapath_baseA = datapath + configA
datapath_baseB = datapath + configB

# dictionary for quantities
quantities = ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
if normalized:
    quantity_names = {"rho": r"$\rho(t)/\rho_{0}$", "press": "$P(t)/P_{0}$",
                      "vel1": r"$v^1(t)/v^1_{0}$", "vel2": r"$v^2(t)/v^2_{0}$",
                      "vel3": r"$v^3(t)/v^3_{0}$", "Bcc1": r"$B^1(t)/B^1_{0}$",
                      "Bcc2": r"$B^2(t)/B^2_{0}$", "Bcc3": r"$B^3(t)/B^3_{0}$"}
else:
    quantity_names = {"rho": "$rho$", "press": "$P$", "vel1": "$v^1$", "vel2": "$v^2$", "vel3": "$v^3$",
                      "Bcc1": "$B^1$", "Bcc2": "$B^2$", "Bcc3": "$B^3$"}

# load initial data for normalization
filepatht0A = datapath_baseA + "/" + configA + ".prim.00000.athdf"
filepatht0B = datapath_baseB + "/" + configB + ".prim.00000.athdf"
datat0A = read_athdf(filepatht0A, quantities=[quantity_to_load])
datat0B = read_athdf(filepatht0B, quantities=[quantity_to_load])
quantity_datat0A = datat0A[quantity_to_load]
quantity_datat0B = datat0B[quantity_to_load]

for timestep in times_to_look_at:
    # load timestep data
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

    # get line out at constant radius
    phi_index = 0; radius_index = 0; theta_index = int(quantity_dataA.shape[1]/2)
    titstr = "Time ind: " + timestep + "\nTime: {} $[GM/c^3]$".format(simulation_timeA)

    # get specifics for type of plot and plot both simulations
    if radius:
        radial_valuesA = dataA['x1v']
        radial_valuesB = dataB['x1v']
        radial_dataA = quantity_dataA[phi_index, theta_index, :]
        radial_dataB = quantity_dataB[phi_index, theta_index, :]
        titstr += "\n" + quantity_to_load + " as a function of radius"
        if normalized:
            new_radial_dataA = radial_dataA/quantity_datat0A[phi_index, theta_index, :]
            new_radial_dataB = radial_dataB/quantity_datat0B[phi_index, theta_index, :]
            plt.plot(radial_valuesA, new_radial_dataA, label="Constant Beta", linestyle="--", color='red')
            plt.plot(radial_valuesB, new_radial_dataB, label="Constant B", color='purple')
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/both_sims/direction_quantity_profiles/rprofile/"\
                      + quantity_to_load + "/Normalized/"
            titstr += "\nNormalized to initial value"
        else:
            plt.plot(radial_valuesA, radial_dataA, label="Constant Beta", linestyle="--", color='red')
            plt.plot(radial_valuesB, radial_dataB, label="Constant B", color='purple')
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/both_sims/direction_quantity_profiles/rprofile/" + quantity_to_load + "/"
        plt.xlabel("Radius")

    if theta:
        theta_valuesA = dataA['x2v']
        theta_valuesB = dataB['x2v']
        theta_dataA = quantity_dataA[phi_index, :, radius_index]
        theta_dataB = quantity_dataB[phi_index, :, radius_index]
        titstr += "\n" + quantity_to_load + " as a function of theta"
        if normalized:
            new_theta_dataA = theta_dataA/quantity_datat0A[phi_index, :, radius_index]
            new_theta_dataB = theta_dataB/quantity_datat0B[phi_index, :, radius_index]
            plt.plot(theta_valuesA, new_theta_dataA, label='Constant Beta', linestyle='--', color='red')
            plt.plot(theta_valuesB, new_theta_dataB, label="Constant B", color='purple')
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/both_sims/direction_quantity_profiles/thetaprofile/"\
                      + quantity_to_load + "/Normalized/"
            titstr += "\nNormalized to initial value"
        else:
            plt.plot(theta_valuesA, theta_dataA, label='Constant Beta', ls='--', color='red')
            plt.plot(theta_valuesB, theta_dataB, label='Constant B', color='purple')
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/both_sims/direction_quantity_profiles/thetaprofile/" + quantity_to_load + "/"
        plt.xlabel('Theta')
        plt.xlim([np.pi / 2 - 0.1, np.pi / 2 + 0.1])

    if phi:
        phi_valuesA = dataA['x3v']
        phi_valuesB = dataB['x3v']
        phi_dataA = quantity_dataA[:, theta_index, radius_index]
        phi_dataB = quantity_dataB[:, theta_index, radius_index]
        titstr += "\n" + quantity_to_load + " as a function of phi"
        if normalized:
            new_phi_dataA = phi_dataA/quantity_datat0A[:, theta_index, radius_index]
            new_phi_dataB = phi_dataB/quantity_datat0B[:, theta_index, radius_index]
            plt.plot(phi_valuesA, new_phi_dataA, label='Constant Beta', linestyle='--', color='red')
            plt.plot(phi_valuesB, new_phi_dataB, label="Constant B", color='purple')
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/both_sims/direction_quantity_profiles/phiprofile/"\
                      + quantity_to_load + "/Normalized/"
            titstr += "\nNormalized to initial value"
        else:
            plt.plot(phi_valuesA, phi_dataA, label='Constant Beta', ls='--', color='red')
            plt.plot(phi_valuesB, phi_dataB, label='Constant B', color='purple')
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/both_sims/direction_quantity_profiles/phiprofile/" + quantity_to_load + "/"
        plt.xlabel('Phi')

    # final plot specifications
    plt.ylabel(quantity_names[quantity_to_load])
    plt.title(titstr)
    plt.legend()
    plt.tight_layout()

    #save figure
    if not os.path.isdir(filedir):
        os.makedirs(filedir)
    figname = quantity_to_load + "_at_timestep_{}".format(timestep)
    plt.savefig(filedir + figname)
    #plt.show()
    plt.gca().clear()
print(quantity_to_load)