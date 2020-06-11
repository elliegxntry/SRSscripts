"""
This script takes the average of each quantity over the whole volume
No calculations, takes things from the .hst file
Does both configurations
*** change values in variable name, plot, ylabel and figname
INPUTS:
    - mass
    - momentum (in each direction)
    - kinetic energy (in each direction)
    - mechanical energy (in each direction)
    - total energy
OUTPUTS:
    - Plots as a function of time
"""

# import packages
import numpy as np
import matplotlib.pyplot as plt
import os

# specifications
variable_name = "mass"

# path to load data
datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
configA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2"
configB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2"
datapath_baseA = datapath + configA
datapath_baseB = datapath + configB

# path to save data
filenameA = datapath_baseA + "/" + configA + ".hst"
filenameB = datapath_baseB + "/" + configB + ".hst"
fig_save_path = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/time_quantity_profiles/hstplots/"
if not os.path.isdir(fig_save_path):
    os.makedirs(fig_save_path)

data_valuesA = np.loadtxt(filenameA, skiprows=2)
time_valuesA = data_valuesA[:, 0]
time_shape = time_valuesA.shape
mass_valuesA = data_valuesA[:, 2].reshape(time_shape)
rad_mom_values_A = data_valuesA[:, 3]
theta_mom_values_A = data_valuesA[:, 4]
phi_mom_values_A = data_valuesA[:, 5]
rad_KE_values_A = data_valuesA[:, 6]
theta_KE_values_A = data_valuesA[:, 7]
phi_KE_values_A = data_valuesA[:, 8]
total_energy_values_A = data_valuesA[:, 9]
rad_ME_values_A = data_valuesA[:, 10]
theta_ME_values_A = data_valuesA[:, 11]
phi_ME_values_A = data_valuesA[:, 12]

data_valuesB = np.loadtxt(filenameB, skiprows=2)
time_valuesB = data_valuesB[:, 0]
mass_valuesB = data_valuesB[:, 2]
rad_mom_values_B = data_valuesB[:, 3]
theta_mom_values_B = data_valuesB[:, 4]
phi_mom_values_B = data_valuesB[:, 5]
rad_KE_values_B = data_valuesB[:, 6]
theta_KE_values_B = data_valuesB[:, 7]
phi_KE_values_B = data_valuesB[:, 8]
total_energy_values_B = data_valuesB[:, 9]
rad_ME_values_B = data_valuesB[:, 10]
theta_ME_values_B = data_valuesB[:, 11]
phi_ME_values_B = data_valuesB[:, 12]

plt.plot(time_valuesA, mass_valuesA, label="Constant Beta")
plt.plot(time_valuesB, mass_valuesB, label="Constant B")
plt.yscale("log")
plt.xlabel("Time [GM/c^3]")
plt.ylabel(variable_name)
figname = "vavg_" + variable_name
plt.legend()
plt.savefig(fig_save_path + figname)
plt.show()
