import numpy as np
import matplotlib.pyplot as plt

datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
configA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2"
configB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2"
datapath_baseA = datapath + configA
datapath_baseB = datapath + configB


filenameA = datapath_baseA + "/" + configA + ".hst"
filenameB = datapath_baseB + "/" + configB + ".hst"

data_valuesA = np.loadtxt(filenameA, skiprows=2)
time_valuesA = data_valuesA[:, 0]
mass_valuesA = data_valuesA[:, 2]
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
plt.ylabel("Mass")
figname = "Volumetrically averaged mass"
plt.legend()
plt.savefig(figname)
plt.show()
