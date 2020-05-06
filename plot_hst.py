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
data_valuesB = np.loadtxt(filenameB, skiprows=2)
time_valuesA = data_valuesA[:, 0]
mass_valuesA = data_valuesA[:, 2]
time_valuesB = data_valuesB[:, 0]
mass_valuesB = data_valuesB[:, 2]

plt.plot(time_valuesA, mass_valuesA, label="Constant Beta")
plt.plot(time_valuesB, mass_valuesB, label="Constant B")
plt.yscale("log")
plt.xlabel("Time [GM/c^3]")
plt.ylabel("Mass")
figname = "Volumetrically averaged mass"
plt.legend()
plt.savefig(figname)
plt.show()
