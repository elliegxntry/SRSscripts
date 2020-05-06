import matplotlib.pyplot as plt
import numpy as np
import os

datafileA = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/ConstantBeta/"
datafileB = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/ConstantB/"


radius = 1.52
mdot_pathA = datafileA + "mdot-data_r{}.txt".format(radius)
mdot_pathB = datafileB + "mdot-data_r{}.txt".format(radius)


with open(mdot_pathA,"r") as file:
    mdot_dataA = np.loadtxt(file, skiprows=1)
    timesA = mdot_dataA[100:, 0]
    mdotA = mdot_dataA[100:, 1]

with open(mdot_pathB,"r") as file:
    mdot_dataB = np.loadtxt(file, skiprows=1)
    timesB = mdot_dataB[100:, 0]
    mdotB = mdot_dataB[100:, 1]


plt.plot(timesA, mdotA, ls="--", color='blue', label="Constant Beta", marker="*")
plt.plot(timesB, mdotB, color='red', label="Constant B", marker="*")
# plt.gca().axhline(0, ls="--", color="black")
plt.xlabel("Time [GM/c^3]")
plt.ylabel("Mass Flux [c^3/G]")
plt.title("Mass Flux through r{}".format(radius))
plt.legend()
plt.tight_layout()
filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/mdot_profiles/Comparison"
if not os.path.isdir(filedir):
    os.mkdir(filedir)

filename = "mdot_r{}.png".format(radius)
print(filename)
print(filedir)
plt.savefig(filedir + filename)
plt.show()
plt.close()