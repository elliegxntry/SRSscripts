import matplotlib.pyplot as plt
import numpy as np
import os

dist = "Beta"

datafile = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/Constant" + dist + "/"

radius = 1.52
mdot_path = datafile + "mdot-data_r{}.txt".format(radius)

with open(mdot_path,"r") as file:
    mdot_data = np.loadtxt(file, skiprows=1)
    times = mdot_data[100:, 0]
    mdot = mdot_data[100:, 1]

plt.plot(times, mdot, ls="", marker="*")
# plt.gca().axhline(0, ls="--", color="black")
plt.xlabel("Time [GM/c^3]")
plt.ylabel("Mass Flux [c^3/G]")
plt.title("Constant " + dist + " Mass Flux through r{}".format(radius))
plt.tight_layout()
filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/mdot_profiles/Constant" + dist + "/specific/"
if not os.path.isdir(filedir):
    os.mkdir(filedir)

filename = "mdot_r{}.png".format(radius)
print(filename)
print(filedir)
plt.savefig(filedir + filename)
plt.show()
plt.close()