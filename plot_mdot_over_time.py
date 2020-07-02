"""
Plots mdot at a specific radius and specific distribution over time
INPUTS:
    - distribution
    - radius
    - takes mdot data from the mdot-data_r_.txt file
OUTPUTS:
    - plot of mdot over time
"""

# import packages
import matplotlib.pyplot as plt
import numpy as np
import os

# Specifications
dist = "B"
radius = 1.52

# Path to load and save data
datafile = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2/"
mdot_path = datafile + "mdot-data_r{}.txt".format(radius)
filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/time_quantity_profiles/mdot_profiles/Constant" + dist + "/specific/"

with open(mdot_path, "r") as file:
    mdot_data = np.loadtxt(file, skiprows=1)
    times = mdot_data[100:, 0]
    mdot = mdot_data[100:, 1]

# Plot data
plt.figure()
plt.plot(times, mdot, ls="", marker="*")
#plt.gca().axhline(0, ls="--", color="black")
plt.xlabel("Time [GM/c^3]")
plt.ylabel("Mass Flux [c^3/G]")
plt.title("Constant " + dist + " Mass Flux through r{}".format(radius))
plt.tight_layout()
if not os.path.isdir(filedir):
    os.mkdir(filedir)

filename = "mdot_r{}.png".format(radius)
#print(filename)
#print(filedir)
plt.savefig(filedir + filename)
plt.show()
plt.close()