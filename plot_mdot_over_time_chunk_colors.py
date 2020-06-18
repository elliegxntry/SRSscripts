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
filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/time_quantity_profiles/mdot_profiles/Constant" + dist + "/chunks/"
if not os.path.isdir(filedir):
    os.makedirs(filedir)

with open(mdot_path, "r") as file:
    mdot_data = np.loadtxt(file, skiprows=1)
    chunk1_times = mdot_data[0:100, 0]
    chunk2_times = mdot_data[100:200, 0]
    chunk3_times = mdot_data[200:300, 0]
    chunk4_times = mdot_data[300:400, 0]
    chunk5_times = mdot_data[400:500, 0]
    chunk6_times = mdot_data[500:600, 0]
    chunk7_times = mdot_data[600:671, 0]
    chunk1_mdot = mdot_data[0:100, 1]
    chunk2_mdot = mdot_data[100:200, 1]
    chunk3_mdot = mdot_data[200:300, 1]
    chunk4_mdot = mdot_data[300:400, 1]
    chunk5_mdot = mdot_data[400:500, 1]
    chunk6_mdot = mdot_data[500:600, 1]
    chunk7_mdot = mdot_data[600:671, 1]

# Plot data
plt.figure()
colormap = plt.cm.plasma
plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.plasma(np.linspace(0, 8))))
plt.plot(chunk1_times, chunk1_mdot, ls="", marker="*", label="Chunk 1")
plt.plot(chunk2_times, chunk2_mdot, ls="", marker="*", label="Chunk 2")
plt.plot(chunk3_times, chunk3_mdot, ls="", marker="*", label="Chunk 3")
plt.plot(chunk4_times, chunk4_mdot, ls="", marker="*", label="Chunk 4")
plt.plot(chunk5_times, chunk5_mdot, ls="", marker="*", label="Chunk 5")
plt.plot(chunk6_times, chunk6_mdot, ls="", marker="*", label="Chunk 6")
plt.plot(chunk7_times, chunk7_mdot, ls="", marker="*", label="Chunk 7")
#plt.xlim([0, 500])
#plt.gca().axhline(0, ls="--", color="black")
plt.xlabel("Time [GM/c^3]")
plt.ylabel("Mass Flux [c^3/G]")
plt.legend()
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