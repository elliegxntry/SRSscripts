"""
This script takes all the mdots calculated from tavg_mdot.py in each chunk and graphs them on the same plot
INPUTS:
    - Distribution
    - time averaged chunks from another script
OUTPUT:
    - plot of all chunks over radius
"""

# import packages
import numpy as np
import matplotlib.pyplot as plt
import os

# Things to modify script
dist = "Beta"

# load paths
data_load_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/mdotAsR/time/Constant" + dist + "/"
fig_save_path = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/tavg_mdot/all_chunks/"
if not os.path.isdir(fig_save_path):
    os.makedirs(fig_save_path)
filestem = "time_mdot_over_r_"
datapath = data_load_path + filestem

# load each chunk from file
chunk1 = np.loadtxt(datapath + "chunk1.txt")
chunk2 = np.loadtxt(datapath + "chunk2.txt")
chunk3 = np.loadtxt(datapath + "chunk3.txt")
chunk4 = np.loadtxt(datapath + "chunk4.txt")
chunk5 = np.loadtxt(datapath + "chunk5.txt")
chunk6 = np.loadtxt(datapath + "chunk6.txt")
chunk7 = np.loadtxt(datapath + "chunk7.txt")
print(chunk1.shape)

# take mean of each chunk
chunk1_mean = np.mean(chunk1, axis=0)
chunk2_mean = np.mean(chunk2, axis=0)
chunk3_mean = np.mean(chunk3, axis=0)
chunk4_mean = np.mean(chunk4, axis=0)
chunk5_mean = np.mean(chunk5, axis=0)
chunk6_mean = np.mean(chunk6, axis=0)
chunk7_mean = np.mean(chunk7, axis=0)

# plot chunks
plt.figure()
colormap = plt.cm.plasma
plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.plasma(np.linspace(0, 8))))
plt.plot(chunk1_mean, label="t=0 GM/c^3 to t=500 GM/c^3")
plt.plot(chunk2_mean, label="t=500 GM/c^3 to t=1000 GM/c^3")
plt.plot(chunk3_mean, label="t=1000 GM/c^3 to t=1500 GM/c^3")
plt.plot(chunk4_mean, label="t=1500 GM/c^3 to t=2000 GM/c^3")
plt.plot(chunk5_mean, label="t=2000 GM/c^3 to t=2500 GM/c^3")
plt.plot(chunk6_mean, label="t=2500 GM/c^3 to t=3000 GM/c^3")
plt.plot(chunk7_mean, label="t=3000 GM/c^3 to t=3350 GM/c^3")
plt.xlabel("Radius [GM/c^2]")
plt.ylabel("Mass Flux [c^3/G]")
plt.legend()
plt.title("Time averaged mass flux over different time intervals" + "\n" + "Constant " + dist)
plt.tight_layout()

# save figure
plt.savefig(fig_save_path + "Constant" + dist)
plt.show()
plt.close()
