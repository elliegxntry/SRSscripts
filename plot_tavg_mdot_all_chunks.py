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
dist = "B"

# load paths
data_load_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/mdotAsR/time/Constant" + dist + "/"
fig_save_path = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/tavg_mdot/all_chunks/"
if not os.path.isdir(fig_save_path):
    os.makedirs(fig_save_path)
filestem = "time_mdot_over_r_"
datapath = data_load_path + filestem

# get time average for each chunk from file
chunk1 = np.mean(datapath + "chunk1.txt", axis=0)
chunk2 = np.mean(datapath + "chunk2.txt", axis=0)
chunk3 = np.mean(datapath + "chunk3.txt", axis=0)
chunk4 = np.mean(datapath + "chunk4.txt", axis=0)
chunk5 = np.mean(datapath + "chunk5.txt", axis=0)
chunk6 = np.mean(datapath + "chunk6.txt", axis=0)
chunk7 = np.mean(datapath + "chunk7.txt", axis=0)

# plot chunks
plt.figure()
#plt.colormaps('viridis')
plt.plot(chunk1, label="t=0 GM/c^3 to t=500 GM/c^3")
plt.plot(chunk2, label="t=500 GM/c^3 to t=1000 GM/c^3")
plt.plot(chunk3, label="t=1000 GM/c^3 to t=1500 GM/c^3")
plt.plot(chunk4, label="t=1500 GM/c^3 to t=2000 GM/c^3")
plt.plot(chunk5, label="t=2000 GM/c^3 to t=2500 GM/c^3")
plt.plot(chunk6, label="t=2500 GM/c^3 to t=3000 GM/c^3")
plt.plot(chunk7, label="t=3000 GM/c^3 to t=3350 GM/c^3")
plt.xlabel("Radius [GM/c^2]")
plt.ylabel("Mass Flux [c^3/G]")
plt.xlim([0,20])
plt.title("Time averaged mass flux over different time intervals" + "\n" + "Constant " + dist)
plt.tight_layout()
plt.savefig(fig_save_path + "Constant" + dist)
plt.show()
plt.close()
