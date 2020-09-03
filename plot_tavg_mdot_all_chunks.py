"""s
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
data_load_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/1.1.1-torus2_b-gz2_a0beta500tor" + dist + \
            "_br32x32x64rl2x2/tavg_mdot_through_shells_over_r/"
fig_save_path = "C:/Users/Ellie/Downloads/nerd/SRSPlots/1.1.1-torus2_b-gz2_a0beta500tor" + dist\
                + "_br32x32x64rl2x2/tavg_mdot/all_times"
if not os.path.isdir(fig_save_path):
    os.makedirs(fig_save_path)
filestem = "time_mdot_over_r_"
datapath = data_load_path + filestem

# load each chunk from file
t0t500 = np.loadtxt(datapath + "t0t500.txt")
t500t1000 = np.loadtxt(datapath + "t500t1000.txt")
t1000t1500 = np.loadtxt(datapath + "t1000t1500.txt")
t1500t2000 = np.loadtxt(datapath + "t1500t2000.txt")
t2000t2500 = np.loadtxt(datapath + "t2000t2500.txt")
t2500t3000 = np.loadtxt(datapath + "t2500t3000.txt")
t3000t3500 = np.loadtxt(datapath + "t3000t3500.txt")
t3500t4000 = np.loadtxt(datapath + "t3500t4000.txt")
t4000t4500 = np.loadtxt(datapath + "t4000t4500.txt")
t4500t5000 = np.loadtxt(datapath + "t4500t5000.txt")
t5000t5500 = np.loadtxt(datapath + "t5000t5500.txt")
t5500t5945 = np.loadtxt(datapath + "t5500t5945.txt")
print(t0t500.shape)

# take mean of each chunk
t0t500_mean = np.mean(t0t500, axis=0)
t500t1000_mean = np.mean(t500t1000, axis=0)
t1000t1500_mean = np.mean(t1000t1500, axis=0)
t1500t2000_mean = np.mean(t1500t2000, axis=0)
t2000t2500_mean = np.mean(t2000t2500, axis=0)
t2500t3000_mean = np.mean(t2500t3000, axis=0)
t3000t3500_mean = np.mean(t3000t3500, axis=0)
t3500t4000_mean = np.mean(t3500t4000, axis=0)
t4000t4500_mean = np.mean(t4000t4500, axis=0)
t4500t5000_mean = np.mean(t4500t5000, axis=0)
t5000t5500_mean = np.mean(t5000t5500, axis=0)
t5500t5945_mean = np.mean(t5500t5945, axis=0)

# plot chunks
plt.figure()
colormap = plt.cm.plasma
plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.plasma(np.linspace(0, 8))))
plt.plot(t0t500_mean, label="t=0 GM/c^3 to t=500 GM/c^3")
plt.plot(t500t1000_mean, label="t=500 GM/c^3 to t=1000 GM/c^3")
plt.plot(t1000t1500_mean, label="t=1000 GM/c^3 to t=1500 GM/c^3")
plt.plot(t1500t2000_mean, label="t=1500 GM/c^3 to t=2000 GM/c^3")
plt.plot(t2000t2500_mean, label="t=2000 GM/c^3 to t=2500 GM/c^3")
plt.plot(t2500t3000_mean, label="t=2500 GM/c^3 to t=3000 GM/c^3")
plt.plot(t3000t3500_mean, label="t=3000 GM/c^3 to t=3350 GM/c^3")
plt.plot(t3500t4000_mean, label="t=3500 GM/c^3 to t=4000 GM/c^3")
plt.plot(t4000t4500_mean, label="t=4000 GM/c^3 to t=4500 GM/c^3")
plt.plot(t4500t5000_mean, label="t=4500 GM/c^3 to t=5000 GM/c^3")
plt.plot(t5000t5500_mean, label="t=5000 GM/c^3 to t=5500 GM/c^3")
plt.plot(t5500t5945_mean, label="t=5500 GM/c^3 to t=5945 GM/c^3")
plt.xlabel("Radius [GM/c^2]")
plt.ylabel("Mass Flux [c^3/G]")
plt.legend()
plt.title("Time averaged mass flux over different time intervals" + "\n" + "Constant " + dist)
#plt.tight_layout()

# save figure
plt.savefig(fig_save_path + "Constant" + dist)
plt.show()
plt.close()
