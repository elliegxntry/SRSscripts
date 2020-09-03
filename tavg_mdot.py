"""
This script takes a time average of each chunk of mdot and plots it
INPUTS:
    - Initial distribution
    - which chunk to do
OUTPUTS:
    - A plot of the time averaged chunk as a function of radius
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# Things to modify script - see if statements for times in chunks
dist = "Beta"
t0t500 = True
t500t1000 = False
t1000t1500 = False
t1500t2000 = False
t2000t2500 = False
t2500t3000 = False
t3000t3500 = False
t3500t4000 = False
t4000t4500 = False
t4500t5000 = False
t5000t5500 = False
t5500t5945 = False
all_time = False

ylim_lower = -0.02
ylim_higher = 0

#pull the data and set up where to save it
data_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/1.1.1-torus2_b-gz2_a0beta500tor" + dist + \
            "_br32x32x64rl2x2/mass_flux_through_shells_over_r/"
data_save_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/1.1.1-torus2_b-gz2_a0beta500tor" + dist + \
            "_br32x32x64rl2x2/tavg_mdot_through_shells_over_r/"
if not os.path.isdir(data_save_path):
    os.makedirs(data_save_path)
filestem = "mass_flux_through_shells_over_r_at_t"

#Changes based on desired time
if t0t500:
    timesteps = np.arange(0, 101)
    filename = "time_mdot_over_r_t0t500.txt"
    figname = "tavg_mdot_over_r_t0t500"
    timestep_name = "t=0 GM/c^3 to t=500 GM/c^3"
if t500t1000:
    timesteps = np.arange(100, 201)
    filename = "time_mdot_over_r_t500t1000.txt"
    figname = "tavg_mdot_over_r_t500t1000"
    timestep_name = "t=500 GM/c^3 to t=1000 GM/c^3"
if t1000t1500:
    timesteps = np.arange(200, 301)
    filename = "time_mdot_over_r_t1000t1500.txt"
    figname = "tavg_mdot_over_r_t1000t1500"
    timestep_name = "t=1000 GM/c^3 to t=1500 GM/c^3"
if t1500t2000:
    timesteps = np.arange(300, 401)
    filename = "time_mdot_over_r_t1500t2000.txt"
    figname = "tavg_mdot_over_r_t1500t2000"
    timestep_name = "t=1500 GM/c^3 to t=2000 GM/c^3"
if t2000t2500:
    timesteps = np.arange(400, 501)
    filename = "time_mdot_over_r_t2000t2500.txt"
    figname = "tavg_mdot_over_r_t2000t2500"
    timestep_name = "t=2000 GM/c^3 to t=2500 GM/c^3"
if t2500t3000:
    timesteps = np.arange(500, 601)
    filename = "time_mdot_over_r_t2500t3000.txt"
    figname = "tavg_mdot_over_r_t2500t3000"
    timestep_name = "t=2500 GM/c^3 to t=3000 GM/c^3"
if t3000t3500:
    timesteps = np.arange(600, 701)
    filename = "time_mdot_over_r_t3000t3500.txt"
    figname = "tavg_mdot_over_r_t3000t3500"
    timestep_name = "t=3000 GM/c^3 to t=3500 GM/c^3"
if t3500t4000:
    timesteps = np.arange(700, 801)
    filename = "time_mdot_over_r_t3500t4000.txt"
    figname = "tavg_mdot_over_r_t3500t4000"
    timestep_name = "t=3500 GM/c^3 to t=4000 GM/c^3"
if t4000t4500:
    timesteps = np.arange(800, 901)
    filename = "time_mdot_over_r_t4000t4500.txt"
    figname = "tavg_mdot_over_r_t4000t4500"
    timestep_name = "t=4000 GM/c^3 to t=4500 GM/c^3"
if t4500t5000:
    timesteps = np.arange(900, 1001)
    filename = "time_mdot_over_r_t4500t5000.txt"
    figname = "tavg_mdot_over_r_t4500t5000"
    timestep_name = "t=4500 GM/c^3 to t=5000 GM/c^3"
if t5000t5500:
    timesteps = np.arange(1000, 1101)
    filename = "time_mdot_over_r_t5000t5500.txt"
    figname = "tavg_mdot_over_r_t5000t5500"
    timestep_name = "t=5000 GM/c^3 to t=5500 GM/c^3"
if t5500t5945:
    timesteps = np.arange(1100, 1189)
    filename = "time_mdot_over_r_t5500t5945.txt"
    figname = "tavg_mdot_over_r_t5500t5945"
    timestep_name = "t=5500 GM/c^3 to t=5945 GM/c^3"
if all_time:
    timesteps = np.arange(0, 1189)
    filename = "time_mdot_over_r_alltime.txt"
    figname = "tavg_mdot_over_r_alltime"
    timestep_name = "t=0 GM/c^3 to t=5945 GM/c^3"

# create empty array to put time data in
mdot_data = np.array([])

# add section of mdot at each radius for each timestep
for time in timesteps:
    time = "{:05d}".format(int(time))
    mdot_path = data_path + filestem + time + ".txt"
    mdot_data_at_t = np.loadtxt(mdot_path, skiprows=1)
    mdot_data_at_t = mdot_data_at_t.reshape((1, mdot_data_at_t.size))
    #print(mdot_data_at_t.shape)
    if mdot_data.size == 0:
        mdot_data = mdot_data_at_t
    else:
        mdot_data = np.append(mdot_data, mdot_data_at_t, axis=0)
    #print(mdot_data.shape)

# Save new mdot data and average over time
np.savetxt(data_save_path + filename, mdot_data)
tavg_mdot_over_r = np.mean(mdot_data, axis=0)
#print(tavg_mdot_over_r.shape)

# Plot time average of mass flux over radius
plt.plot(tavg_mdot_over_r)
plt.xlabel("Radius [GM/c^2]")
plt.ylabel("Mass Flux [c^3/G]")
plt.title("Time averaged mass flux through shells" + "\n" + "Constant " + dist + "\n" + timestep_name)
plt.xlim([0, 20])
plt.ylim([ylim_lower, ylim_higher])
plt.tight_layout()
fig_save_path = "C:/Users/Ellie/Downloads/nerd/SRSPlots/1.1.1-torus2_b-gz2_a0beta500tor" + dist\
                + "_br32x32x64rl2x2/tavg_mdot/"
if not os.path.isdir(fig_save_path):
    os.makedirs(fig_save_path)
plt.savefig(fig_save_path + figname)
plt.show()
plt.close()