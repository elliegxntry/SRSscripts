import numpy as np
import matplotlib.pyplot as plt
import os

# Things to modify script - see if statements for times in chunks
dist = "B"
chunk1 = False
chunk2 = True
chunk3 = False
chunk4 = False
chunk5 = False
chunk6 = False
chunk7 = False
all_time = False

all_chunks = True
only_one = False

#pull the data and set up where to save it
data_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/mdotAsR/Constant" + dist + "/"
data_save_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/mdotAsR/time/Constant" + dist + "/"
if not os.path.isdir(data_save_path):
    os.makedirs(data_save_path)
filestem = "mass_flux_over_r_at_t"

#Changes based on desired time
if chunk1:
    timesteps = np.arange(0, 101)
    filename = "time_mdot_over_r_chunk1.txt"
    figname = "tavg_mdot_over_r_chunk1"
    timestep_name = "t=0 GM/c^3 to t=500 GM/c^3"
if chunk2:
    timesteps = np.arange(100, 201)
    filename = "time_mdot_over_r_chunk2.txt"
    figname = "tavg_mdot_over_r_chunk2"
    timestep_name = "t=500 GM/c^3 to t=1000 GM/c^3"
if chunk3:
    timesteps = np.arange(200, 301)
    filename = "time_mdot_over_r_chunk3.txt"
    figname = "tavg_mdot_over_r_chunk3"
    timestep_name = "t=1000 GM/c^3 to t=1500 GM/c^3"
if chunk4:
    timesteps = np.arange(300, 401)
    filename = "time_mdot_over_r_chunk4.txt"
    figname = "tavg_mdot_over_r_chunk4"
    timestep_name = "t=1500 GM/c^3 to t=2000 GM/c^3"
if chunk5:
    timesteps = np.arange(400, 501)
    filename = "time_mdot_over_r_chunk5.txt"
    figname = "tavg_mdot_over_r_chunk5"
    timestep_name = "t=2000 GM/c^3 to t=2500 GM/c^3"
if chunk6:
    timesteps = np.arange(500, 601)
    filename = "time_mdot_over_r_chunk6.txt"
    figname = "tavg_mdot_over_r_chunk6"
    timestep_name = "t=2500 GM/c^3 to t=3000 GM/c^3"
if chunk7:
    timesteps = np.arange(600, 671)
    filename = "time_mdot_over_r_chunk7.txt"
    figname = "tavg_mdot_over_r_chunk7"
    timestep_name = "t=3000 GM/c^3 to t=3350 GM/c^3"
if all_time:
    timesteps = np.arange(0,671)
    filename = "time_mdot_over_r_alltime.txt"
    figname = "tavg_mdot_over_r_alltime"
    timestep_name = "t=0 GM/c^3 to t=3350 GM/c^3"

# create empty array to put time data in
mdot_data = np.array([])

# add section of mdot at each radius for each timestep
for time in timesteps:
    time = "{:05d}".format(int(time))
    mdot_path = data_path + filestem + time + ".txt"
    mdot_data_at_t = np.loadtxt(mdot_path, skiprows=1)
    mdot_data_at_t = mdot_data_at_t.reshape((1,mdot_data_at_t.size))
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
plt.title("Time averaged mass flux" + "\n" + "Constant " + dist + "\n" + timestep_name)
plt.xlim([0, 20])
plt.ylim([-0.05, 0])
plt.tight_layout()
fig_save_path = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/tavg_mdot/Constant" + dist + "/"
if not os.path.isdir(fig_save_path):
    os.makedirs(fig_save_path)
plt.savefig(fig_save_path + figname)
plt.show()
plt.close()