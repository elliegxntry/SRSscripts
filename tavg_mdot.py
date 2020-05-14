import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import csv
sys.path.append("C:/Users/Ellie/Downloads/nerd/scripts/modules/")
import new_athena_read

dist = "B"
data_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/mdotAsR/Constant" + dist + "/"
data_save_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/mdotAsR/time/Constant" + dist + "/"
filestem = "mass_flux_over_r_at_t"
timesteps = np.arange(0, 671)

for time in timesteps:
    time = "{:05d}".format(int(time))
    mdot_path = data_path + filestem + time
    with open(mdot_path, "r") as file:
        mdot_data = csv.reader(file, skiprows=1)
        radius = mdot_data[100:, 0]
        mdot = mdot_data[100:, 1]
    with open(data_save_path, "w"):
        np.savetxt(time, radius, mdot)
with open(data_save_path, "r") as file:
    new_data = np.loadtxt(file)
tavg_mdot = np.mean(new_data, axis=2)
plt.figure()
plt.plot(new_data[1], tavg_mdot)
plt.xlabel("Radius")
plt.ylabel("Mass Flux")
plt.title("Time averaged mass flux")
figname = "tavg_mdot_over_r"
fig_save_path = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/tavg_mdot/"
if not os.path.isdirs(fig_save_path):
    os.makedirs(fig_save_path)
plt.savefig(fig_save_path + figname)
plt.show()
plt.close()