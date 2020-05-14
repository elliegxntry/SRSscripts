import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import csv
sys.path.append(".../GRvis-master/modules")
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

for r in radius:
    with open(data_save_path, "r") as file:
        new_data = np.loadtxt(file)
    r = new_data[1]

