# import packages
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('GRvis-master/scripts/modules/')
from raw_data_utils import read_athdf

#specifications
radius = 10
quantity = "press"
dist = "B"

quantities = ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
quantity_names = {"rho":r"$\rho(t)/\rho_{0,max}$", "press":"$P$", "vel1":"$v^1$", "vel2":"$v^2$",
                  "vel3":"$v^3$", "Bcc1":"$B^1$", "Bcc2":"$B^2$",
                  "Bcc3":"$B^3$"}

initial_data_path = "C:/Users/Ellie/Downloads/nerd/SRSData/1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2/"
initial_filename = "1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2.prim.00000.athdf"

#load initial data
initial_data = read_athdf(initial_data_path + initial_filename, quantities=[quantity])
quantity_max = np.max(initial_data[quantity])
print(quantity_max)

# get initial values
quantity_norm_values = {"rho":"Density", "press":"Pressure", "vel1":"Radial velocity", "vel2":"Theta velocity",
                  "vel3":"Azimuthal velocity", "Bcc1":"Radial magnetic field", "Bcc2":"Theta magnetic field",
                  "Bcc3":"Azimuthal magnetic field"}

# open data
datafile = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/ConstantB/"
path = datafile + quantity + "-data_r{}.txt".format(radius)
with open(path) as file:
    data = np.loadtxt(file, skiprows=1)
    time_data = data[:, 0]
    quantity_data = data[:, 1]

# Plot from time
plt.plot(time_data, quantity_data/quantity_max, ls="", marker="*")
plt.xlabel("Time [GM/c^3]")
plt.ylabel(quantity_names[quantity])
plt.title("Constant " + dist +" at r{}".format(radius) + "\nNormalized to initial maximum value")
filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/time_profiles/ConstantB/normalized/"
if not os.path.isdir(filedir):
    os.makedirs(filedir)

filename = quantity + "_r{}.png".format(radius)
print(filename)
print(filedir)
plt.savefig(filedir + filename)
plt.show()
