import numpy as np
import matplotlib.pyplot as plt
import os

datafile = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/ConstantB/"
radius = 10
quantities = ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
quantity_names = {"rho":"Density", "press":"Pressure", "vel1":"Radial velocity", "vel2":"Theta velocity",
                  "vel3":"Azimuthal velocity", "Bcc1":"Radial magnetic field", "Bcc2":"Theta magnetic field",
                  "Bcc3":"Azimuthal magnetic field"}
quantity = "rho"
path = datafile + quantity + "-data_r{}.txt".format(radius)

with open(path) as file:
    data = np.loadtxt(file, skiprows=1)
    time_data = data[:,0]
    quantity_data = data[:,1]

plt.plot(time_data, quantity_data, ls="", marker="*")
plt.xlabel("Time [GM/c^3]")
plt.ylabel(quantity_names[quantity] )
plt.title("Constant B at r{}".format(radius))
filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/time_profiles/ConstantB/"
if not os.path.isdir(filedir):
    os.mkdir(filedir)

filename = quantity + "_r{}.png".format(radius)
print(filename)
print(filedir)
plt.savefig(filedir + filename)
#plt.show()
