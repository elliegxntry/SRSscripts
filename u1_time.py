import numpy as np
import matplotlib.pyplot as plt
import os

datafile = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/ConstantB/"
radius = 10
quantities = ['u1', 'rhou1', 'bsq', 'beta', 'test']
#transformed is into frame of bh
quantity_names = {"u1":"transformed radial velocity", "rhou1":"transformed radius", "bsq":"Magnetic Energy", "beta":"what it sounds like", "test":"what"}
quantity = "u1"
path = datafile + quantity + "-data_r{}.txt".format(radius)

with open(path) as file:
    data = np.loadtxt(file, skiprows=1)
    time_data = data[:,0]
    quantity_data = data[:,1]

plt.plot(time_data, quantity_data, ls="", marker="*")
plt.xlabel("Time [GM/c^3]")
plt.ylabel(quantity)
plt.title("Constant B at r{}".format(radius))
filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/time_profiles/ConstantB/"
if not os.path.isdir(filedir):
    os.mkdir(filedir)

filename = quantity + "_r{}.png".format(radius)
print(filename)
print(filedir)
plt.savefig(filedir + filename)
#plt.show()
