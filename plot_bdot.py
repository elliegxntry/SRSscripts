import matplotlib.pyplot as plt
import numpy as np
import os

dist = "B"
radius = 1.52
do_Bcc1 = False
do_Bcc2 = False
do_Bcc3 = True

if do_Bcc1:
    datafile = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/Constant" + dist + "/Bcc1/"
if do_Bcc2:
    datafile = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/Constant" + dist + "/Bcc2/"
if do_Bcc3:
    datafile = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/Constant" + dist + "/Bcc3/"
bdot_path = datafile + "bdot-data_r{}.txt".format(radius)

with open(bdot_path,"r") as file:
    bdot_data = np.loadtxt(file, skiprows=1)
    times = bdot_data[100:, 0]
    bdot = bdot_data[100:, 1]

plt.plot(times, bdot, ls="", marker="*")
# plt.gca().axhline(0, ls="--", color="black")
plt.xlabel("Time [GM/c^3]")
if do_Bcc1:
    plt.ylabel("Radial Magnetic Flux [c^3/G]")
    filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/bdot_profiles/Bcc1/Constant" + dist + "/"
    if not os.path.isdir(filedir):
        os.makedirs(filedir)
if do_Bcc2:
    plt.ylabel("Vertical Magnetic Flux [c^3/G]")
    filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/bdot_profiles/Bcc2/Constant" + dist + "/"
    if not os.path.isdir(filedir):
        os.makedirs(filedir)
if do_Bcc3:
    plt.ylabel("Azimuthal Magnetic Flux [c^3/G]")
    filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/bdot_profiles/Bcc3/Constant" + dist + "/"
    if not os.path.isdir(filedir):
        os.makedirs(filedir)
plt.title("Constant " + dist + " Magnetic Flux through r{}".format(radius))
plt.tight_layout()


filename = "bdot_r{}.png".format(radius)
print(filename)
print(filedir)
plt.savefig(filedir + filename)
plt.show()
plt.close()