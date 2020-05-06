import matplotlib.pyplot as plt
import numpy as np
import os
import new_athena_read as ath

dist = "Beta"
datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2/"
config = "1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2"
savepath = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/"
radius = 5
times = np.arange(0, 671)

quantities = ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
quantity_names = {"rho":"Density", "press":"Pressure", "vel1":"Radial velocity", "vel2":"Theta velocity",
                  "vel3":"Azimuthal velocity", "Bcc1":"Radial magnetic field", "Bcc2":"Theta magnetic field",
                  "Bcc3":"Azimuthal magnetic field"}
quantity = "Bcc1"

for timestep in times:
    timestep = "{:05d}".format(int(timestep))
    filepath = datapath + "/" + config + ".prim." + timestep + ".athdf"
    # Read data
    with open(filepath) as file:
        data = ath.athdf(filepath, quantities=[quantity])

    quantity_data = data[quantity]

    # indexes
    phi_index = 0
    radius_index = 0
    theta_index = int(quantity_data.shape[1] / 2)

    times = quantity_data[0]
    rdata = quantity_data[phi_index, theta_index, :]
    phidata = quantity_data[:, theta_index, radius_index]
    thetadata = quantity_data[phi_index, :, radius_index]

mean_data = np.mean(rdata, phidata, thetadata)

plt.plot(timestep, mean_data, color='purple')
plt.title("Mean " + quantity_names[quantity] + " oscillations" + "\nr={}".format(radius))
plt.xlabel("Time [GM/c^3")
plt.ylabel(quantity + " mean")
figname = quantity + "_mean"
plt.savefig(savepath + figname)
if not os.path.isdir(savepath + figname):
    os.mkdir(savepath + figname)
plt.show()
