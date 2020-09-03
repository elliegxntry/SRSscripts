"""
This script will calculate the total vertical flux at a given radius for the entire azimuthal
domain, and plot each radius at each timestep.
BAD SCRIPT
INPUTS:
    -
OUTPUTS:
    -
TO DO:
    -
"""

#import packages
import sys
sys.path.append('C:/Users/Ellie/Downloads/nerd/scripts/GRvis/scripts/modules/')
import numpy as np
import matplotlib.pyplot as plt
import os
import raw_data_utils
import setup_manager
from raw_data_utils import read_athdf
sys.path.append('C:/Users/Ellie/Downloads/nerd/scripts/GRvis/scripts/modules/m_athay/')
from kerrmetric import kerr,fourvector

# specifications
dist = "Beta"
timesteps = np.arange(0, 1)

#load data
config = "1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2"
datapath_base = "C:/Users/Ellie/Downloads/nerd/SRSData/" + config
mdot_reduced_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/" + config + '/vert_mdot_over_r/' + dist + "/"
if not os.path.isdir(mdot_reduced_path):
    os.makedirs(mdot_reduced_path)
fig_save_path = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/vertical_mdot_over_r/"
if not os.path.isdir(fig_save_path):
    os.makedirs(fig_save_path)

for time in timesteps:
    # load data
    timestep = "{:05d}".format(int(time))
    fname = datapath_base + "/" + config + ".prim." + timestep + ".athdf"
    quantities = ['rho', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
    print("Loading time {}".format(time))

    # define variables from data
    data = read_athdf(fname, quantities=quantities)
    r = data['x1v']
    theta = data['x2v']
    phi = data['x3v']
    R = data['x1f']

    # define lengths for reshaping
    nx1 = len(r)
    nx2 = len(theta)
    nx3 = len(phi)

    # reshape variables
    dp = (2 * np.pi / nx3) * np.ones(nx3)
    dt = (np.pi / nx2) * np.ones(nx2)
    Dp = dp.reshape(nx3, 1, 1)
    Dt = dt.reshape(1, nx2, 1)
    dr = np.ones(nx1)

    # tbh i don't know what this is for
    for i in range(nx1):
        dr[i] = R[i + 1] - R[i]

    # reshape and load more variables, including GR variables
    Dr = dr.reshape(1, 1, nx1)
    rho = data['rho']
    ks = kerr(r, theta)
    f = fourvector(r, theta, data['vel1'], data['vel2'], data['vel3'], data['Bcc1'], data['Bcc2'], data['Bcc3'])

    # perform mflux calculation
    k = Dt * Dp * ks.g
    Mflux = f.u1 * k * rho
    total_mflux = np.sum(np.sum(Mflux, axis=0), axis=0)
    output_data = np.array([r, total_mflux])

    # save data to file
    filename = "mdot_over_r_t{}".format(timestep)
    if not os.path.isfile(mdot_reduced_path):
        header = "radius, mdot"
    print(mdot_reduced_path)
    with open(mdot_reduced_path + filename, "a") as f:
        np.savetxt(f, output_data, header=header)

    # Clean up mdot by removing duplicates and sorting
    with open(mdot_reduced_path + filename, "r") as f:
        mdotdata = np.loadtxt(f, skiprows=1)
    unique_times, inds = np.unique(mdotdata[:, 0], return_index=True)
    new_mdotdata = np.transpose(np.stack([unique_times, mdotdata[:, 1][inds]]))
    with open(mdot_reduced_path + filename, "w") as f:
        np.savetxt(f, new_mdotdata, header="radius, mdot")

    #plot mass flux
    figname = "mdot_over_r_t{}".format(timestep)
    plt.figure()
    plt.plot(mdot_reduced_path + filename, skiprows=1)
    plt.xlabel("Radius [GM/c^3]")
    plt.ylabel("Mass flux through midplane")
    titstring = "Constant: " + dist + "\nVertical mass flux through the midplane" \
                                      "\nAzimuthally summed" \
                                      "\nt={}".format(timestep)
    plt.title(titstring)
    plt.tight_layout()
    plt.savefig(fig_save_path + figname)
    plt.show()
    plt.close()