"""
Calculates a mass flux through a given radius at each timestep
calculation is GR
INPUTS:
    - radius to calculate mdot at
    - initial distribution
    - timesteps to perform calculation on
OUTPUTS:
    - mdot for each timestep in an array
        - each array is separated by radius
"""

#import packages
import numpy as np
import os
import sys
sys.path.append("C:/Users/Ellie/Downloads/nerd/scripts/modules/")
import new_athena_read as read
from kerrmetric import kerr,fourvector

# Specifications
dist = "Beta"
radius = 1.52
time_steps = np.arange(440, 441)

# Paths to load and save
datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
config = "1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2"
datapath_base = datapath + config
data_load_path = datapath_base + config + "/"
data_save_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2/"
mdot_path = data_save_path + "mdot-data_r{}.txt".format(radius)
header = ""

# Calculate mdot at each timestep
for time in time_steps:
    # load data
    timestep = "{:05d}".format(int(time))
    fname = datapath_base + "/" + config + ".prim." + timestep + ".athdf"
    quantities = ['rho', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
    #print("Loading time {}".format(time))

    # define variables from data
    p = read.athdf(fname, quantities=quantities, x1_min=radius, x1_max=radius)
    r = p['x1v']
    theta = p['x2v']
    phi = p['x3v']
    R = p['x1f']

    # define lengths for reshaping
    nx1 = len(r)
    nx2 = len(theta)
    nx3 = len(phi)

    # reshape variables
    dp = (2*np.pi/nx3)*np.ones(nx3)
    dt = (np.pi/nx2)*np.ones(nx2)
    Dp = dp.reshape(nx3,1,1)
    Dt = dt.reshape(1,nx2,1)
    dr = np.ones(nx1)

    #tbh i don't know what this is for
    for i in range(nx1):
        dr[i] = R[i+1]-R[i]

    # reshape and load more variables, including GR variables
    Dr = dr.reshape(1,1,nx1)
    rho = p['rho']
    ks = kerr(r, theta)
    f = fourvector(r, theta, p['vel1'], p['vel2'], p['vel3'], p['Bcc1'], p['Bcc2'], p['Bcc3'])

    # perform mflux calculation
    k = Dt*Dp*ks.g
    Mflux = f.u1*k*rho
    total_mflux = np.sum(np.sum(Mflux, axis=0), axis=0)
    code_time = np.array(p['Time']).reshape(1,)
    output_data = np.array([code_time, total_mflux]).reshape((1,2))

    # save data to file
    if not os.path.isfile(mdot_path):
        header = "time, mdot"
    with open(mdot_path, "a") as f:
        np.savetxt(f, output_data, header=header)

# Clean up mdot by removing duplicates and sorting
with open(mdot_path, "r") as f:
    mdotdata = np.loadtxt(f, skiprows=1)
unique_times, inds = np.unique(mdotdata[:, 0], return_index=True)
new_mdotdata = np.transpose(np.stack([unique_times, mdotdata[:, 1][inds]]))
with open(mdot_path, "w") as f:
    np.savetxt(f, new_mdotdata, header="time, mdot")