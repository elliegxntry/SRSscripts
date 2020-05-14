#!/usr/bin/env python
# coding: utf-8

import argparse
import warnings
import numpy as np
import os
import sys
sys.path.append("C:/Users/Ellie/Downloads/nerd/scripts/modules")
import new_athena_read as read

datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
dist = "Beta"

config = "1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2"
datapath_base = datapath + config

data_save_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/Constant" + dist + "/"

data_load_path = datapath_base + config + "/"
radius = 10
phi = 1
theta = np.pi/2.0 - .01
quantity = 'rho'
save_path = data_save_path + quantity + "-data_r{}.txt".format(radius)

# time_steps = np.arange(0, 400)
time_steps = np.arange(0, 448)
header = ""

for time in time_steps:
    timestep = "{:05d}".format(int(time))
    filename = datapath_base + "/" + config + ".prim." + timestep + ".athdf"
    print("Loading time {}".format(time))
    p = read.athdf(filename, quantities=[quantity], x1_min=radius, x1_max=radius, x2_min=theta, x2_max=theta, x3_min=phi, x3_max=phi)

    data = p[quantity]
    # print(data)

    code_time = np.array(p['Time']).reshape(1,)
    output_data = np.array([code_time, data]).reshape((1,2))
    if not os.path.isfile(save_path):
        header="time, " + quantity
    with open(save_path, "a") as f:
        np.savetxt(f, output_data, header=header)

# clean up mdot by removing duplicates and sorting
with open(save_path, "r") as f:
    data = np.loadtxt(f, skiprows=1)

unique_times, inds = np.unique(data[:,0], return_index=True)
new_data = np.transpose(np.stack([unique_times, data[:,1][inds]]))
with open(save_path, "w") as f:
    np.savetxt(f, new_data, header="time, " + quantity)

