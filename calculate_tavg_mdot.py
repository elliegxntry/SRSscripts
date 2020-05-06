#!/usr/bin/env python
# coding: utf-8

import argparse
import warnings
import numpy as np
import os
import sys
sys.path.append("C:/Users/Ellie/Downloads/nerd/athena-public-version/vis/python/modules")
import new_athena_read as read
from kerrmetric import kerr,fourvector
import matplotlib.pyplot as plt

datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
dist = "Beta"
config = "1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2"
datapath_base = datapath + config
data_load_path = datapath_base + config + "/"

data_save_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/mdotAsR/Constant" + dist + "/"
if not os.path.isdir(data_save_path):
    os.makedirs(data_save_path)
time_path = data_save_path + "mdot_for_time"
mdot_path = data_save_path + "mdot-data_allr"

radius_steps = np.arange(2, 4)
time_steps = np.arange(0, 3)
header = ""

for radius in radius_steps:
    for time in time_steps:
        timestep = "{:05d}".format(int(time))
        fname = datapath_base + "/" + config + ".prim." + timestep + ".athdf"
        quantities=['rho', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
        p = read.athdf(fname, quantities=quantities)
        r = p['x1v']
        theta = p['x2v']
        phi = p['x3v']
        R = p['x1f']
        nx1 = len(r)
        nx2 = len(theta)
        nx3 = len(phi)
        dp = (2*np.pi/nx3)*np.ones(nx3)
        dt = (np.pi/nx2)*np.ones(nx2)
        Dp = dp.reshape(nx3, 1, 1)
        Dt = dt.reshape(1, nx2, 1)
        dr = np.ones(nx1)
        for i in range(nx1):
            dr[i] = R[i+1]-R[i]
        Dr = dr.reshape(1, 1, nx1)
        rho = p['rho']
        ks = kerr(r, theta)
        f = fourvector(r, theta, p['vel1'], p['vel2'], p['vel3'], p['Bcc1'], p['Bcc2'], p['Bcc3'])
        k = Dt*Dp*ks.g
        Mflux = f.u1*k*rho
        total_mflux = np.sum(np.sum(Mflux, axis=0), axis=0)
        code_time = np.array(p['Time']).reshape(1,)
        time_data = np.array([code_time, total_mflux]).reshape(1, 2)
        print("time data:")
        print(time_data)
        if not os.path.isfile(time_path):
            header = "time, mdot"
        with open(time_path, "a") as f:
            np.savetxt(f, time_data, header=header, fmt='%s')
    #perform average calculations
    with open(time_path, "r") as f:
        mdot_avg = np.mean(f, axis=1)
        print(mdot_avg)

# clean up mdot by removing duplicates and sorting
with open(mdot_path, "r") as f:
    mdotdata = np.loadtxt(f, skiprows=1)

unique_times, inds = np.unique(mdotdata[:,0], return_index=True)
new_mdotdata = np.transpose(np.stack([unique_times, mdotdata[:,1][inds]]))
with open(mdot_path, "w") as f:
    np.savetxt(f, new_mdotdata, header="time, mdot")