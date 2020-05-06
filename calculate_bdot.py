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
dist = "B"
radius = 1.52
time_steps = np.arange(0, 671)

do_Bcc1 = False
do_Bcc2 = False
do_Bcc3 = True

config = "1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2"
datapath_base = datapath + config
if do_Bcc1:
    data_save_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/Constant" + dist + "/Bcc1/"
    if not os.path.isdir(data_save_path):
        os.makedirs(data_save_path)
if do_Bcc2:
    data_save_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/Constant" + dist + "/Bcc2/"
    if not os.path.isdir(data_save_path):
        os.makedirs(data_save_path)
if do_Bcc3:
    data_save_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/Constant" + dist + "/Bcc3/"
    if not os.path.isdir(data_save_path):
        os.makedirs(data_save_path)

data_load_path = datapath_base + config + "/"
bdot_path = data_save_path + "bdot-data_r{}.txt".format(radius)

header = ""

for time in time_steps:
    timestep = "{:05d}".format(int(time))
    fname = datapath_base + "/" + config + ".prim." + timestep + ".athdf"
    quantities=['rho','vel1','vel2','vel3','Bcc1','Bcc2','Bcc3']
    print("Loading time {}".format(time))
    p = read.athdf(fname, quantities=quantities, x1_min=radius, x1_max=radius)
    r = p['x1v']
    theta = p['x2v']
    phi = p['x3v']
    R = p['x1f']
    nx1 = len(r)
    nx2 = len(theta)
    nx3 = len(phi)
    #what is p in dp?
    dp = (2*np.pi/nx3)*np.ones(nx3)
    dt = (np.pi/nx2)*np.ones(nx2)
    Dp = dp.reshape(nx3,1,1)
    Dt = dt.reshape(1,nx2,1)
    dr = np.ones(nx1)
    for i in range(nx1):
        dr[i] = R[i+1]-R[i]
    Dr = dr.reshape(1,1,nx1)
    rho = p['rho']
    Bcc1 = p['Bcc1']
    Bcc2 = p['Bcc2']
    Bcc3 = p['Bcc3']
    ks = kerr(r,theta)
    f = fourvector(r,theta,p['vel1'],p['vel2'],p['vel3'],p['Bcc1'],p['Bcc2'],p['Bcc3'])
    k = Dt*Dp*ks.g
    if do_Bcc1:
        Bflux = f.u1*k*Bcc1
    if do_Bcc2:
        Bflux = f.u1*k*Bcc2
    if do_Bcc3:
        Bflux = f.u1*k*Bcc3
    total_bflux = np.sum(np.sum(Bflux, axis=0), axis=0)
    code_time = np.array(p['Time']).reshape(1,)
    output_data = np.array([code_time, total_bflux]).reshape((1,2))
    if not os.path.isfile(bdot_path):
        header = "time, bdot"
    with open(bdot_path, "a") as f:
        np.savetxt(f, output_data, header=header)

# clean up bdot by removing duplicates and sorting
with open(bdot_path, "r") as f:
    bdotdata = np.loadtxt(f, skiprows=1)

unique_times, inds = np.unique(bdotdata[:, 0], return_index=True)
new_bdotdata = np.transpose(np.stack([unique_times, bdotdata[:, 1][inds]]))
with open(bdot_path, "w") as f:
    np.savetxt(f, new_bdotdata, header="time, bdot")
