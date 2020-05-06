import numpy as np
import matplotlib.pyplot as plt

datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
configB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2"
configA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2"
datapath_baseA = datapath + configA
datapath_baseB = datapath + configB

filenameA = datapath_baseA + "/" + configA + ".hst"
filenameB = datapath_baseB + "/" + configB + ".hst"

mean_values: object = hst.data[:, 12]
time_values = hst.data[:,0]
for i in range (1, time_values)
    time = time_values[i]
    mean_values.step_str: str = "{:05d}".format(i)
    data = load_athd(filenameA + step_str + "athdf")
    plot(data["rho"]["phi", "theta", :]**2 - mean_values[i]