"""
Extracts growth rate of fourier transform from mri_modes.py and fits an exponential curve to the data
INPUTS:
    - Quantity (usually a bfield)
    - radius
    - mode to calculate and plot
OUTPUTS:
    - Plot of the data with a curve fit
    - prints growth rate and covariance (in python and in title of plot)

If line won't work/there's an error in size, try:
    - adding a p0 variable (from print(param)) to the curve_fit function from a similar plot
    - limit times
    - add more here - google
"""

# Import packages
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import sys
sys.path.append("C:/Users/Ellie/Downloads/nerd/scripts/modules")

# Quantity dictionary
quantities = ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
quantity_names = {"rho": "Density", "press": "Pressure", "vel1": "Radial velocity", "vel2": "Theta velocity",
                  "vel3": "Azimuthal velocity", "Bcc1": "Radial magnetic field", "Bcc2": "Theta magnetic field",
                  "Bcc3": "Azimuthal magnetic field"}

# Specifications
quantity = "Bcc3"
radius = 1.52
m = 2 #mode
dist = "B"
start = 0
end = -1

# Paths to load and save
datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/Constant" + dist + "/Fourier/"
fourier_path = "fourier-amplitude-" + quantity + "_data_r{}.txt".format(radius)
figsavepath = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/mode_amplitudes/Constant" + dist + "/m" + str(m) + "/"
figname = quantity + "m{}_r{}.png".format(m, radius)

# Load fourier data
with open(datapath + fourier_path) as file:
    data = np.loadtxt(file, skiprows=1)

# Define where to load data from in array and normalize times
times = (data[:, 0])[start:end]
norm_times = times / times[-1]
const_amp = data[:, 1]
m1_amp = data[:, 2]
mN_amp = (data[:, m + 1])[start:end]

# Create exponential function
def exponential(norm_times, amp, gamma, shift, offset):
    return amp * np.exp(gamma * norm_times - shift) + offset

# Using the curve_fit function, calculate the equation for the exponential growth
    # param is gamma/the rate of growth, param_cov is the accuracy/covariance of the estimate
param, param_cov = curve_fit(exponential, norm_times, mN_amp, p0=[1.32786845e-11, 1.69920332e+01, 6.37176322e+00, 1.87400919e-09])
print("amp = " + str(param[0]))
print("shift = " + str(param[2]))
print("offset = " + str(param[3]))
#print("Exponent coefficient: ")
print(param)
#print("Covariance: ")
#print(param_cov)

# define growth rate from calculated equation
norm_fit = (param[0] * np.exp(param[1] * norm_times - param[2])) + param[3]
norm_gamma = param[1]
gamma = norm_gamma/times[-1]
fit = param[0] * np.exp(gamma * times - param[2]) + param[3]
print("Normalized gamma:")
print(norm_gamma)
print("Gamma: ")
print(gamma)

# Plot curve fit
plt.plot(times, mN_amp, marker='+', ls="", color='blue', label='Data')
plt.plot(times, fit, ls="--", color='red', label='Best Fit')
# plt.plot(norm_times, mN_amp, marker='+', ls="", color='blue', label='Data')
# plt.plot(norm_times, fit, ls="--", color='red', label='Best Fit')
plt.legend()
title = "Constant " + dist + " Mode m{} Amplitude".format(m) + "\nRadius: {}".format(radius) + "\n" + quantity_names[quantity]
values = "\nNormalized best fit growth rate: {:.2f}".format(norm_gamma) + "\nActual best fit growth rate: {:.2E} [c^3/GM]".format(gamma)
plt.xlabel("Time [GM/c^3]")
plt.ylabel("Mode Amplitude [arbitrary units]")
plt.title(title + values)
if not os.path.isdir(figsavepath):
    os.mkdirs(figsavepath)
#print(figsavepath)
#print(figname)
plt.savefig(figsavepath + figname)
plt.tight_layout()
plt.show()
