import argparse
import warnings
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ----- Input Parameters ------
# What to plot
radius = 10.0
quantity = "Bcc1"
m_to_plot = 1
# Which simulation to run
config = "1.1.1-torus2_b-gz2"
specs = "a0beta500torB_br32x32x64rl2x2"
specs = "a0beta500torBeta_br32x32x64rl2x2"
# where to get data from
# Ellie will need to set data_save_path to datapath, fourier_path to path
data_save_path = "/run/media/amha5924/data/" + config + "_" + specs + "/data/reduced/"
fourier_path = "fourier-amplitude-" + quantity + "_data_r{}.txt".format(radius)



# ---- Some constants that don't need to be modified
# (except perhaps once, e.g. to change figsavedir)
quantity_names = {"rho": "Density", "press": "Pressure", "vel1": "Radial velocity", "vel2": "Theta velocity", "vel3": "Azimuthal velocity", "Bcc1": "Radial magnetic field", "Bcc2": "Theta magnetic field", "Bcc3": "Azimuthal magnetic field"}
specstr = {"a0beta500torB_br32x32x64rl2x2":"Constant Magnetic Field B","a0beta500torBeta_br32x32x64rl2x2":"Constant Beta"}
specname = {"a0beta500torB_br32x32x64rl2x2":"ConstantB","a0beta500torBeta_br32x32x64rl2x2":"ConstantBeta"}

# where to save figures
fig_save_dir = "/home/amha5924/research/projects/gr-torus/figures/fourier_amplitude/" + quantity + "/" + specname[specs] + "/"
# figsavedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/mode_amplitude/ConstantBeta/"
fig_name = "m{}_r{}.png".format(m_to_plot, radius)

def exponential(norm_times, amp, gamma, offset):
    return amp * np.exp(gamma * norm_times) + offset
# --------------------------------------------

# Load the data from file
with open(data_save_path + fourier_path) as file:
    data = np.loadtxt(file, skiprows=1)

# Normalize times to be between 0 and 1
times = data[:, 0]
norm_times = times / times[-1]
m_to_plot_amp = data[:, 1 + m_to_plot]

# Fit m amplitude to an exponential
param, param_cov = curve_fit(exponential, norm_times, m_to_plot_amp)

print("Exponent coefficient: ")
print(param)
print("Covariance: ")
print(param_cov)

# Create best fit line
best_fit_norm = (param[0] * np.exp(param[1] * norm_times))
normalized_gamma = param[1]
actual_gamma = normalized_gamma/times[-1]
best_fit = (param[0] * np.exp(actual_gamma * times))

# Plot data and best fit line
plt.plot(times, m_to_plot_amp, marker='+', ls="", color='blue', label='Data')
plt.plot(times, best_fit, ls="--", color='red', label='Best Fit')
# plt.plot(norm_times, m_to_plot_amp, marker='+', ls="", color='blue', label='Data')
# plt.plot(norm_times, best_fit_norm, ls="--", color='red', label='Best Fit')
plt.legend()
plt.xlabel("Time [GM/c^3]")
plt.ylabel("Mode Amplitude [arbitrary units]")
# Create title string that contains all information about the plot
# This should be cropped when putting into papers/presentations
# (the info will be in the caption instead of the title) but it is
# nice to see it on the plot for reference
title_str = specstr[specs]
title_str += "\n Mode {} Amplitude".format(m_to_plot)
title_str += "\n Radius: {}".format(radius)
title_str += "\n " + quantity_names[quantity]
title_str += "\n Normalized best fit growth rate: {:.2f}".format(normalized_gamma)
title_str += "\n Actual best fit growth rate: {:.2E} [c^3/GM]".format(actual_gamma)
plt.title(title_str)
plt.tight_layout()
if not os.path.isdir(fig_save_dir):
    os.mkdir(fig_save_dir)
print("Saving figure " + fig_save_dir + fig_name)
plt.savefig(fig_save_dir + fig_name)
plt.show()
