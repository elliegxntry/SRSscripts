"""
This script calculates the total mass loss from the .hst file, and compares it with the
mass loss at 3 points in the simulation.
INPUTS:
    - mdot data from calculate_mdot.py (1.52, 20, 40)
    - .hst file
    - initial distribution
    - timesteps to calculate
    -
OUTPUTS:
    - File with mass loses
"""

# import packages
import numpy as np
import os

# specifications
dist = "B"

# path bases
path_to_raw_data = "C:/Users/Ellie/Downloads/nerd/SRSData/"
path_to_reduced_data = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/"
config = "1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2"
filedir_base = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/"

# load mdot data
mdot_inner_edge_data = np.loadtxt(path_to_reduced_data + config + "/mdot-data_r1.52.txt")
mdot_outer_torus_data = np.loadtxt(path_to_reduced_data + config + "/mdot-data_r20.txt")
mdot_outer_edge_data = np.loadtxt(path_to_reduced_data + config + "/mdot-data_r40.txt")

# load .hst data
hst_data = np.loadtxt(path_to_raw_data + config + "/" + config + ".hst", skiprows=2)
time = hst_data[:, 0]
total_mass = hst_data[:, 1]

# net mass loss
for mass_step in total_mass:
    change_in_mass = total_mass[mass_step + 1] - total_mass[mass_step]