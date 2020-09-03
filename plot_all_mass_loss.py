"""
This script calculates the mass loss over the whole disk (hst file),
 the mass loss over the inner edge, the edge of the torus and the
outer edge of the simulation
MUST CALCULATE MDOT IN calculate_mdot.py FIRST!!!!!!!!! OR ELSE IT WON'T WORK!!!! DON'T FORGET!!!!
INPUTS:
    -
OUTPUTS:
    -
"""

# import packages
import matplotlib.pyplot as plt
import numpy as np

# specifications
dist = "B"
#prolly some if statements

# path bases
path_to_raw_data = "C:/Users/Ellie/Downloads/nerd/SRSData/"
path_to_reduced_data = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/"
config = "1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2"
filedir_base = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/time_quantity_profiles/"

# load mdot data
mdot_inner_edge_data = np.loadtxt(path_to_reduced_data + config + "/mdot-data_r1.52.txt")
mdot_outer_torus_data = np.loadtxt(path_to_reduced_data + config + "/mdot-data_r20.txt")
mdot_outer_edge_data = np.loadtxt(path_to_reduced_data + config + "/mdot-data_r40.txt")

# load .hst data
hst_data = np.loadtxt(path_to_raw_data + config + "/" + config + ".hst", skiprows=2)
hst_time = hst_data[:, 0]
total_mass = hst_data[:, 1]

# call mdot data
mdot_inner_edge = mdot_inner_edge_data[:, 1]
mdot_outer_torus = mdot_outer_torus_data[:, 1]
mdot_outer_edge = mdot_outer_edge_data[:, 1]
time_inner_edge = mdot_inner_edge_data[:, 0]
time_outer_torus= mdot_outer_torus_data[:, 0]
time_outer_edge = mdot_outer_edge_data[:, 0]

title_string = "Mass loss comparison"

plt.figure()
plt.plot(time_inner_edge, mdot_inner_edge, color="pink", label='inner edge (r=1.52)')
plt.plot(time_outer_torus, mdot_outer_torus, color='blue', label='outer torus(r=20.0)')
plt.plot(time_outer_edge, mdot_outer_edge, color='purple', label='outer edge (r=40.0)')
plt.plot(hst_time, total_mass, color='red', label='total mass loss')
plt.legend()
plt.title(title_string + "\n Constant: " + dist)
plt.tight_layout()
figname = "mass_loss"
plt.show()
#plt.savefig(filedir_base + figname + "/" + dist)
plt.close()