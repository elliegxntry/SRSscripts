# import Python modules
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append("C:/Users/Ellie/Downloads/nerd/scripts/modules/")
import new_athena_read

datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
config = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2"
datapath_base = datapath + config

times_to_look_at = np.arange(0, 600, 100)
quantities_to_load = ['vel2', 'Bcc2']

# dictionary for quantities
quantities = ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
quantity_names = {"rho":"Density", "press":"Pressure", "vel1":"Radial velocity", "vel2":"Theta velocity",
                  "vel3":"Azimuthal velocity", "Bcc1":"Radial magnetic field", "Bcc2":"Theta magnetic field",
                  "Bcc3":"Azimuthal magnetic field"}

for timestep in times_to_look_at:
    timestep = "{:05d}".format(int(timestep))
    filepath = datapath_base + "/" + config + ".prim." + timestep + ".athdf"
    print("Loading time step {}".format(timestep))
    data = athena_read.athdf(filepath, quantities=quantities_to_load)

    # get variables
    simulation_time = data["Time"]
    theta_values = data['x2v']
    #print(theta_values[int(theta_values.size/2)])

    radius_in_codeunits = 4
    radiusind = (np.abs(data['x1v'] - radius_in_codeunits)).argmin()
    #print(radiusind)

    # get line out at variable theta
    phi_index = 0; radius_index = radiusind; theta_indexA = int(theta_values.size/2)
    vel2 = data['vel2'][phi_index, :, radius_index]
    Bcc2 = data['Bcc2'][phi_index, :, radius_index]
    plt.plot(theta_values, vel2/np.max(np.abs(vel2)), label = "Theta velocity", linestyle = "--")
    plt.plot(theta_values, Bcc2/np.max(np.abs(Bcc2)), label = "Theta Magnetic field")
    plt.xlabel("Theta")
    plt.xlim([np.pi/2-0.1, np.pi/2+0.1])
    plt.title("Time: {}\n radius: {}\n".format(simulation_time, radius_in_codeunits) + config)
    plt.legend()
    plt.tight_layout()

    figname = "v2b2corr_r{}".format(radius_in_codeunits) + "_at_timestep_{}".format(timestep)
    filedir ="C:/Users/Ellie/Downloads/nerd/SRSProfiles/v2b2corr_thetaprofile_c01/"
    if not os.path.isdir(filedir):
        os.mkdir(filedir)
    #print(filedir)
    #print("Saving figure " + filedir + figname)
    plt.savefig(filedir + figname)
    plt.show()
    plt.gca().clear()

