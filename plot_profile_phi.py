# Other Python modules
import numpy as np
import matplotlib.pyplot as plt
import os

# Athena++ modules
from scripts import athena_read

datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
configA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2"
configB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2"
datapath_baseA = datapath + configA
datapath_baseB = datapath + configB

# Timesteps NOT time numbers, this is what shows up in the file name
# np.arange(start, stop, step) for even steps
# np.linspace(start, stop, number of samples*, endpoint (t/f)*) for evenly spaced steps over the interval
times_to_look_at = np.arange(0, 10)
#times_to_look_at = np.arange(0, 200, 100)
#times_to_look_at = np.concatenate((np.arange(0, 26), np.arange(27, 448)))

# dictionary, input "rho" or "Bcc1" or smth, get out info in the corresponding value
quantities = ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
quantity_names = {"rho":"Density", "press":"Pressure", "vel1":"Radial velocity", "vel2":"Theta velocity",
                  "vel3":"Azimuthal velocity", "Bcc1":"Radial magnetic field", "Bcc2":"Theta magnetic field",
                  "Bcc3":"Azimuthal magnetic field"}
quantity_to_load = "Bcc2"

for timestep in times_to_look_at:
    timestep = "{:05d}".format(int(timestep))
    filepathA = datapath_baseA + "/" + configA + ".prim." + timestep + ".athdf"
    filepathB = datapath_baseB + "/" + configB + ".prim." + timestep + ".athdf"
    print("Loading time step {}".format(timestep))
    dataA = athena_read.athdf(filepathA, quantities=[quantity_to_load])
    dataB = athena_read.athdf(filepathB, quantities=[quantity_to_load])

    # print(data.keys())
    # print("Simulation time is {}".format(data["Time"]))

    #"Time" is how that variable is titled in the data, so the capital is important
    simulation_timeA = dataA["Time"]
    simulation_timeB = dataB["Time"]

    quantity_dataA = dataA[quantity_to_load]
    quantity_dataB = dataB[quantity_to_load]
    phi_valuesA = dataA['x3v']
    phi_valuesB = dataB['x3v']
    #print(phi_valuesA[int(phi_valuesA.size)])

    radius_in_codeunits = 2
    radiusind = (np.abs(dataA['x1v'] - radius_in_codeunits)).argmin()
    # print(radiusind)

    # get line out at constant radius
    phi_index = 0; radius_index = radiusind; theta_indexA = int(quantity_dataA.shape[1]/2)
    phi_index = 0; radius_index = radiusind; theta_indexB = int(quantity_dataB.shape[1]/2)
    phi_dataA = quantity_dataA[:, theta_indexA, radius_index]
    phi_dataB = quantity_dataB[:, theta_indexB, radius_index]
    plt.plot(phi_valuesA, phi_dataA, label = "Constant Beta", linestyle = "--")
    # plt.plot(phi_valuesB, phi_dataB, label = "Constant B")
    plt.xlabel("Phi")
    plt.ylabel(quantity_names[quantity_to_load])
    plt.title("Time: {}".format(simulation_timeA))
    plt.legend()
    plt.tight_layout()

    figname = quantity_to_load + "_at_timestep_{}".format(timestep)
    filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/" + quantity_to_load + "_r2phiprofile/"
    if not os.path.isdir(filedir):
        os.mkdir(filedir)
    print(filedir)
    print("Saving figure " + filedir + figname)
    plt.savefig(filedir + figname)
    plt.show()
    plt.gca().clear()