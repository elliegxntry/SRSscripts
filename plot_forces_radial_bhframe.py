# Other Python modules
import numpy as np
import matplotlib.pyplot as plt
import os
import kerrmetric as kerr

# Athena++ modules
from scripts import athena_read

times_to_look_at = np.arange(0, 671)

datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
configA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2"
configB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2"
datapath_baseA = datapath + configA
datapath_baseB = datapath + configB


# dictionary, input "rho" or "Bcc1" or smth, get out info in the corresponding value
quantities = ['press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
quantity_names = {"rho":"Density", "press":"Pressure", "vel1":"Radial velocity", "vel2":"Theta velocity",
                  "vel3":"Azimuthal velocity", "Bcc1":"Radial magnetic field", "Bcc2":"Theta magnetic field",
                  "Bcc3":"Azimuthal magnetic field"}

for timestep in times_to_look_at:
    timestep = "{:05d}".format(int(timestep))
    filepathA = datapath_baseA + "/" + configA + ".prim." + timestep + ".athdf"
    filepathB = datapath_baseB + "/" + configB + ".prim." + timestep + ".athdf"
    print("Loading time step {}".format(timestep))
    dataA = athena_read.athdf(filepathA, quantities=quantities)
    dataB = athena_read.athdf(filepathB, quantities=quantities)

    # print(data.keys())
    # print("Simulation time is {}".format(data["Time"]))


# "Time" is how that variable is titled in the data, so the capital is important
    simulation_timeA = dataA["Time"]
    simulation_timeB = dataB["Time"]

    pressdataA = dataA['press']
    pressdataB = dataB['press']
    r = dataA["x1v"]
    theta = dataA["x2v"]
    metric = kerr.kerr(r, theta)
    fourvectorA = kerr.fourvector(r, theta, dataA['vel1'], dataA['vel2'], dataA['vel3'], dataA['Bcc1'], dataA['Bcc2'],
                                 dataA['Bcc3'])
    pmagdataA = fourvectorA.pmag
    fourvectorB = kerr.fourvector(r, theta, dataB['vel1'], dataB['vel2'], dataB['vel3'], dataB['Bcc1'], dataB['Bcc2'],
                                  dataB['Bcc3'])
    pmagdataB = fourvectorB.pmag

    radius_in_codeunits = 5
    radiusind = (np.abs(dataA['x1v'] - radius_in_codeunits)).argmin()
   # print(radiusind)

    # get line out at constant radius
    phi_index = 0; radius_index = 0; theta_indexA = int(pressdataA.shape[1]/2)
    phi_index = 0; radius_index = 0; theta_indexB = int(pressdataB.shape[1]/2)
    press_radial_dataA = pressdataA[phi_index, theta_indexA, :]
    press_radial_dataB = pressdataB[phi_index, theta_indexB, :]
    pmag_radial_dataA = pmagdataA[phi_index, theta_indexA, :]
    pmag_radial_dataB = pmagdataB[phi_index, theta_indexB, :]
    pressforceA = -np.gradient(press_radial_dataA, dataA['x1v'], edge_order=2)
    pressforceB = -np.gradient(press_radial_dataB, dataA['x1v'], edge_order=2)
    pmagforceA = -np.gradient(pmag_radial_dataA, dataA['x1v'], edge_order=2)
    pmagforceB = -np.gradient(pmag_radial_dataB, dataB['x1v'], edge_order=2)


    #plot constant Beta
    plt.plot(dataA['x1v'], pressforceA, label = "Gas Pressure Force", linestyle = "--", color='blue')
    plt.plot(dataA['x1v'], pmagforceA, label = "Magnetic Pressure Force", linestyle = "--", color='red')
    plt.xlabel("Radius")
    plt.ylabel("Radial Forces")
    plt.title("Time: {}".format(simulation_timeA))
    plt.tight_layout()
    plt.legend()


    figname = "forces_at_timestep_{}".format(timestep)
    filedir ="C:/Users/Ellie/Downloads/nerd/SRSProfiles/forces_rprofile/Beta/"
    if not os.path.isdir(filedir):
        os.mkdir(filedir)
    #print(filedir)
    #print("Saving figure " + filedir + figname)
    plt.savefig(filedir + figname)
    #plt.show()
    plt.gca().clear()

    #plot constant B
    plt.plot(dataB['x1v'], pressforceB, label = "Gas Pressure Force", linestyle = "--", color='blue')
    plt.plot(dataB['x1v'], pmagforceB, label = "Magnetic Pressure Force", linestyle = "--", color='red')
    plt.xlabel("Radius")
    plt.ylabel("Radial Forces")
    plt.title("Time: {}".format(simulation_timeA))
    plt.tight_layout()
    plt.legend()


    figname = "forces_at_timestep_{}".format(timestep)
    filedir ="C:/Users/Ellie/Downloads/nerd/SRSProfiles/forces_rprofile/B/"
    if not os.path.isdir(filedir):
        os.mkdir(filedir)
    #print(filedir)
    #print("Saving figure " + filedir + figname)
    plt.savefig(filedir + figname)
    #plt.show()
    plt.gca().clear()