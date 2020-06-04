"""
This is the most useful forces script - mag tension is transformed to local frame
script is mostly controlled by boolean operators
calculates and plots the forces at each timestep
INPUTS:
    - Times to look at
    - which forces to calculate and plot - mag, gas or total
    - with lines at inner radius of torus and at max density or not
    - limit plot to torus or not
OUTPUTS:
    - Plot with specified forces
"""

# import packages
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
sys.path.append("C:/Users/Ellie/Downloads/nerd/scripts/modules/")
import new_athena_read

# specifications
times_to_look_at = np.arange(206, 671)
mag_force = False
nonmag_force = False
total_force = True
lines = True
specific = False

# path to load data - save paths at end of script
datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
configA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2"
configB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2"
datapath_baseA = datapath + configA
datapath_baseB = datapath + configB

# dictionary for quantities necessary
quantities = ['press', 'Bcc1', 'Bcc2', 'Bcc3']
quantity_names = {"rho":"Density", "press":"Pressure", "vel1":"Radial velocity", "vel2":"Theta velocity",
                  "vel3":"Azimuthal velocity", "Bcc1":"Radial magnetic field", "Bcc2":"Theta magnetic field",
                  "Bcc3":"Azimuthal magnetic field"}

#calculate forces at local frame for each timestep
for timestep in times_to_look_at:
    # get data from file
    timestep = "{:05d}".format(int(timestep))
    filepathA = datapath_baseA + "/" + configA + ".prim." + timestep + ".athdf"
    filepathB = datapath_baseB + "/" + configB + ".prim." + timestep + ".athdf"
    print("Loading time step {}".format(timestep))
    dataA = new_athena_read.athdf(filepathA, quantities=quantities)
    dataB = new_athena_read.athdf(filepathB, quantities=quantities)

    # Get variables from data
    simulation_timeA = dataA["Time"]
    simulation_timeB = dataB["Time"]
    r = dataA["x1v"]
    theta = dataA["x2v"]
    pressdataA = dataA['press']
    pressdataB = dataB['press']
    pmagdataA = 0.5 * (dataA['Bcc1'] ** 2 + dataA['Bcc2'] ** 2 + dataA['Bcc3'] ** 2)
    pmagdataB = 0.5 * (dataB['Bcc1'] ** 2 + dataB['Bcc2'] ** 2 + dataB['Bcc3'] ** 2)

    # get line out at constant radius
    phi_index = 0; radius_index = 0; theta_indexA = int(pressdataA.shape[1] / 2)
    phi_index = 0; radius_index = 0; theta_indexB = int(pressdataB.shape[1] / 2)

    # calculate gas pressure
    if nonmag_force:
        drpressA = np.gradient(dataA['press'], r, edge_order=2, axis=2)
        drpressB = np.gradient(dataB['press'], r, edge_order=2, axis=2)

        pressforceA = - 0.5 * (drpressA[phi_index, theta_indexA, :] + drpressA[phi_index, theta_indexA - 1, :])
        pressforceB = - 0.5 * (drpressB[phi_index, theta_indexB, :] + drpressB[phi_index, theta_indexB - 1, :])

    if mag_force:
        # calculate magnetic pressure
        pmag_radial_dataA = 0.5 * (pmagdataA[phi_index, theta_indexA, :] + pmagdataA[phi_index, theta_indexA - 1, :])
        pmag_radial_dataB = 0.5 * (pmagdataB[phi_index, theta_indexB, :] + pmagdataB[phi_index, theta_indexB - 1, :])

        Bcc1_radial_dataA = dataA['Bcc1'][phi_index, theta_indexA, :]
        Bcc1_radial_dataB = dataB['Bcc1'][phi_index, theta_indexB, :]
        Bcc2_radial_dataA = dataA['Bcc2'][phi_index, theta_indexA, :]
        Bcc2_radial_dataB = dataB['Bcc2'][phi_index, theta_indexB, :]
        Bcc3_radial_dataA = dataA['Bcc3'][phi_index, theta_indexA, :]
        Bcc3_radial_dataB = dataB['Bcc3'][phi_index, theta_indexB, :]

        pmagforceA = -np.gradient(pmag_radial_dataA, dataA['x1v'], edge_order=2)
        pmagforceB = -np.gradient(pmag_radial_dataB, dataB['x1v'], edge_order=2)

        # set up variables for magnetic tension
        drBrA = np.gradient(Bcc1_radial_dataA, dataA['x1v'], edge_order=2) #derivative of Bcc1 with respect to radius
        drBrB = np.gradient(Bcc1_radial_dataB, dataB['x1v'], edge_order=2) #derivative of Bcc1 with respect to radius
        dthetaBrA = np.gradient(dataA['Bcc1'], dataA['x2v'], edge_order=2, axis=1) #derivative of Bcc1 with respect to theta
        dthetaBrB = np.gradient(dataB['Bcc1'], dataB['x2v'], edge_order=2, axis=1) #derivative of Bcc1 with respect to theta
        dphiBrA = np.gradient(dataA['Bcc1'], dataA['x3v'], edge_order=2, axis=0) #derivative of Bcc1 with respect to phi
        dphiBrB = np.gradient(dataB['Bcc1'], dataB['x3v'], edge_order=2, axis=0) #derivative of Bcc1 with respect to phi

        drBr_radial_dataA = drBrA
        drBr_radial_dataB = drBrB
        dthetaBr_radial_dataA = dthetaBrA[phi_index, theta_indexA, :]
        dthetaBr_radial_dataB = dthetaBrB[phi_index, theta_indexB, :]
        dphiBr_radial_dataA = dphiBrA[phi_index, theta_indexA, :]
        dphiBr_radial_dataB = dphiBrB[phi_index, theta_indexB, :]

        #big ass equation (magnetic tension)
        magtensionA = Bcc1_radial_dataA * drBr_radial_dataA\
                      + np.divide(Bcc2_radial_dataA * dthetaBr_radial_dataA +
                                  Bcc3_radial_dataA * dphiBr_radial_dataA -
                                  Bcc2_radial_dataA**2 - Bcc3_radial_dataA**2, dataA['x1v'])

        total_mag_forceA = magtensionA + pmagforceA

        magtensionB = Bcc1_radial_dataB * drBr_radial_dataB \
                     + np.divide(Bcc2_radial_dataB * dthetaBr_radial_dataB +
                                 Bcc3_radial_dataB * dphiBr_radial_dataB -
                                 Bcc2_radial_dataB ** 2 - Bcc3_radial_dataB ** 2, dataB['x1v'])

        total_mag_forceB = magtensionB + pmagforceB

    # Repeats of the above calculations but all together
    if total_force:
        # gas force
        drpressA = np.gradient(dataA['press'], dataA['x1v'], edge_order=2, axis=2)
        drpressB = np.gradient(dataB['press'], dataB['x1v'], edge_order=2, axis=2)

        pressforceA = - 0.5 * (drpressA[phi_index, theta_indexA, :] + drpressA[phi_index, theta_indexA - 1, :])
        pressforceB = - 0.5 * (drpressB[phi_index, theta_indexB, :] + drpressB[phi_index, theta_indexB - 1, :])

        # mag forces
        pmag_radial_dataA = 0.5 * (pmagdataA[phi_index, theta_indexA, :] + pmagdataA[phi_index, theta_indexA - 1, :])
        pmag_radial_dataB = 0.5 * (pmagdataB[phi_index, theta_indexB, :] + pmagdataB[phi_index, theta_indexB - 1, :])

        Bcc1_radial_dataA = dataA['Bcc1'][phi_index, theta_indexA, :]
        Bcc1_radial_dataB = dataB['Bcc1'][phi_index, theta_indexB, :]
        Bcc2_radial_dataA = dataA['Bcc2'][phi_index, theta_indexA, :]
        Bcc2_radial_dataB = dataB['Bcc2'][phi_index, theta_indexB, :]
        Bcc3_radial_dataA = dataA['Bcc3'][phi_index, theta_indexA, :]
        Bcc3_radial_dataB = dataB['Bcc3'][phi_index, theta_indexB, :]

        pmagforceA = -np.gradient(pmag_radial_dataA, dataA['x1v'], edge_order=2)
        pmagforceB = -np.gradient(pmag_radial_dataB, dataB['x1v'], edge_order=2)

        drBrA = np.gradient(Bcc1_radial_dataA, dataA['x1v'], edge_order=2)  # derivative of Bcc1 with respect to radius
        drBrB = np.gradient(Bcc1_radial_dataB, dataB['x1v'], edge_order=2)  # derivative of Bcc1 with respect to radius
        dthetaBrA = np.gradient(dataA['Bcc1'], dataA['x2v'], edge_order=2,
                                axis=1)  # derivative of Bcc1 with respect to theta
        dthetaBrB = np.gradient(dataB['Bcc1'], dataB['x2v'], edge_order=2,
                                axis=1)  # derivative of Bcc1 with respect to theta
        dphiBrA = np.gradient(dataA['Bcc1'], dataA['x3v'], edge_order=2,
                              axis=0)  # derivative of Bcc1 with respect to phi
        dphiBrB = np.gradient(dataB['Bcc1'], dataB['x3v'], edge_order=2,
                              axis=0)  # derivative of Bcc1 with respect to phi

        drBr_radial_dataA = drBrA
        drBr_radial_dataB = drBrB
        dthetaBr_radial_dataA = dthetaBrA[phi_index, theta_indexA, :]
        dthetaBr_radial_dataB = dthetaBrB[phi_index, theta_indexB, :]
        dphiBr_radial_dataA = dphiBrA[phi_index, theta_indexA, :]
        dphiBr_radial_dataB = dphiBrB[phi_index, theta_indexB, :]

        # big ass equation (magnetic tension)
        magtensionA = Bcc1_radial_dataA * drBr_radial_dataA \
                      + np.divide(Bcc2_radial_dataA * dthetaBr_radial_dataA +
                                  Bcc3_radial_dataA * dphiBr_radial_dataA -
                                  Bcc2_radial_dataA ** 2 - Bcc3_radial_dataA ** 2, dataA['x1v'])

        total_mag_forceA = magtensionA + pmagforceA

        magtensionB = Bcc1_radial_dataB * drBr_radial_dataB \
                      + np.divide(Bcc2_radial_dataB * dthetaBr_radial_dataB +
                                  Bcc3_radial_dataB * dphiBr_radial_dataB -
                                  Bcc2_radial_dataB ** 2 - Bcc3_radial_dataB ** 2, dataB['x1v'])

        total_mag_forceB = magtensionB + pmagforceB

        total_forceA = magtensionA + pressforceA + pmagforceA
        total_forceB = magtensionB + pressforceB + pmagforceB

    #plot constant Beta
    plt.figure()
    if mag_force:
        plt.plot(dataA['x1v'], pmagforceA, label ="Magnetic Pressure Force", linestyle="--", color='red')
        plt.plot(dataA['x1v'], magtensionA, label="Magnetic Tension Force", linestyle="--", color='purple')
        plt.plot(dataA['x1v'], total_mag_forceA, label='Total Magnetic Force', linestyle=":", color='black')
        if lines:
            plt.axvline(5, -100, 100, label='Inner Radius', linestyle='dotted', color='green')
            plt.axvline(10, -100, 100, label='Highest Density', linestyle='dotted', color='green')
    if nonmag_force:
        plt.plot(dataA['x1v'], pressforceA, label = "Gas Pressure Force", linestyle="--", color='blue')
        if lines:
            plt.axvline(5, -100, 100, label='Inner Radius', linestyle='dotted', color='green')
            plt.axvline(10, -100, 100, label='Highest Density', linestyle='dotted', color='green')
    if total_force:
        plt.plot(dataA['x1v'], pmagforceA, label="Magnetic Pressure Force", linestyle="--", color='red')
        plt.plot(dataA['x1v'], magtensionA, label="Magnetic Tension Force", linestyle="--", color='purple')
        plt.plot(dataA['x1v'], total_mag_forceA, label='Total Magnetic Force', linestyle=":", color='black')
        plt.plot(dataA['x1v'], pressforceA, label="Gas Pressure Force", linestyle="--", color='blue')
        plt.plot(dataA['x1v'], total_forceA, label='Total Force', linestyle="-.", color='black')
        if lines:
            plt.axvline(5, -100, 100, label='Inner Radius', linestyle='dotted', color='green')
            plt.axvline(10, -100, 100, label='Highest Density', linestyle='dotted', color='green')
    if specific:
        plt.xlim(5, 20)
    plt.xlabel("Radius")
    plt.ylabel("Radial Forces")
    plt.title("Time: {}".format(simulation_timeA) + "\nConstant Beta")
    plt.tight_layout()
    plt.legend()

    if mag_force:
        if specific:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/forces_rprofile_local/specific/Beta/mag/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        else:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/forces_rprofile_local/Beta/mag/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        # plt.show()
        plt.gca().clear()
        plt.close()
    if nonmag_force:
        if specific:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/forces_rprofile_local/Beta/gas/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        else:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/forces_rprofile_local/specific/Beta/gas/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        # plt.show()
        plt.gca().clear()
        plt.close()
    if total_force:
        if specific:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/forces_rprofile_local/specific/Beta/total/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        else:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/forces_rprofile_local/Beta/total/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        #plt.show()
        plt.gca().clear()
        plt.close()


    #plot constant B
    plt.figure()
    if mag_force:
        plt.plot(dataB['x1v'], pmagforceB, label="Magnetic Pressure Force", linestyle="--", color='red')
        plt.plot(dataB['x1v'], magtensionB, label="Magnetic Tension Force", linestyle="--", color='purple')
        plt.plot(dataB['x1v'], total_mag_forceB, label='Total Magnetic Force', linestyle=":", color='black')
        if lines:
            plt.axvline(5, -100, 100, label='Inner Radius', linestyle='dotted', color='green')
            plt.axvline(10, -100, 100, label='Highest Density', linestyle='dotted', color='green')
    if nonmag_force:
        plt.plot(dataB['x1v'], pressforceB, label="Gas Pressure Force", linestyle="--", color='blue')
        if lines:
            plt.axvline(5, -100, 100, label='Inner Radius', linestyle='dotted', color='green')
            plt.axvline(10, -100, 100, label='Highest Density', linestyle='dotted', color='green')
    if total_force:
        plt.plot(dataB['x1v'], pmagforceB, label="Magnetic Pressure Force", linestyle="--", color='red')
        plt.plot(dataB['x1v'], magtensionB, label="Magnetic Tension Force", linestyle="--", color='purple')
        plt.plot(dataB['x1v'], total_mag_forceB, label='Total Magnetic Force', linestyle=":", color='black')
        plt.plot(dataB['x1v'], pressforceB, label="Gas Pressure Force", linestyle="--", color='blue')
        plt.plot(dataB['x1v'], total_forceB, label='Total Force', linestyle="-.", color='black')
        if lines:
            plt.axvline(5, -100, 100, label='Inner Radius', linestyle='dotted', color='green')
            plt.axvline(10, -100, 100, label='Highest Density', linestyle='dotted', color='green')
    if specific:
        plt.xlim(5, 20)
    plt.xlabel("Radius")
    plt.ylabel("Radial Forces")
    plt.title("Time: {}".format(simulation_timeA) + "\nConstant B")
    plt.tight_layout()
    plt.legend()

    if mag_force:
        if specific:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/forces_rprofile_local/specific/B/mag/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        else:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/forces_rprofile_local/B/mag/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        #plt.show()
        plt.gca().clear()
        plt.close()
    if nonmag_force:
        if specific:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/forces_rprofile_local/specific/B/gas/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        else:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/forces_rprofile_local/B/gas/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        #plt.show()
        plt.gca().clear()
        plt.close()
    if total_force:
        if specific:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/forces_rprofile_local/specific/B/total/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        else:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSProfiles/forces_rprofile_local/B/total/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        #plt.show()
        plt.gca().clear()
        plt.close()