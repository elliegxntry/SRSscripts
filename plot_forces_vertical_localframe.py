# Other Python modules
# import packages
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append("C:/Users/Ellie/downloads/nerd/scripts/modules")
import new_athena_read

# specifications
times_to_look_at = np.arange(0, 1)
mag_force = False
nonmag_force = False
total_force = True
lines = True
specific = True

# Path to load data
datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
configA = "1.1.1-torus2_b-gz2_a0beta500torBeta_br32x32x64rl2x2"
configB = "1.1.1-torus2_b-gz2_a0beta500torB_br32x32x64rl2x2"
datapath_baseA = datapath + configA
datapath_baseB = datapath + configB

# dictionary for quantities
quantities = ['press', 'Bcc1', 'Bcc2', 'Bcc3']
quantity_names = {"rho":"Density", "press":"Pressure", "vel1":"Radial velocity", "vel2":"Theta velocity",
                  "vel3":"Azimuthal velocity", "Bcc1":"Radial magnetic field", "Bcc2":"Theta magnetic field",
                  "Bcc3":"Azimuthal magnetic field"}

# Calculate forces
for timestep in times_to_look_at:
    # load data
    timestep = "{:05d}".format(int(timestep))
    filepathA = datapath_baseA + "/" + configA + ".prim." + timestep + ".athdf"
    filepathB = datapath_baseB + "/" + configB + ".prim." + timestep + ".athdf"
    #print("Loading time step {}".format(timestep))
    dataA = new_athena_read.athdf(filepathA, quantities=quantities)
    dataB = new_athena_read.athdf(filepathB, quantities=quantities)

    # define variables
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

    if nonmag_force:
        dthetapressA = np.gradient(dataA['press'], dataA['x2v'], edge_order=2, axis=1)
        dthetapressB = np.gradient(dataB['press'], dataB['x2v'], edge_order=2, axis=1)

        pressforceA = - 0.5 * (dthetapressA[phi_index, theta_indexA, :] + dthetapressA[phi_index, theta_indexA - 1, :])
        pressforceB = - 0.5 * (dthetapressB[phi_index, theta_indexB, :] + dthetapressB[phi_index, theta_indexB - 1, :])

    if mag_force:
        dthetapmagA = np.gradient(pmagdataA, dataA['x2v'], edge_order=2,
                                   axis=1)  # derivative of Bcc1 with respect to theta
        dthetapmagB = np.gradient(pmagdataB, dataB['x2v'], edge_order=2,
                                   axis=1)  # derivative of Bcc1 with respect to theta

        Bcc1_radial_dataA = dataA['Bcc1'][phi_index, theta_indexA, :]
        Bcc1_radial_dataB = dataB['Bcc1'][phi_index, theta_indexB, :]
        Bcc2_radial_dataA = dataA['Bcc2'][phi_index, theta_indexA, :]
        Bcc2_radial_dataB = dataB['Bcc2'][phi_index, theta_indexB, :]
        Bcc3_radial_dataA = dataA['Bcc3'][phi_index, theta_indexA, :]
        Bcc3_radial_dataB = dataB['Bcc3'][phi_index, theta_indexB, :]

        pmagforceA = dthetapmagA[phi_index, theta_indexA, :]
        pmagforceB = dthetapmagB[phi_index, theta_indexB, :]

        dthetaBthetaA = np.gradient(dataA['Bcc2'], dataA['x2v'], edge_order=2,
                                axis=1)  # derivative of Bcc1 with respect to theta
        dthetaBthetaB = np.gradient(dataB['Bcc2'], dataB['x2v'], edge_order=2,
                                axis=1)  # derivative of Bcc1 with respect to theta
        dphiBthetaA = np.gradient(dataA['Bcc2'], dataA['x3v'], edge_order=2,
                              axis=0)  # derivative of Bcc1 with respect to phi
        dphiBthetaB = np.gradient(dataB['Bcc2'], dataB['x3v'], edge_order=2,
                              axis=0)  # derivative of Bcc1 with respect to phi
        drBthetaA = np.gradient(dataA['Bcc2'], dataA['x1v'], edge_order=2, axis=2)
        drBthetaB = np.gradient(dataB['Bcc2'], dataB['x1v'], edge_order=2, axis=2)

        dthetaBtheta_radial_dataA = dthetaBthetaA[phi_index, theta_indexA, :]
        dthetaBtheta_radial_dataB = dthetaBthetaB[phi_index, theta_indexB, :]
        dphiBtheta_radial_dataA = dphiBthetaA[phi_index, theta_indexA, :]
        dphiBtheta_radial_dataB = dphiBthetaB[phi_index, theta_indexB, :]
        drBtheta_radial_dataA = drBthetaA[phi_index, theta_indexA, :]
        drBtheta_radial_dataB = drBthetaB[phi_index, theta_indexA, :]


        #big ass equation (magnetic tension)
        magtensionA = Bcc1_radial_dataA * drBtheta_radial_dataA + np.divide(Bcc2_radial_dataA * dthetaBtheta_radial_dataA +
                                  Bcc3_radial_dataA * dphiBtheta_radial_dataA +
                                  Bcc2_radial_dataA * Bcc1_radial_dataA, dataA['x1v'])

        total_mag_forceA = magtensionA + pmagforceA

        magtensionB = Bcc1_radial_dataB * drBtheta_radial_dataB + np.divide(Bcc2_radial_dataB * dthetaBtheta_radial_dataB +
                                  Bcc3_radial_dataB * dphiBtheta_radial_dataB +
                                  Bcc2_radial_dataB * Bcc1_radial_dataB, dataB['x1v'])

        total_mag_forceB = magtensionB + pmagforceB

    if total_force:
        dthetapressA = np.gradient(dataA['press'], dataA['x2v'], edge_order=2,
                                   axis=1) / r  # derivative of Bcc1 with respect to theta
        dthetapressB = np.gradient(dataB['press'], dataB['x2v'], edge_order=2,
                                   axis=1) / r  # derivative of Bcc1 with respect to theta

        pressforceA = - 0.5 * (dthetapressA[phi_index, theta_indexA, :] + dthetapressA[phi_index, theta_indexA - 1, :])
        pressforceB = - 0.5 * (dthetapressB[phi_index, theta_indexB, :] + dthetapressB[phi_index, theta_indexB - 1, :])

        dthetapmagA = np.gradient(pmagdataA, dataA['x2v'], edge_order=2,
                                  axis=1)  # derivative of Bcc1 with respect to theta
        dthetapmagB = np.gradient(pmagdataB, dataB['x2v'], edge_order=2,
                                  axis=1)  # derivative of Bcc1 with respect to theta

        Bcc1_radial_dataA = dataA['Bcc1'][phi_index, theta_indexA, :]
        Bcc1_radial_dataB = dataB['Bcc1'][phi_index, theta_indexB, :]
        Bcc2_radial_dataA = dataA['Bcc2'][phi_index, theta_indexA, :]
        Bcc2_radial_dataB = dataB['Bcc2'][phi_index, theta_indexB, :]
        Bcc3_radial_dataA = dataA['Bcc3'][phi_index, theta_indexA, :]
        Bcc3_radial_dataB = dataB['Bcc3'][phi_index, theta_indexB, :]

        pmagforceA = dthetapmagA[phi_index, theta_indexA, :]
        pmagforceB = dthetapmagB[phi_index, theta_indexB, :]

        dthetaBthetaA = np.gradient(dataA['Bcc2'], dataA['x2v'], edge_order=2,
                                    axis=1)  # derivative of Bcc1 with respect to theta
        dthetaBthetaB = np.gradient(dataB['Bcc2'], dataB['x2v'], edge_order=2,
                                    axis=1)  # derivative of Bcc1 with respect to theta
        dphiBthetaA = np.gradient(dataA['Bcc2'], dataA['x3v'], edge_order=2,
                                  axis=0)  # derivative of Bcc1 with respect to phi
        dphiBthetaB = np.gradient(dataB['Bcc2'], dataB['x3v'], edge_order=2,
                                  axis=0)  # derivative of Bcc1 with respect to phi
        drBthetaA = np.gradient(dataA['Bcc2'], dataA['x1v'], edge_order=2, axis=2)
        drBthetaB = np.gradient(dataB['Bcc2'], dataB['x1v'], edge_order=2, axis=2)

        dthetaBtheta_radial_dataA = dthetaBthetaA[phi_index, theta_indexA, :]
        dthetaBtheta_radial_dataB = dthetaBthetaB[phi_index, theta_indexB, :]
        dphiBtheta_radial_dataA = dphiBthetaA[phi_index, theta_indexA, :]
        dphiBtheta_radial_dataB = dphiBthetaB[phi_index, theta_indexB, :]
        drBtheta_radial_dataA = drBthetaA[phi_index, theta_indexA, :]
        drBtheta_radial_dataB = drBthetaB[phi_index, theta_indexA, :]

        # big ass equation (magnetic tension)
        magtensionA = Bcc1_radial_dataA * drBtheta_radial_dataA + np.divide(
            Bcc2_radial_dataA * dthetaBtheta_radial_dataA +
            Bcc3_radial_dataA * dphiBtheta_radial_dataA +
            Bcc2_radial_dataA * Bcc1_radial_dataA, dataA['x1v'])

        total_mag_forceA = magtensionA + pmagforceA
        total_forceA = total_mag_forceA + pressforceA

        magtensionB = Bcc1_radial_dataB * drBtheta_radial_dataB + np.divide(
            Bcc2_radial_dataB * dthetaBtheta_radial_dataB +
            Bcc3_radial_dataB * dphiBtheta_radial_dataB +
            Bcc2_radial_dataB * Bcc1_radial_dataB, dataB['x1v'])

        total_mag_forceB = magtensionB + pmagforceB
        total_forceB = total_mag_forceB + pressforceB


    #plot constant Beta
    plt.figure()
    if mag_force:
        plt.plot(dataA['x1v'], pmagforceA, label="Magnetic Pressure Force", linestyle="--", color='red')
        plt.plot(dataA['x1v'], magtensionA, label="Magnetic Tension Force", linestyle="--", color='purple')
        plt.plot(dataA['x1v'], total_mag_forceA, label='Total Magnetic Force', linestyle=":", color='black')
        if lines:
            plt.axvline(5, -100, 100, label='Inner Radius', linestyle='dotted', color='green')
            plt.axvline(10, -100, 100, label='Highest Density', linestyle='dotted', color='green')
    if nonmag_force:
        plt.plot(dataA['x1v'], pressforceA, label="Gas Pressure Force", linestyle="--", color='blue')
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
    plt.ylabel("Vertical Forces")
    plt.title("Time: {}".format(simulation_timeA) + "\nConstant Beta")
    plt.tight_layout()
    plt.legend()

    if mag_force:
        figname = "forces_at_timestep_{}".format(timestep)
        if specific:
            filedir ="C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/vertical_forces_rprofile_local/specific/Beta/mag/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            #print(filedir)
            #print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        else:
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/vertical_forces_rprofile_local/Beta/mag/"
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
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/vertical_forces_rprofile_local/specific/Beta/gas/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        else:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/vertical_forces_rprofile_local/Beta/gas/"
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
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/vertical_forces_rprofile_local/specific/Beta/total/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        else:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/vertical_forces_rprofile_local/Beta/total/"
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
    plt.ylabel("Vertical Forces")
    plt.title("Time: {}".format(simulation_timeA) + "\nConstant B")
    plt.tight_layout()
    plt.legend()

    if mag_force:
        if specific:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/vertical_forces_rprofile_local/specific/B/mag/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        else:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/vertical_forces_rprofile_local/B/mag/"
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
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/vertical_forces_rprofile_local/specific/B/gas/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        else:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/vertical_forces_rprofile_local/B/gas/"
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
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/vertical_forces_rprofile_local/specific/B/total/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        else:
            figname = "forces_at_timestep_{}".format(timestep)
            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/direction_quantity_profiles/vertical_forces_rprofile_local/B/total/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            # print(filedir)
            # print("Saving figure " + filedir + figname)
            plt.savefig(filedir + figname)
        #plt.show()
        plt.gca().clear()
        plt.close()