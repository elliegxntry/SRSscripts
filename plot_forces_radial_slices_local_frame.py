"""
does both midplane and vertical
"""

# import packages
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import numpy as np
from modules import new_athena_read

times = np.arange(0, 1)
dist = "Beta"
do_midplane = True
do_average = False
color_log = False
r_max = 20
if do_midplane:
    xlims = [-r_max, r_max]
else:
    xlims = [0, r_max]
ylims = [-r_max, r_max]
phi_ind = 1
theta_ind = 63
cmap = "viridis"
vmin = None
vmax = None
gas_press = False
mag_press = False
mag_tension = True


# load data
sim_str = "Constant_" + dist + "_"
config = "1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2"
data_file = "C:/Users/Ellie/Downloads/nerd/SRSData/"
datapath_base = data_file + config
title_string = "Constant " + dist + "\n"

for timestep in times:
    timestep = "{:05d}".format(int(timestep))
    filepath = datapath_base + "/" + config + ".prim." + timestep + ".athdf"
    print("Loading time step {}".format(timestep))
    if gas_press:
        filename = sim_str + "_gas_press" + timestep
    elif mag_press:
        filename = sim_str + "_mag_press" + timestep

    elif mag_tension:
        filename = sim_str + "_mag_tension" + timestep

    if do_midplane:
        # define midplane variable
        orientation = "Midplane"
        # Read data
        data = new_athena_read.athdf(filepath, j_min=theta_ind-1, j_max=theta_ind+2)
        print(data['rho'].shape)

        # Extract basic coordinate information
        coordinates = data['Coordinates']
        time = data['Time']
        r_vals = data['x1v']
        theta_vals = data['x2v']
        phi_vals = data['x3v']
        r_face = data['x1f']
        theta_face = data['x2f']
        phi_face = data['x3f']
        nx1 = len(r_vals)
        nx2 = len(theta_vals)
        nx3 = len(phi_vals)

        title_string += 'Midplane slice at theta={:.2f} (index={:d})\n'.format(theta_vals[1], theta_ind)

        # Create scalar grid
        r_grid, phi_grid = np.meshgrid(r_face, phi_face)
        x_grid = r_grid * np.cos(phi_grid)
        y_grid = r_grid * np.sin(phi_grid)

        # perform midplane calculations
        if gas_press:
            # print(data['press'].shape)
            drpress = np.gradient(np.squeeze(data['press'][:, 1, :]), data['x1v'], edge_order=2, axis=1)
            # print(drpress.shape)
            vals = - drpress
            type = "gas_press_force"

        elif mag_press:
            pmagdata = np.squeeze(0.5 * (data['Bcc1'][:, 1, :] ** 2 + data['Bcc2'][:, 1, :] ** 2 + data['Bcc3'][:, 1, :] ** 2))
            pmagforce = -np.gradient(pmagdata, data['x1v'], edge_order=2, axis=1)

            vals = pmagforce
            type = "mag_press_force"

        elif mag_tension:
            Bcc1_radial_data = data['Bcc1'][:, 1, :]
            Bcc2_radial_data = data['Bcc2'][:, 1, :]
            Bcc3_radial_data = data['Bcc3'][:, 1, :]

            drBr = np.gradient(Bcc1_radial_data, data['x1v'], edge_order=2, axis=1)
            dthetaBr = np.gradient(data['Bcc1'], data['x2v'], edge_order=2, axis=1)[:, 1, :]  # derivative of Bcc1 with respect to theta
            dphiBr = np.gradient(data['Bcc1'][:, 1, :], data['x3v'], edge_order=2, axis=0)  # derivative of Bcc1 with respect to phi

            # big ass equation (magnetic tension)
            magtension = Bcc1_radial_data * drBr \
                         + np.divide(Bcc2_radial_data * dthetaBr \
                                     + np.divide(Bcc3_radial_data * dphiBr, np.sin(data['x2v'][1])) \
                                     - Bcc2_radial_data ** 2 - Bcc3_radial_data ** 2, data['x1v'])

            vals = magtension
            type = "mag_tension"

    else:
        # define vertical variable
        orientation = "Vertical"
        # Read data
        data = new_athena_read.athdf(filepath, k_min=phi_ind-1, k_max=phi_ind+2)
        print(data['rho'].shape)

        # Extract basic coordinate information
        coordinates = data['Coordinates']
        time = data['Time']
        r_vals = data['x1v']
        theta_vals = data['x2v']
        r_face = data['x1f']
        theta_face = data['x2f']
        phi_face = data['x3f']
        nx1 = len(r_vals)
        nx2 = len(theta_vals)

        title_string += 'Vertical slice at phi={:.2f} (index={:d})\n'.format(data['x3v'][1], phi_ind)

        # Create scalar grid
        r_grid, theta_grid = np.meshgrid(r_face, theta_face)
        r_grid_v, theta_grid_v = np.meshgrid(r_vals, theta_vals)
        x_grid = r_grid * np.sin(theta_grid)
        y_grid = r_grid * np.cos(theta_grid)

        # perform vertical calculations
        if gas_press:
            #print(data['press'].shape)
            drpress = np.gradient(np.squeeze(data['press'][1, :, :]), data['x1v'], edge_order=2, axis=1)
            #print(drpress.shape)
            vals = - drpress
            type = "gas_press_force"

        elif mag_press:
            pmagdata = np.squeeze(0.5 * (data['Bcc1'][1, :, :] ** 2 + data['Bcc2'][1, :, :] ** 2 + data['Bcc3'][1, :, :] ** 2))
            pmagforce = -np.gradient(pmagdata, data['x1v'], edge_order=2, axis=1)

            vals = pmagforce
            type = "mag_press_force"

        elif mag_tension:
            Bcc1_radial_data = data['Bcc1'][1, :, :]
            Bcc2_radial_data = data['Bcc2'][1, :, :]
            Bcc3_radial_data = data['Bcc3'][1, :, :]

            drBr = np.gradient(Bcc1_radial_data, data['x1v'], edge_order=2, axis=1)
            #print(Bcc1_radial_data.shape)
            #print(drBr.shape)
            #print(data['x3v'].shape)
            #print(data['Bcc1'].shape)
            dthetaBr = np.gradient(Bcc1_radial_data, data['x2v'], edge_order=2, axis=0)  # derivative of Bcc1 with respect to theta
            dphiBr = np.gradient(data['Bcc1'], data['x3v'], edge_order=1, axis=0)[1, :, :]  # derivative of Bcc1 with respect to phi

            print(Bcc3_radial_data.shape)
            print(dphiBr.shape)
            # big ass equation (magnetic tension)
            magtension = Bcc1_radial_data * drBr \
                          + np.divide(Bcc2_radial_data * dthetaBr
                                      + np.divide(Bcc3_radial_data * dphiBr, np.sin(theta_grid_v))
                                      - Bcc2_radial_data ** 2 - Bcc3_radial_data ** 2, data['x1v'])

            vals = magtension
            type = "mag_tension"

    # Determine colormapping properties
    if color_log:
        norm = colors.LogNorm()
    else:
        norm = colors.Normalize()

    title_string += type + "\n$t={:.2f}~[GM/c^3]$ (timestep=".format(data["Time"]) + timestep + ")"

    # Make plot
    plt.figure()
    plt.title(title_string)
    print(x_grid.shape)
    print(y_grid.shape)
    im = plt.pcolormesh(x_grid, y_grid, vals, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
    plt.gca().set_aspect('equal')
    plt.xlim(xlims)
    plt.ylim(ylims)
    plt.xlabel(r'$x / r$')
    plt.ylabel(r'$y / r$')
    plt.colorbar(im)
    plt.tight_layout()

    filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Slices/" + type + orientation + "/Constant" + dist + "/"
    #print(filedir)
    if not os.path.isdir(filedir):
        os.makedirs(filedir)
    plt.savefig(filedir + filename, bbox_inches='tight')
    plt.show()
    plt.close()