# import packages
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import numpy as np
from modules import new_athena_read

times = np.arange(0, 1)
dist = "Beta"
do_average = False
color_log = False
r_max = 15
xlims = [-r_max, r_max]
ylims = [-r_max, r_max]
cmap = "viridis"
vmin = None
vmax = None
mag_press = False
mag_tension = False
gas_press = True

# load data
sim_str = "Constant_" + dist + "_"
config = "1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2"
data_file = "C:/Users/Ellie/Downloads/nerd/SRSData/"
datapath_base = data_file + config

for timestep in times:
    timestep = "{:05d}".format(int(timestep))
    filepath = datapath_base + "/" + config + ".prim." + timestep + ".athdf"
    print("Loading time step {}".format(timestep))
    if mag_press:
        filename = sim_str + "_mag_press" + timestep
    if gas_press:
        filename = sim_str + "_gas_press" + timestep
    if mag_tension:
        filename = sim_str + "_mag_tension" + timestep

    # find theta values
    def theta_func(xmin, xmax, _, nf):
        x2_vals = np.linspace(xmin, xmax, nf)
        theta_vals = x2_vals + (1.0 - h) / 2.0 * np.sin(2.0 * x2_vals)
        return theta_vals

    # Read data
    data = new_athena_read.athdf(filepath)

    # Extract basic coordinate information
    coordinates = data['Coordinates']
    time = data['Time']
    r = data['x1v']
    theta = data['x2v']
    phi = data['x3v']
    r_face = data['x1f']
    theta_face = data['x2f']
    phi_face = data['x3f']
    nx1 = len(r)
    nx2 = len(theta)
    nx3 = len(phi)

    # define variables
    pressdata = data['press']
    pmagdata = 0.5 * (data['Bcc1'] ** 2 + data['Bcc2'] ** 2 + data['Bcc3'] ** 2)

    # Create scalar grid
    theta_face_extended = np.concatenate((theta_face, 2.0*np.pi - theta_face[-2::-1]))
    r_grid, theta_grid = np.meshgrid(r_face, theta_face_extended)
    x_grid = r_grid * np.sin(theta_grid)
    y_grid = r_grid * np.cos(theta_grid)

    # perform calculations
    if gas_press:
        drpress = np.gradient(data['press'], data['x1v'], edge_order=2, axis=2)

        pressforce = - 0.5 * (drpress[phi, theta, :] + drpress[phi, theta - 1, :])

        # Perform slicing/averaging of scalar data
        if do_average:
            vals_right = np.mean(pressforce)
            vals_left = vals_right
        else:
            vals_right = 0.5 * (pressforce
                                [-1, :, :] + pressforce[0, :, :])
            vals_left = 0.5 * (pressforce[int((nx3 / 2)) - 1, :, :]
                               + pressforce[int(nx3 / 2), :, :])

            # Join scalar data through boundaries
            vals = np.vstack((vals_right, vals_left[::-1, :]))

        # Determine colormapping properties
        if color_log:
            norm = colors.LogNorm()
        else:
            norm = colors.Normalize()

        # plot figure
        plt.figure()
        im = plt.pcolormesh(x_grid, y_grid, vals, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
        plt.gca().set_aspect('equal')
        plt.xlim(xlims)
        plt.ylim(ylims)

        if not do_average:
            plt.xlabel(r'$x / r$')
            plt.ylabel(r'$y / r$')
        else:
            plt.xlabel(r'$x / r$')
            plt.ylabel(r'$z / r$')

        plt.colorbar(im)

        filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Slices/forcesSlices/gas_press/" + dist + "/"
        if not os.path.isdir(filedir):
          os.makedirs(filedir)
        print(filedir)
        plt.savefig(filedir + filename, bbox_inches='tight')
        plt.show()
        plt.close()

        if mag_press:
            # mag forces
            pmag_radial_data = 0.5 * (pmagdata[phi, theta, :] + pmagdata[phi, theta - 1, :])
            Bcc1_radial_data = data['Bcc1'][phi, theta, :]
            Bcc2_radial_data = data['Bcc2'][phi, theta, :]
            Bcc3_radial_data = data['Bcc3'][phi, theta, :]

            pmagforce = -np.gradient(pmag_radial_data, data['x1v'], edge_order=2)

            if do_average:
                vals_right = np.mean(pmagforce)
                vals_left = vals_right
            else:
                vals_right = 0.5 * (pmagforce
                                    [-1, :, :] + pmagforce[0, :, :])
                vals_left = 0.5 * (pmagforce[int((nx3 / 2)) - 1, :, :]
                                   + pmagforce[int(nx3 / 2), :, :])

                # Join scalar data through boundaries
                vals = np.vstack((vals_right, vals_left[::-1, :]))

            # Determine colormapping properties
            if color_log:
                norm = colors.LogNorm()
            else:
                norm = colors.Normalize()

            # plot figure
            plt.figure()
            im = plt.pcolormesh(x_grid, y_grid, vals, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
            plt.gca().set_aspect('equal')
            plt.xlim(xlims)
            plt.ylim(ylims)

            if not do_average:
                plt.xlabel(r'$x / r$')
                plt.ylabel(r'$y / r$')
            else:
                plt.xlabel(r'$x / r$')
                plt.ylabel(r'$z / r$')

            plt.colorbar(im)

            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Slices/forcesSlices/mag_press/" + dist + "/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            print(filedir)
            plt.savefig(filedir + filename, bbox_inches='tight')
            plt.show()
            plt.close()

        if mag_tension:
            Bcc1_radial_data = data['Bcc1'][phi, theta, :]
            Bcc2_radial_data = data['Bcc2'][phi, theta, :]
            Bcc3_radial_data = data['Bcc3'][phi, theta, :]

            drBr = np.gradient(Bcc1_radial_data, data['x1v'], edge_order=2)  # derivative of Bcc1 with respect to radius
            dthetaBr = np.gradient(data['Bcc1'], data['x2v'], edge_order=2,
                                    axis=1)  # derivative of Bcc1 with respect to theta
            dphiBr = np.gradient(data['Bcc1'], data['x3v'], edge_order=2,
                                  axis=0)  # derivative of Bcc1 with respect to phi

            drBr_radial_data = drBr
            dthetaBr_radial_data = dthetaBr[phi, theta, :]
            dphiBr_radial_data = dphiBr[phi, theta, :]

            # big ass equation (magnetic tension)
            magtension = Bcc1_radial_data * drBr_radial_data \
                          + np.divide(Bcc2_radial_data * dthetaBr_radial_data +
                                      Bcc3_radial_data * dphiBr_radial_data -
                                      Bcc2_radial_data ** 2 - Bcc3_radial_data ** 2, data['x1v'])

            if do_average:
                vals_right = np.mean(magtension)
                vals_left = vals_right
            else:
                vals_right = 0.5 * (magtension
                                    [-1, :, :] + magtension[0, :, :])
                vals_left = 0.5 * (magtension[int((nx3 / 2)) - 1, :, :]
                                   + magtension[int(nx3 / 2), :, :])

                # Join scalar data through boundaries
                vals = np.vstack((vals_right, vals_left[::-1, :]))

            # Determine colormapping properties
            if color_log:
                norm = colors.LogNorm()
            else:
                norm = colors.Normalize()

            # plot figure
            plt.figure()
            im = plt.pcolormesh(x_grid, y_grid, vals, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
            plt.gca().set_aspect('equal')
            plt.xlim(xlims)
            plt.ylim(ylims)

            if not do_average:
                plt.xlabel(r'$x / r$')
                plt.ylabel(r'$y / r$')
            else:
                plt.xlabel(r'$x / r$')
                plt.ylabel(r'$z / r$')

            plt.colorbar(im)

            filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Slices/forcesSlices/mag_tension/" + dist + "/"
            if not os.path.isdir(filedir):
                os.makedirs(filedir)
            print(filedir)
            plt.savefig(filedir + filename, bbox_inches='tight')
            plt.show()
            plt.close()