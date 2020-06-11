# import packages
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import numpy as np
from scripts import athena_read

times = np.arange(0, 1)
dist = "Beta"
quantity = "rho"
do_midplane = True
do_average = False
color_log = False
r_max = 15
xlims = [-r_max, r_max]
ylims = [-r_max, r_max]
cmap = "viridis"
vmin = 1e0
vmax = -1e-3

# load data
sim_str = "Constant_" + dist + "_"
config = "1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2"
data_file = "C:/Users/Ellie/Downloads/nerd/SRSData/"
datapath_base = data_file + config

for timestep in times:
    timestep = "{:05d}".format(int(timestep))
    filepath = datapath_base + "/" + config + ".prim." + timestep + ".athdf"
    print("Loading time step {}".format(timestep))

    # Main function
    filename = sim_str + quantity + "_" + timestep

    def theta_func(xmin, xmax, _, nf):
        x2_vals = np.linspace(xmin, xmax, nf)
        theta_vals = x2_vals + (1.0 - h) / 2.0 * np.sin(2.0 * x2_vals)
        return theta_vals

    # Read data
    data = athena_read.athdf(filepath, quantities=[quantity])

    # Extract basic coordinate information
    coordinates = data['Coordinates']
    r = data['x1v']
    theta = data['x2v']
    phi = data['x3v']
    r_face = data['x1f']
    theta_face = data['x2f']
    phi_face = data['x3f']
    nx1 = len(r)
    nx2 = len(theta)
    nx3 = len(phi)

    # Create scalar grid
    if not do_midplane:
        r_grid, phi_grid = np.meshgrid(r_face, phi_face)
        x_grid = r_grid * np.cos(phi_grid)
        y_grid = r_grid * np.sin(phi_grid)
    else:
        theta_face_extended = np.concatenate((theta_face, 2.0*np.pi - theta_face[-2::-1]))
        r_grid, theta_grid = np.meshgrid(r_face, theta_face_extended)
        x_grid = r_grid * np.sin(theta_grid)
        y_grid = r_grid * np.cos(theta_grid)

    # Perform slicing/averaging of scalar data
    if not do_midplane:
        if nx2 % 2 == 0:
            vals = np.mean(data[quantity][:, int(nx2/2)-1:int(nx2/2)+1, :], axis=1)
        else:
            vals = data[quantity][:, int(nx2/2), :]
        if do_average:
            vals = np.repeat(np.mean(vals, axis=0, keepdims=True), nx3, axis=0)
    else:
        if do_average:
            vals_right = np.mean(data[quantity], axis=0)
            vals_left = vals_right
        else:
            vals_right = 0.5 * (data[quantity]
                                [-1, :, :] + data[quantity][0, :, :])
            vals_left = 0.5 * (data[quantity][int((nx3/2))-1, :, :]
                               + data[quantity][int(nx3 / 2), :, :])

        # Join scalar data through boundaries
        vals = np.vstack((vals_right, vals_left[::-1, :]))

    # Determine colormapping properties

    if color_log:
        norm = colors.LogNorm()
    else:
        norm = colors.Normalize()

    # Make plot
    plt.figure()
    plt.title("Constant " + dist + ": " + quantity + " at timestep " + timestep)
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

    filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Slices/" + quantity + "Slices/B/Sensitive/"
    if not os.path.isdir(filedir):
      os.makedirs(filedir)
    #print(filedir)
    #plt.savefig(filedir + filename, bbox_inches='tight')
    plt.show()
    plt.close()