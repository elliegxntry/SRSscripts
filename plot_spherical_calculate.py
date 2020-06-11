# import packages
import sys
sys.path.append("C:/Users/Ellie/Downloads/nerd/athena-public-version/vis/python/modules")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import kerrmetric as kerr
import numpy as np
import os
import new_athena_read

# specifications
times = np.arange(0, 671)
dist = "Beta"
quantity = "vel1"
do_vertical = False
do_diff = True
do_average = False
color_log = False
r_max = 15
xlims = [-r_max, r_max]
ylims = [-r_max, r_max]
cmap = "RdBu"
vmin = None
vmax = None

raw_quantities = ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
calc_quantities = ['u1', 'rhou1', 'bsq', 'beta', 'test']

# path to load data
config = "1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2"
data_file = "C:/Users/Ellie/Downloads/nerd/SRSData/"
sim_str = "Constant_" + dist + "_"
datapath_base = data_file + config

def get_calculated_data(data_in, calc_q):
    r = data_in["x1v"]
    theta = data_in["x2v"]
    metric = kerr.kerr(r, theta)
    fourvector = kerr.fourvector(r, theta, data_in['vel1'], data_in['vel2'], data_in['vel3'], data_in['Bcc1'], data_in['Bcc2'],
                                 data_in['Bcc3'])
    if calc_q == 'u1':
        data_in[calc_q] = fourvector.u1
    elif calc_q == 'rhou1':
        data_in[calc_q] = data_in['rho'] * fourvector.u1
    elif calc_q == 'bsq':
        data_in[calc_q] = fourvector.bsq
    elif calc_q == 'beta':
        data_in[calc_q] = data_in['press'] / fourvector.pmag
    elif calc_q == 'test':
        data_in[calc_q] = 2.0 * data_in['press'] / (data_in['Bcc3'] ** 2)
    return data_in

if do_diff:
    datat0_name = data_file + config + "/" + config + ".prim.00000.athdf"
    if quantity in raw_quantities:
        datat0 = athena_read.athdf(datat0_name, quantities=[quantity])
    elif quantity in calc_quantities:
        datat0 = athena_read.athdf(datat0_name, quantities=raw_quantities)
        datat0 = get_calculated_data(datat0, quantity)

for timestep in times:
    timestep = "{:05d}".format(int(timestep))
    filepath = datapath_base + "/" + config + ".prim." + timestep + ".athdf"
    #print("Loading time step {}".format(timestep))

    # Main function
    filename = sim_str + quantity
    if do_diff:
        filename += "_diff"
    filename += "_" + timestep

    def theta_func(xmin, xmax, _, nf):
        x2_vals = np.linspace(xmin, xmax, nf)
        theta_vals = x2_vals + (1.0 - h) / 2.0 * np.sin(2.0 * x2_vals)
        return theta_vals

    if quantity in raw_quantities:
        # Read data
        data = athena_read.athdf(filepath, quantities=[quantity])
    elif quantity in calc_quantities:
        data = athena_read.athdf(filepath, quantities=raw_quantities)
        data = get_calculated_data(data, quantity)
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
    if not do_vertical:
        r_grid, phi_grid = np.meshgrid(r_face, phi_face)
        x_grid = r_grid * np.cos(phi_grid)
        y_grid = r_grid * np.sin(phi_grid)
    else:
        theta_face_extended = np.concatenate((theta_face, 2.0*np.pi - theta_face[-2::-1]))
        r_grid, theta_grid = np.meshgrid(r_face, theta_face_extended)
        x_grid = r_grid * np.sin(theta_grid)
        y_grid = r_grid * np.cos(theta_grid)

    # Perform slicing/averaging of scalar data
    if not do_vertical:
        if nx2 % 2 == 0:
            if do_diff:
                vals = np.mean(data[quantity][:, int(nx2 / 2) - 1:int(nx2 / 2) + 1, :]
                               - datat0[quantity][:, int(nx2/2)-1:int(nx2/2) + 1, :], axis=1)
            else:
                vals = np.mean(data[quantity][:, int(nx2/2)-1:int(nx2/2)+1, :], axis=1)
        else:
            if do_diff:
                vals = data[quantity][:, int(nx2/2), :] - datat0[quantity][:, int(nx2/2), :]
            else:
                vals = data[quantity][:, int(nx2/2), :]
        if do_average:
            vals = np.repeat(np.mean(vals, axis=0, keepdims=True), nx3, axis=0)
    else:
        if do_average:
            if do_diff:
                vals_right = np.mean(data[quantity] - datat0[quantity], axis=0)
                vals_left = vals_right
            else:
                vals_right = np.mean(data[quantity], axis=0)
                vals_left = vals_right
        else:
            vals_right = 0.5 * (data[quantity]
                                [-1, :, :] + data[quantity][0, :, :])
            vals_left = 0.5 * (data[quantity][int((nx3/2))-1, :, :]
                               + data[quantity][int(nx3 / 2), :, :])
            if do_diff:
                vals_right = vals_right - 0.5 * (datat0[quantity]
                                [-1, :, :] + datat0[quantity][0, :, :])
                vals_left = vals_left - 0.5 * (datat0[quantity][int((nx3/2))-1, :, :]
                               + datat0[quantity][int(nx3 / 2), :, :])

        # Join scalar data through boundaries
        vals = np.vstack((vals_right, vals_left[::-1, :]))

    # Determine colormapping properties
    if color_log:
        norm = colors.LogNorm()
    else:
        norm = colors.Normalize()

    # Make plot
    plt.figure()
    title = "Constant " + dist + ": " + quantity
    if do_diff:
        title += "-" + quantity + "(t=0)"
    title += " at timestep " + timestep
    plt.title(title)
    im = plt.pcolormesh(x_grid, y_grid, vals, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)

    plt.gca().set_aspect('equal')
    plt.xlim(xlims)
    plt.ylim(ylims)

    if not do_vertical:
        plt.xlabel(r'$x / r$')
        plt.ylabel(r'$y / r$')
    else:
        plt.xlabel(r'$x / r$')
        plt.ylabel(r'$z / r$')

    plt.colorbar(im)

    filedir = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Slices/" + quantity
    if do_diff:
        filedir += "Diff"
    if do_average:
        filedir += "Average"
    if do_vertical:
        filedir += "Vertical"
    else:
        filedir += "Midplane"
    filedir += "Slices/" + dist + "/"
    if not os.path.isdir(filedir):
        os.makedirs(filedir)
    #print(filedir)
    plt.savefig(filedir + filename, bbox_inches='tight')
    plt.show()
    plt.close()