#! /usr/bin/env python

"""
Script for plotting vertical (r,theta) or midplane (r,phi) slices of data in
spherical coordinates.

Run "plot_spherical.py -h" to see description of inputs.

See documentation on athena_read.athdf() for important notes about reading files
with mesh refinement.

The -c <colormap> option is highly recommended to change the default. Consider
"RdBu_r" or "gist_heat" for example.

Users are encouraged to make their own versions of this script for improved
results by adjusting figure size, spacings, tick locations, axes labels, etc.
The script must also be modified to plot any functions of the quantities in the
file, including combinations of multiple quantities.

"""

# Python standard modules
import argparse
import warnings

# Other Python modules
import h5py
import numpy as np

# Athena++ modules
import athena_read


# Main function
def main(**kwargs):
    # Load Python plotting modules
    if kwargs['output_file'] != 'show':
        import matplotlib
        matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    # Determine refinement level to use
    if kwargs['level'] is not None:
        level = kwargs['level']
    else:
        level = None

    # Determine if vector quantities should be read
    quantities = [kwargs['quantity']]

    # Define grid compression in theta-direction
    h = kwargs['theta_compression'] if kwargs['theta_compression'] is not None else 1.0

    def theta_func(xmin, xmax, _, nf):
        x2_vals = np.linspace(xmin, xmax, nf)
        theta_vals = x2_vals + (1.0 - h) / 2.0 * np.sin(2.0 * x2_vals)
        return theta_vals

    # Read data
    data = athena_read.athdf(kwargs['data_file'], quantities=quantities, level=level,
                             face_func_2=theta_func)

    # Extract basic coordinate information
    with h5py.File(kwargs['data_file'], 'r') as f:
        coordinates = f.attrs['Coordinates']
    r = data['x1v']
    theta = data['x2v']
    phi = data['x3v']
    nx1 = len(r)
    nx2 = len(theta)
    nx3 = len(phi)

    print(nx1)
    print(nx2)
    print(nx3)

    # Set radial extent
    if kwargs['r_max'] is not None:
        r_max = kwargs['r_max']
    else:
        r_max = data['x1f'][-1]

    # Account for logarithmic radial coordinate
    if kwargs['logr']:
        r = np.log10(r)
        r_max = np.log10(r_max)

    # Create scalar grid
    if kwargs['midplane']:
        # phi doesn't go all the way to 2pi, so extend it
        phi_extended = \
            np.concatenate((phi[-1:] - 2.0 * np.pi, phi, phi[:1] + 2.0 * np.pi))
        phi_extended -= 0.5 * 2.0 * np.pi / nx3
        r_grid, phi_grid = np.meshgrid(r, phi_extended)
        x_grid = r_grid * np.cos(phi_grid)
        y_grid = r_grid * np.sin(phi_grid)
    else:
        theta_extended = np.concatenate((-theta[0:1], theta, 2.0 * np.pi - theta[::-1],
                                         2.0 * np.pi + theta[0:1]))
        theta_extended_corrected = theta_extended - 0.5 * (theta[1] - theta[0])
        r_grid, theta_grid = np.meshgrid(r, theta_extended_corrected)
        x_grid = r_grid * np.sin(theta_grid)
        y_grid = r_grid * np.cos(theta_grid)


    # Perform slicing/averaging of scalar data
    if kwargs['midplane']:
        if nx2 % 2 == 0:
            vals = np.mean(data[kwargs['quantity']][:, int(nx2 / 2 - 1):int(nx2 / 2 + 1), :], axis=1)
        else:
            vals = data[kwargs['quantity']][:, int(nx2 / 2), :]
        if kwargs['average']:
            vals = np.repeat(np.mean(vals, axis=0, keepdims=True), nx3, axis=0)
    else:
        if kwargs['average']:
            vals_right = np.mean(data[kwargs['quantity']], axis=0)
            vals_left = vals_right
        else:
            vals_right = 0.5 * (data[kwargs['quantity']]
                                [-1, :, :] + data[kwargs['quantity']][0, :, :])
            vals_left = 0.5 * (data[kwargs['quantity']][int((nx3/2)-1), :, :]
                               + data[kwargs['quantity']][int(nx3 / 2), :, :])

    # Join scalar data through boundaries
    if kwargs['midplane']:
        vals = np.vstack((vals[-1:, :], vals, vals[:1, :]))
    else:
        vals = np.vstack((vals_left[:1, :], vals_right,
                          vals_left[::-1, :], vals_right[:1, :]))


    # Determine colormapping properties
    cmap = plt.get_cmap(kwargs['colormap'])
    vmin = kwargs['vmin']
    vmax = kwargs['vmax']
    if kwargs['logc']:
        norm = colors.LogNorm()
    else:
        norm = colors.Normalize()

    # Make plot
    plt.figure()
    im = plt.pcolormesh(x_grid, y_grid, vals, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
    plt.gca().set_aspect('equal')
    plt.xlim((-r_max, r_max))
    plt.ylim((-r_max, r_max))
    r_string = r'\log_{10}(r)\ ' if kwargs['logr'] else r'r\ '
    angle_string_x = r'\cos(\phi)' if kwargs['midplane'] else r'\sin(\theta)'
    angle_string_y = r'\sin(\phi)' if kwargs['midplane'] else r'\cos(\theta)'
    plt.xlabel('$' + r_string + angle_string_x + '$')
    plt.ylabel('$' + r_string + angle_string_y + '$')
    plt.colorbar(im)
    if kwargs['output_file'] == 'show':
        plt.show()
    else:
        plt.savefig(kwargs['output_file'], bbox_inches='tight')


# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('data_file',
                        help='name of input file, possibly including path')
    parser.add_argument('quantity',
                        help='name of quantity to be plotted')
    parser.add_argument('output_file',
                        help=('name of output to be (over)written, possibly including '
                              'path; use "show" to show interactive plot instead'))
    parser.add_argument('-m',
                        '--midplane',
                        action='store_true',
                        help=('flag indicating plot should be midplane (r,phi) rather '
                              'than (r,theta)'))
    parser.add_argument('-a', '--average',
                        action='store_true',
                        help='flag indicating phi-averaging should be done')
    parser.add_argument('-l',
                        '--level',
                        type=int,
                        default=None,
                        help=('refinement level to be used in plotting (default: max '
                              'level in file)'))
    parser.add_argument('-r', '--r_max',
                        type=float,
                        default=None,
                        help='maximum radial extent of plot')
    parser.add_argument('--logr',
                        action='store_true',
                        help='flag indicating data should be logarithmically in radius')
    parser.add_argument('-c',
                        '--colormap',
                        default=None,
                        help=('name of Matplotlib colormap to use instead of default'))
    parser.add_argument('--vmin',
                        type=float,
                        default=None,
                        help=('data value to correspond to colormap minimum; use '
                              '--vmin=<val> if <val> has negative sign'))
    parser.add_argument('--vmax',
                        type=float,
                        default=None,
                        help=('data value to correspond to colormap maximum; use '
                              '--vmax=<val> if <val> has negative sign'))
    parser.add_argument('--logc',
                        action='store_true',
                        help='flag indicating data should be colormapped logarithmically')
    parser.add_argument('--theta_compression',
                        type=float,
                        default=None,
                        help=('compression parameter h in '
                              'theta = pi*x_2 + (1-h)/2 * sin(2*pi*x_2)'))
    args = parser.parse_args()
    main(**vars(args))
