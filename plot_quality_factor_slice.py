import numpy as np
import sys
sys.path.append('../modules/')
from raw_data_utils import *
from metric_utils import *
from calculate_data_utils import *
from plotting_utils import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import pickle

def plot_quality_factor_slice(config, specs, raw_data_path, reduced_data_path, time, **kwargs):
    """
    A script to calculate and plot a slice of the
    toroidal quality factor at a given time.
    Plotting capabilities are taken from Athena++'s plot_spherical.py.
    Equation for Qphi taken from White, Quataert, and Blaes 2019 Eq. 1.

    INPUTS:
    - config: the configuration of the Athena++ executable (LH naming convention)
    - specs: the specifications of the desired simulation (LH naming convention)
    - raw_data_path: the directory where the raw data is housed
        ...it does not include the file names.
        .athdf filename will be constructed from config and specs.
    - reduced_data_path: the directory where the reduced data should be saved,
        which is needed both to write to the file (if calculating from scratch)
        and to load from (if not calculating from scratch)
    - time: the timestep to be loaded, i.e. an integer.

    TO DO:
    - add option to bypass config/specs
    - add option to time average
    - add option to load from file instead of calculating
    - make legend over time nicer
    - build save figure path

    """

    load_from_file = kwargs.get("load_from_file", True)
    write_to_file = kwargs.get("write_to_file", True)

    metric = None
    coords = None

    if not os.path.exists(reduced_data_path):
        os.makedirs(reduced_data_path)

    # ----- First have to get data ------
    filename = raw_data_path + config + "_" + specs + ".prim.{:05}.athdf".format(time)
    raw_data = read_athdf(filename, quantities=["rho", "vel1", "vel2", "vel3", "press", "Bcc1", "Bcc2", "Bcc3"])
    out_vel = (raw_data["vel1"], raw_data["vel2"], raw_data["vel3"])
    out_B = (raw_data["Bcc1"], raw_data["Bcc2"], raw_data["Bcc3"])
    press = raw_data["press"]
    density = raw_data["rho"]

    sim_time = raw_data["Time"]
    x1v = raw_data["x1v"]
    x2v = raw_data["x2v"]; nx2 = x2v.size
    x3v = raw_data["x3v"]; nx3 = x3v.size

    metric = kerrschild(x1v, x2v, x3v)
    coords = (x1v, x2v, x3v)

    four_velocity = metric.get_four_velocity_from_output(out_vel)
    proj_b = metric.get_proj_bfield_from_outputB_fourV(out_B, four_velocity)

    quality_factor_phi = calculate_quality_factor_phi((proj_b, press, density, coords))
    data = {"quantity":quality_factor_phi, "x1v":x1v, "x2v":x2v, "x3v":x3v, "x1f":raw_data["x1f"], "x2f":raw_data["x2f"], "x3f":raw_data["x3f"]}

    if write_to_file:
        reduced_data_fp = reduced_data_path + "quality_factor_phi_at_t{:05}.txt".format(time)
        pickle.dump(quality_factor_phi, open(reduced_data_fp, "wb"))

    kwargs['title_str'] = config + "\n Toroidal quality factor slice\n $t=${:.2f}$~GM/c^3$".format(sim_time)
    kwargs['cbar_label'] = r"$Q_\phi$"

    # plt.figure()
    # plt.plot(x1v, quality_factor_phi[0, int(nx2/2), :])
    # plt.title("Toroidal quality factor radial profile at midplane\n $t=${:.2f}$~GM/c^3$".format(sim_time))
    # plt.xlabel("$r$")
    # plt.ylabel(r"$Q_\phi$")
#
    # plt.show()
    plot_slice(data, **kwargs)
    plt.savefig("/work/06165/ahankla/stampede2/torus_vert-flux/figures/qf.png")


if __name__ == '__main__':

    # ----- inputs ----------
    # config = "1.1.1-torus2_b-gz2"
    # specs = "a0beta500torB_br32x32x64rl2x2"
    # times = [801] # np.arange(0, 801)
    # specs = "a0beta500torBeta_br32x32x64rl2x2"
    # times = [26, 801] # np.arange(0, 801)
    config = "vert-flux"
    specs = "initial"
    times = [0, 50, 100]

    # ------ kwargs --------
    kwargs = {}
    kwargs['output_file'] = 'show'
    kwargs["midplane"] = False
    kwargs['r_max'] = None
    kwargs['logr'] = False
    kwargs['average'] = False
    kwargs['colormap'] = 'viridis'
    kwargs['vmin'] = None
    kwargs['vmax'] = None
    kwargs['logc'] = False
    kwargs['stream_average'] = False
    kwargs["stream"] = None

    # build paths
    # raw_data_path = "/mnt/c/Users/liaha/research/projects/ellie/raw_data/" + specs + "/"
    # reduced_data_path = "/mnt/c/Users/liaha/research/projects/ellie/reduced_data/" + specs + "/"
    # XX build savefig path.
    raw_data_path = "/scratch/06165/ahankla/torus_vert-flux/data_raw/vert-flux_initial/"
    reduced_data_path = "/work/06165/ahankla/stampede2/torus_vert-flux/data_reduced/vert-flux_initial/"

    for time in times:
        plot_quality_factor_slice(config, specs, raw_data_path, reduced_data_path, time, **kwargs)
