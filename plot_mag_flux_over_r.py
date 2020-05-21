import numpy as np
import sys
sys.path.append('C:/Users/Ellie/Downloads/nerd/scripts/GRvis-master/modules/')
from raw_data_utils import *
from metric_utils import *
from calculate_data_utils import *
import matplotlib.pyplot as plt
from setup_manager import *
import argparse

def plot_mag_flux_over_r(setup_manager, times, **kwargs):
    """
    Currently, a very bare-bones method to calculate and plot the
    magnetic flux over spherical shells at given times over radius.

    INPUTS:
    - setup_manager: contains config, specs, project, and path information.
    - times: the timesteps to be loaded, i.e. an array of integers.

    TO DO:
    - add option to bypass config/specs
    - add option to time average
    - add option to load from file instead of calculating
    - make legend over time nicer
    """
    load_from_file = kwargs.get("load_from_file", True)
    write_to_file = kwargs.get("write_to_file", True)

    metric = None
    coords = None

    raw_data_path = setup_manager.path_to_raw_data
    reduced_data_path = setup_manager.path_to_reduced_data
    figure_path = setup_manager.path_to_figures
    if not os.path.exists(reduced_data_path):
        os.makedirs(reduced_data_path)
    if not os.path.exists(figure_path):
        os.makedirs(figure_path)

    figure_name = figure_path + "mag_flux_over_r_t"

    for time in times:
        figure_name += "-{:05d}".format(time)
        # ----- First have to get data ------
        filename = raw_data_path + config + "_" + specs + ".prim.{:05}.athdf".format(time)
        raw_data = read_athdf(filename, quantities=["Bcc1"])
        sim_time = raw_data["Time"]
        x1v = raw_data["x1v"]
        x2v = raw_data["x2v"]
        x3v = raw_data["x3v"]

        # don't need to load metric every time since the coords stay the same
        if metric is None:
            metric = kerrschild(x1v, x2v, x3v)
        if coords is None:
            coords = (x1v, x2v, x3v)

        mag_flux_over_r = calculate_magnetic_flux_over_radial_shells((raw_data["Bcc1"], coords))

        reduced_data_fp = reduced_data_path + "mag_flux_over_r_at_t{:05}.txt".format(time)

        if write_to_file:
            x1vstr = ', '.join(map(str, raw_data["x1v"]))
            np.savetxt(reduced_data_fp, mag_flux_over_r, header=x1vstr)


        label = "$t={:.2f}~GM/c^3$".format(sim_time)

        # ----- Now the actual plotting ------
        # There is no minus sign here according to WSQ 2019 Eq. 5
        plt.plot(raw_data["x1v"], mag_flux_over_r, marker="o", markersize=3, label=label)

    plt.xlabel(r"$r/r_g$")
    plt.ylabel(r"$\Phi$")
    titstr = config
    titstr += "\n Magnetic flux through each radial shell"
    titstr += "\n $t=${:.2f}$~GM/c^3$".format(sim_time)
    plt.legend()
    plt.title(titstr)
    plt.tight_layout()

    figure_name = kwargs.get("figure_name", figure_name)
    if figure_name == "show":
        plt.show()
    else:
        plt.savefig(figure_name + ".png", bbox_inches='tight')



if __name__ == '__main__':
    config = "1.1.1-torus2_b-gz2"
    specs = "a0beta500torBeta_br32x32x64rl2x2"
    project = "ellie"
    configuration = "lia-hp"

    times = [26, 670]
    kwargs = {}
    # kwargs["figure_name"] = "show"

    sim_name = config + "_" + specs
    sm = setup_manager(configuration, project, sim_name)

    plot_mag_flux_over_r(sm, times, **kwargs)

