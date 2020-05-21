import sys
sys.path.append('C:/Users/Ellie/Downloads/nerd/scripts/GRvis-master/modules/')
from setup_manager import *
from plotting_utils import *


if __name__ == '__main__':
    # ----- inputs ----------
    configuration = "lia-hp"
    project = "ellie"
    config = "1.1.1-torus2_b-gz2"
    specs = "a0beta500torB_br32x32x64rl2x2"

    times = [801]
    quantities_to_plot =  ["rho", "press", "vel1", "vel2", "vel3", "Bcc1", "Bcc2", "Bcc3"]

    # ------ kwargs --------
    kwargs = {}
    kwargs["midplane"] = False
    kwargs['r_max'] = 25
    kwargs['logr'] = False
    kwargs['average'] = False
    kwargs['colormap'] = 'viridis'
    kwargs['vmin'] = None
    kwargs['vmax'] = None
    kwargs['logc'] = False
    kwargs['stream_average'] = False
    kwargs["stream"] = None
    kwargs["quantities"] = ["Bcc2", "rho"]

    sim_name = config + "_" + specs
    sm = setup_manager(configuration, project, sim_name)
    for time in times:
        print("plotting time {}".format(time))
        for quantity in quantities_to_plot:
            figure_name = sm.path_to_figures + quantity + "_t{:05d}_slice_".format(int(time)) + sm.sim_name + ".png"
            kwargs['output_file'] = figure_name
            plot_output_quantity_slice(sm, time, quantity, **kwargs)
