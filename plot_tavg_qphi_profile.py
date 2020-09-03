# import packages
import numpy as np
import sys
sys.path.append('modules/')
sys.path.append('C:/Users/Ellie/Downloads/nerd/scripts/GRvis/scripts/modules/')
from raw_data_utils import *
from reduce_data_utils import *
from metric_utils import *
from calculate_data_utils import *
from plotting_utils import *
from setup_manager import *
import matplotlib
matplotlib.use('tkAgg')
import matplotlib.pyplot as plt
import pickle

# attempt (and probably fail) to take a line out of the slice as a function
def plot_tavg_qphi_profile(setup, time, **kwargs):
    # i think that this is getting the already calculated qphi slice and coordinates
    quality_factor_phi_slice = retrieve_quality_factor_phi_slice(setup, time, **kwargs)
    coords = retrieve_coords(setup, time, **kwargs)
    time_dict = retrieve_time_dict(setup)
    sim_time = time_dict[time]

    # i can't figure out what this does
    data = coords
    data["quantity"] = quality_factor_phi_slice

    #create figure to mess with data
    plt.figure()
    x1v = coords['x1v']
    nx2 = coords['x2v'].shape
    qphi_profile = quality_factor_phi_slice.reshape([0], [1])
    print(quality_factor_phi_slice.shape)

if __name__ == '__main__':
    # ----- inputs ----------
    configuration = "ellie"
    project = "ellie"
    config = "1.1.1-torus2_b-gz2"
    specs = "a0beta500torB_br32x32x64rl2x2"
    sim_name = config + "_" + specs
    overwrite_data = False

    setup = setup_manager(configuration, project, sim_name)
    times = get_list_of_time_steps(setup)
    #times = times[::-100]
    times = np.arange(700, 1211)

    # ------ kwargs --------
    kwargs = {}
    kwargs["overwrite_data"] = overwrite_data
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
    kwargs['setup_val'] = setup
    kwargs['quantity_val'] = "quality-factor-phi"