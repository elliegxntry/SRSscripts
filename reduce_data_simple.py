import numpy as np
import h5py
import os
import sys
sys.path.append("/home/amha5924/research/projects/gr-torus/scripts/modules")
from scripts import athena_read
from utils import *


config = "1.1.1-torus2_b-gz2"
specs = "a0beta500torBeta_br32x32x64rl2x2"
project = "gr-torus"

time_steps = np.arange(400, 405)

variables = ['press', 'Bcc1', 'Bcc2', 'Bcc3', 'rho', 'vel1', 'vel2', 'vel3']
vars_to_save = variables + ["invbeta", "qphi", "qtheta"]

data_load_path = "/run/media/amha5924/data/" + config + "_" + specs + "/data/raw/"
data_save_path = "/run/media/amha5924/data/" + config + "_" + specs + "/data/reduced/"

mp_path = data_save_path + "slices_midplane.hdf5"
vert_path = data_save_path + "slices_vertical.hdf5"
vertAvg_path = data_save_path + "slices_vertical-average.hdf5"
ext_path = data_save_path + "extrema.p"
consts_path = data_save_path + "constants.hdf5"

if not os.path.isdir(data_save_path):
    os.path.mkdir(data_save_path)

wrote_constants = False

for time in time_steps:
    # get raw data
    fname = data_load_path + config + "_" + specs + ".prim.{:05d}.athdf".format(time)
    print("-----------------------------------")
    print("Loading time step {}".format(time))
    with h5py.File(fname, 'r') as f:
        level = f.attrs['MaxLevel']
        data_t = athena_read.athdf(fname,
                                   quantities=variables, level=level)
    print("Done.")
    if not wrote_constants:
        save_constants(consts_path, data_t, False)
        wrote_constants = True

    for q in vars_to_save:
        # calculate the variables first
        computed_data = calculate_qdata(q, data_t)

        computed_mp = get_slice(computed_data, True, False)
        computed_vert = get_slice(computed_data, False, False)
        computed_vertAvg = get_slice(computed_data, False, True)

        # now save to path
        save_to_hdf5(mp_path, q, time, computed_mp, True)
        save_to_hdf5(vert_path, q, time, computed_vert, True)
        save_to_hdf5(vertAvg_path, q, time, computed_vertAvg, True)

        # Calculate extrema and save
        header = ["vol_min", "vol_max", "mp_min", "mp_max",
                  "vert_min", "vert_max", "vertAvg_min", "vertAvg_max"]
        [vmin, vmax] = [np.min(computed_data), np.max(computed_data)]
        [mpmin, mpmax] = [np.min(computed_mp), np.max(computed_mp)]
        [vertmin, vertmax] = [np.min(computed_vert), np.max(computed_vert)]
        [vertAvgmin, vertAvgmax] = [np.min(computed_vertAvg), np.max(computed_vertAvg)]
        extrema_data = [vmin, vmax, mpmin, mpmax, vertmin, vertmax, vertAvgmin, vertAvgmax]
        save_extrema(ext_path, q, time, extrema_data, False)

