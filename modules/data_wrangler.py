import gc
import csv
import sys
import numpy as np
import pandas as pd
import h5py
import os
import athena_read
import pickle
import glob
from collections import Mapping, Container
from matplotlib.animation import FFMpegWriter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
matplotlib.rc('text', usetex=True)
#plt.rcParams['animation.ffmpeg_path'] = '/home1/06165/ahankla/ffmpeg_sources/ffmpeg/ffmpeg'


# -------------- TO DO: -------------------
# - enable slicing data at arbitrary phi, theta
# - enable vertical slicing
# - proper error handling
# - grid for slice plotting (will need to wrap around data too)
# - extrema for phi-avg
# - write get_coord_ranges: call load, add if don't have all

class data_wrangler():

    def __init__(self, project, config, specs, to_load=False, **kwargs):
        """Initial load in. Loads in all time slices
)data_load_fp=None, fig_save_fp=None, data_save_fp=None, projhome=None):"""

        # -------------
        # Initialize
        # -------------
        self.config = config
        self.project = project
        self.specs = specs
        self.variables = ['press', 'Bcc1', 'Bcc2', 'Bcc3', 'rho', 'vel1', 'vel2', 'vel3']
        self.varnames = {'press':'$P$', 'rho':r'$\rho$', 'vel1':r'$v_r$', 'vel2':r'$v_\theta$', 'vel3':r'$v_\phi$', 'Bcc1': r'$B_r$', 'Bcc2': r'$B_\theta$', 'Bcc3': r'$B_\phi$', 'press-diff':r'$P(t)-P(0)$', 'rho-diff':r'$\rho(t)-\rho(0)$', 'Bcc3-diff':r'$B_\phi(t)-B_\phi(0)$', 'vel3-diff':r'$v_\phi(t)-v_\phi(0)$', 'mdot': r'$\dot M$', 'rhovr':r'$\rho v_r$', 'invbeta':r'$\beta^{-1}$', 'invbeta1':r'$\beta_r^{-1}$', 'invbeta2':r'$\beta_\theta^{-1}$', 'invbeta3':r'$\beta_\phi^{-1}$', 'Qphi':r'$Q_\phi$', 'Qtheta':r'$Q_\theta$', 'Rxrad':r'$T_{r\phi}^{\mathrm{Rey}}$', 'Rxpol':r'$T_{\theta\phi}^{\mathrm{Rey}}$', 'Mxrad':'$T_{r\phi}^{\mathrm{Max}}$', 'Mxpol':r'$T_{\theta\phi}^{\mathrm{Max}}$', 'alpha':r'$\frac{T_{r\phi}}{P_{\mathrm{gas}}}$'}
        self.qvars_diff = ["Bcc3-diff", "rho-diff", "press-diff", "vel3-diff"]
        self.qvars_calc = ["rhovr", "invbeta", "invbeta1", "invbeta2", "invbeta3", "Qtheta", "Qphi", "Rxrad", "Rxpol", "Mxrad", "Mxpol"] # XX to add: ...

        self.x1v = None; self.x2v = None; self.x3v = None
        self.x1f = None; self.x2f = None; self.x3f = None
        self.nx1 = None; self.nx2 = None; self.nx3 = None
        self.x1s_grid = None; self.x2s_grid = None; self.x3s_grid = None
        self.dphi_mp = None; self.dphi_vert = None;
        self.dtheta_mp = None; self.dtheta_vert = None;
        self.name_parts = None
        self.known_eidstr = ["vol", "x2s", "x3s"] # also defined in init
        self.qshear = None; self.Omega0 = None; self.pres = None; self.gamma = None;
        self.cs2 = None; self.level = 0;
        self.dfloor = sys.float_info.epsilon; self.pfloor = self.dfloor

        # -------------
        # Set/create file paths
        # ------------
        if 'projhome' not in kwargs:
            self.projhome = os.path.expanduser("~")+"/research/projects/"+self.project+"/"
        else:
            self.projhome = kwargs['projhome']
        if not os.path.exists(self.projhome):
            print("creating project home directory " + self.projhome)
            #os.makedirs(self.projhome)
        #
        if 'fig_save_fp' not in kwargs:
            self.fig_save_fp = (self.projhome+"figures/" +
            self.config + "_" + self.specs + "/")
        else:
            self.fig_save_fp = kwargs['fig_save_fp']
        if not os.path.exists(self.fig_save_fp):
            os.makedirs(self.fig_save_fp)
        #
        if 'data_save_fp' not in kwargs:
            self.data_save_fp = self.projhome + "data/reduced/"+self.config+"_"+self.specs+"/"
        else:
            self.data_save_fp = kwargs['data_save_fp']
        if not os.path.exists(self.data_save_fp):
            print("creating data save path "+self.data_save_fp)
            os.makedirs(self.data_save_fp)
        #
        if 'data_load_fp' not in kwargs:
            self.data_load_fp = self.projhome+"data/raw/"+self.config+"_"+self.specs+"/"
        else:
            self.data_load_fp = kwargs['data_load_fp']
        if not os.path.exists(self.data_load_fp):
            print("Error: file path to load data from doesn't exist: " + self.data_load_fp)

        if 'input_fp' not in kwargs:
            self.input_fp = self.projhome + "inputs/"+"athinput." + self.config + "_" + self.specs
        else:
            self.input_fp = kwargs['input_fp']
        if not os.path.exists(self.input_fp):
            print("Error: file path to input file doesn't exist: " + self.input_fp)
        self.consts_path = self.data_save_fp + "constants.hdf5"
        self.extrema_path = self.data_save_fp + "extrema.csv"
        self.ygforce_path = self.data_save_fp + "ygforce.csv"
        self.hst_path = glob.glob(self.data_load_fp + '*.hst')[0]
        self.mp_slice_path = self.data_save_fp + "slices_mp.hdf5"
        self.vert_slice_path = self.data_save_fp + "slices_vert.hdf5"
        self.vert_slice_aAvg_path = self.data_save_fp + "slices_vert_aAvg.hdf5"
        # self.times = self.find_all_t(to_load)

        # ----------------
        print(self.projhome)
        print(self.data_save_fp)
        print(self.data_load_fp)
        print(self.input_fp)
        print(self.hst_path)

    #------------------------------------------------------------------------
    def set_sim_params(self, to_load=False):
        """ read input file to get simulations parameters """
        with open(self.input_fp, 'r') as f:
            lines = f.readlines()
        lines.reverse()
        import difflib
        lines = np.array(lines)
        get_dims = False
        for i in np.arange(0, lines.size):
            line = lines[i].replace(" ", "")
            line = line.replace("#", "=")
            line = line.strip("\n")
            line = line.split("=")
            # physical parameters
            if line[0] == 'p0_over_r0':
                self.pres = np.float(line[1])
            elif line[0] == 'qshear':
                self.qshear = np.float(line[1])
            elif line[0] == 'Omega0':
                self.Omega0 = np.float(line[1])
            elif line[0] == 'gamma':
                self.gamma = np.float(line[1])
            elif line[0] == 'beta':
                self.beta = np.float(line[1])
            elif line[0] == 'iprob':
                self.iprob = np.int(line[1])
            elif line[0] == 'r1':
                self.r1 = np.float(line[1])
            elif line[0] == 'r2':
                self.r2 = np.float(line[1])
            elif line[0] == 'rm':
                self.rm = np.float(line[1])
            elif line[0] == 'angle':
                self.angle = np.float(line[1])
            # numerical parameters
            elif line[0] == 'x1rat':
                self.x1rat = np.float(line[1])
            elif line[0] == 'x3max':
                self.x3max = np.float(line[1])
                get_dims = True
            elif line[0] == 'x3min':
                self.x3min = np.float(line[1])
            elif line[0] == 'x2max':
                self.x2max = np.float(line[1])
            elif line[0] == 'x2min':
                self.x2min = np.float(line[1])
            elif line[0] == 'x1max':
                self.x1max = np.float(line[1])
            elif line[0] == 'x1min':
                self.x1min = np.float(line[1])
            elif line[0] == 'level':
                if np.int(line[1]) > self.level:
                    self.level = np.int(line[1])
            elif line[0] == 'dfloor':
                self.dfloor = np.float(line[1])
            elif line[0] == 'pfloor':
                self.pfloor = np.float(line[1])
            elif get_dims and line[0] == 'nx3':
                self.nx3 = np.int(line[1])*2**self.level
            elif get_dims and line[0] == 'nx2':
                self.nx2 = np.int(line[1])*2**self.level
            elif get_dims and line[0] == 'nx1':
                self.nx1 = np.int(line[1])*2**self.level

        # -- Defaults --
        if self.qshear is None:
            self.qshear = 1.5
        if self.Omega0 is None:
            self.Omega0 = 1.0
        if self.pres is None:
            self.pres = 0.1
        if self.gamma is None:
            self.gamma = 1.0
        if self.x1rat == 1.0:
            self.logr = False
        else:
            self.logr = True

        self.cs2 = self.gamma*self.pres
        self.H = np.sqrt(self.cs2)/self.Omega0 # at r=1

        self.times = self.find_all_t(to_load)
        self.maxt = max(self.tvals) # note: is code-time, not index. self.tvals defined in find_all_t()
        self.load_coord_ranges()

        try:
            # calculate radius where b field switches sign
            if self.iprob == 3 and self.rm is None:
                ra = 0.5*(self.r1**1.5 + self.r2**1.5)
                self.rm = ra**(2./3.)
        except: pass

        print("\n---- Simulation Parameters ----")
        print("Radial range: [{}, {}]".format(self.x1min, self.x1max))
        print("Theta range: [{}, {}]".format(self.x2min, self.x2max))
        print("Azimuthal range: [{}, {}]".format(self.x3min, self.x3max))
        print("Radial Logarithmic ratio: {}".format(self.x1rat))
        print("Base resolution: [{}, {}, {}]".format(self.nx1, self.nx2, self.nx3))
        print("Max refinement level: {}".format(self.level))
        print("Refined (total) resolution: [{}, {}, {}]".format(self.nx1*2**self.level, self.nx2*2**self.level, self.nx3*2**self.level))
        print("Density floor: {}".format(self.dfloor))
        print("Pressure floor: {}".format(self.pfloor))
        print("q: {} ".format(self.qshear))
        print("Omega0: {}".format(self.Omega0))
        print("pres: {}".format(self.pres))
        print("gamma: {}".format(self.gamma))
        print("scale height H: {}".format(self.H))
        print("max time: {}".format(self.maxt))
        try:
            print("beta: {}".format(self.beta))
        except: pass
        try:
            print("r1: {}".format(self.r1))
        except: pass
        try:
            print("rm: {}".format(self.rm))
        except: pass
        try:
            print("r2: {}".format(self.r2))
        except: pass
        try:
            print("angle: {}".format(self.angle))
        except: pass
        print("------------------------------\n")


    #------------------------------------------------------------------------
    def load_coord_ranges(self, force_load=False):
        """ Load coordinate ranges from existing constants file. Return None if doesn't exist.
        XX could need to be edited?
        """
        towrite = False
        with h5py.File(self.consts_path, 'r') as f:
            try:
                self.x1f = f["x1f"][:]
                self.x2f = f["x2f"][:]
                self.x3f = f["x3f"][:]
                self.x1v = f["x1v"][:]
                self.x2v = f["x2v"][:]
                self.x3v = f["x3v"][:]
                self.x1s_grid = f["x1s_grid"][:]
                self.x2s_grid = f["x2s_grid"][:]
                self.x3s_grid = f["x3s_grid"][:]
                self.x2vcorr = f["x2vcorr"][:]
                self.x2fcorr = f["x2fcorr"][:]
                self.dtheta_vert = f["dtheta_vert"][:]
                self.dtheta_mp = f["dtheta_mp"][:]
                self.dphi_vert = f["dphi_vert"][:]
                self.dphi_mp = f["dphi_mp"][:]
                print("loaded all coord ranges")
            except:
                print("need to load coord ranges")
                t0data = self.load_q_data_raw_at_t(["rho"], 0)
                self.add_coord_ranges(t0data)
                towrite = True
        if force_load:
            print("need to load coord ranges")
            t0data = self.load_q_data_raw_at_t(["rho"], 0)
            self.add_coord_ranges(t0data)
            towrite = True

        if towrite:
            with h5py.File(self.consts_path, 'a') as f:
                # write to file
                try:
                    del f["x1v"], f["x2v"], f["x3v"]
                    del f["x1f"], f["x2f"], f["x3f"]
                    del f["x1s_grid"], f["x2s_grid"], f["x3s_grid"]
                    del f["x2fcorr"], f["x2vcorr"]
                    del f["dtheta_mp"], f["dphi_mp"]
                    del f["dtheta_vert"], f["dphi_vert"]
                except:
                    pass
                finally:
                    f.create_dataset("x1v", data=self.x1v)
                    f.create_dataset("x2v", data=self.x2v)
                    f.create_dataset("x3v", data=self.x3v)
                    f.create_dataset("x1f", data=self.x1f)
                    f.create_dataset("x2f", data=self.x2f)
                    f.create_dataset("x3f", data=self.x3f)
                    f.create_dataset("x2fcorr", data=self.x2fcorr)
                    f.create_dataset("x2vcorr", data=self.x2vcorr)
                    f.create_dataset("x1s_grid", data=self.x1s_grid)
                    f.create_dataset("x2s_grid", data=self.x2s_grid)
                    f.create_dataset("x3s_grid", data=self.x3s_grid)
                    f.create_dataset("dtheta_mp", data=self.dtheta_mp)
                    f.create_dataset("dtheta_vert", data=self.dtheta_vert)
                    f.create_dataset("dphi_mp", data=self.dphi_mp)
                    f.create_dataset("dphi_vert", data=self.dphi_vert)
        if self.nx1 is None:
            self.nx1 = self.x1v.size
            self.nx2 = self.x2v.size
            self.nx3 = self.x3v.size
        return [self.x1s_grid, self.x2s_grid, self.x3s_grid]

    #------------------------------------------------------------------------
    def add_coord_ranges(self, data_raw, force_overwrite=False):
        """ Add coordinate ranges and write to constants file
        XX could need to edit? Meshgrid in particular
        """
        # print("adding to coords")
        self.x1v = data_raw['x1v']
        self.x2v = data_raw['x2v']
        self.x3v = data_raw['x3v']
        self.x1f = data_raw['x1f']
        self.x2f = data_raw['x2f']
        self.x3f = data_raw['x3f']

        # Phi, theta don't go all the way to (2)pi, so extend it
        # see athena++ release's "plot_slice.py"
        theta_extended = np.concatenate((-self.x2v[0:1], self.x2v, 2.0 * np.pi - self.x2v[::-1],
                                        2.0 * np.pi + self.x2v[0:1]))
        theta_extended_corrected = theta_extended - 0.5 * (self.x2v[1] - self.x2v[0])
        self.x2vcorr = theta_extended_corrected
        theta_extended2 = np.concatenate((-self.x2f[0:1], self.x2f, 2.0 * np.pi - self.x2f[::-1],
                                        2.0 * np.pi + self.x2f[0:1]))
        theta_extended_corrected2 = theta_extended2 - 0.5 * (self.x2f[1] - self.x2f[0])
        self.x2fcorr = theta_extended_corrected2
        phi_extended = \
            np.concatenate((self.x3v[-1:] - 2.0 * np.pi, self.x3v, self.x3v[:1] + 2.0 * np.pi))
        phi_extended -= 0.5 * 2.0 * np.pi / self.nx3

        # --- calculate spherical slice grids ---
        # theta-phi plane (weird)
        self.x1s_grid = np.meshgrid(theta_extended_corrected, phi_extended)
        # r-phi plane (midplane slice)
        r_grid, phi_grid = np.meshgrid(self.x1v, phi_extended)
        self.x2s_grid = [r_grid*np.cos(phi_grid), r_grid*np.sin(phi_grid)]
        # r-theta plane (vertical slice)
        r_grid, theta_grid = np.meshgrid(self.x1v, theta_extended_corrected)
        self.x3s_grid = [r_grid*np.sin(theta_grid), r_grid*np.cos(theta_grid)]

        # --- calculate differences ---
        dtheta = (np.hstack((self.x2f, 0)) - np.hstack((0, self.x2f)))[1:-1]
        dtheta_ext = np.concatenate((dtheta[0:1], dtheta, dtheta[::-1], dtheta[0:1]))
        dphi = (np.hstack((self.x3f, 0)) - np.hstack((0, self.x3f)))[1:-1]
        dphi_ext = np.concatenate((dphi[-1:], dphi, dphi[1:]))

        # for vertical slices, take dphi[0]. Extend dtheta to every radius
        self.dphi_vert = dphi[0]*np.ones((self.x2v.size*2+2, self.x1v.size))
        self.dtheta_vert = np.tile(dtheta_ext, (self.x1v.size, 1)).transpose()

        # for midplane slices, take midplane value of dtheta. Extend dphi to every radius
        self.dtheta_mp = dtheta[int(self.x2v.size/4)]*np.ones((self.x3v.size+2, self.x1v.size))
        self.dphi_mp = np.tile(dphi_ext, (self.x1v.size, 1)).transpose()

        return [self.x1s_grid, self.x2s_grid, self.x3s_grid]

    #------------------------------------------------------------------------
    def get_name_parts(self):
        name_stem = glob.glob(self.data_load_fp + '*.out1.*.athdf')[0]
        self.name_parts = name_stem.split('.')
        return self.name_parts
    #------------------------------------------------------------------------
    def load_qs_data_hst(self, qs=None):
        data = athena_read.hst(self.hst_path)
        time = data.pop("time")
        if qs is not None:
            data = dict((q, data[q]) for q in qs)
        return [time, data]

    #------------------------------------------------------------------------
    def load_q_data_raw_at_t(self, q, t):
        """Load in variable(s) q from raw .hdf5 files at specific time t.

        Keyword arguments:
        t: the output file number, not the index, i.e. sim_data.out2.00010.hdf5 would be 10.
        q: a list of variables

        Do not write to file; simply load it in.
        """
        #print("loading raw data from .athdf files")
        #print(glob.glob(self.data_load_fp+'*out2.*.athdf'))

        # extract a representative name
        if self.name_parts is None:
            self.get_name_parts()
        # print(self.name_stem)
        ts = "{:05}".format(t)
        self.name_parts[-2] = ts
        name = '.'.join(self.name_parts)

        if not os.path.isfile(name):
            print("ERROR:")
            print("FILE " + name + " DOES NOT EXIST. ")
            print("EXITING...")
            return None
        print("Loading time step {}".format(t))
        with h5py.File(name, 'r') as f:
            level = f.attrs['MaxLevel']
            data_raw_q_t = athena_read.athdf(name,
                                                quantities=q, level=level)
        print("Done.")
        # print(data_raw_q_t.keys())
        # if "press" in data_raw_q_t.keys():
            # print(np.max(data_raw_q_t["press"]))


        # Extract coordinate ranges if not already saved
        # self.add_coord_ranges(data_raw_q_t, True)

        return data_raw_q_t


    #------------------------------------------------------------------------
    def extract_q_extrema_at_t(self, t, qs, data):
        """ only output extrema of variables qs at time t.

        Does not check if extrema exists in file.
        Does not write to file or anything.

        For use in both get_q_extrema_at_t and get_q_extrema_allt
        """
        lines = np.array([])
        for q in qs:
            print(q)
            vol_extrema = [np.nanmin(data[q]), np.nanmax(data[q])]

            x2s_data = self.extract_slice_data_at_t(data[q], "x2s")
            x3s_data = self.extract_slice_data_at_t(data[q], "x3s")

            x2s_extrema = [np.nanmin(x2s_data), np.nanmax(x2s_data)]
            x3s_extrema = [np.nanmin(x3s_data), np.nanmax(x3s_data)]


            self.known_eidstr = ["vol", "x2s", "x3s"] # also defined in init

            # write line to CSV
            row = np.array([int(t), q, vol_extrema[0], vol_extrema[1], x2s_extrema[0], x2s_extrema[1], x3s_extrema[0], x3s_extrema[1]])
            if lines.size < 1:
                lines = np.reshape(row, (1, row.size))
            else:
                lines = np.vstack((lines, row))

        # create pandas dataframe
        if lines.ndim == 1:
            extrema = pd.DataFrame(lines, columns=self.header)
            extrema.set_index(['time', 'variable'])
        else:
            tuples = list(zip(*[lines[:, 0], lines[:, 1]]))
            idx = pd.MultiIndex.from_tuples(tuples, names=['time', 'variable'])

            extrema = pd.DataFrame(lines[:, 2:], index=idx, columns=self.header[2:])
            extrema = extrema.astype(np.float64)
            # print(extrema.dtypes)

        return extrema

    # ---------------------------------------------------
    def save_qs_slices(self, qs=None, **kwargs):
        """ save both vertical and midplane slices for variable q"""
        if qs is None: qs = self.variables
        force_load = kwargs.get("force_load", False)
        t0data = kwargs.get("t0data", None)

        for t in self.tinds:
            print("Loading slices from t={}".format(t))
            vertdata = None
            mpdata = {}; vertdata = {}
            with h5py.File(self.mp_slice_path, 'a') as f:
                qsdata = None
                for q in qs:
                    if q not in f:
                        f.create_group(q)
                    if str(t) not in f[q]:
                        if qsdata is None:
                            qsdata = self.load_q_data_raw_at_t(qs, t)
                        # slice data
                        mpdata[q] = self.extract_slice_data_at_t(qsdata[q], "x2s")
                        vertdata[q] = self.extract_slice_data_at_t(qsdata[q], "x3s")
                        f[q].create_dataset(str(t), data=mpdata[q])
                        # f[q][str(t)].attrs["Time"] = qsdata.time

                    # do diff variables
                    if t0data is not None and q in self.qvars_diff:
                        qd = q + "-diff"
                        if qd not in f:
                            f.create_group(qd)
                        if str(t) not in f[qd]:
                            if q not in mpdata.keys() or q not in vertdata.keys():
                                if qsdata is None:
                                    qsdata = self.load_q_data_raw_at_t(qs, t)
                                mpdata[q] = self.extract_slice_data_at_t(qsdata[q], "x2s")
                                vertdata[q] = self.extract_slice_data_at_t(qsdata[q], "x3s")
                            # slice data
                            mpdata[qd] = mpdata[q] - self.extract_slice_data_at_t(t0data[q], "x2s")
                            vertdata[qd] = vertdata[q] - self.extract_slice_data_at_t(t0data[q], "x3s")
                            f[qd].create_dataset(str(t), data=mpdata[qd])
                            f[qd][str(t)].attrs["Time"] = qsdata.time

            with h5py.File(self.vert_slice_path, 'a') as f:
                if q not in f:
                    f.create_group(q)
                if str(t) not in f[q]:
                    if vertdata[q] is None:
                        if qsdata is None:
                            qsdata = self.load_q_data_raw_at_t([q], t)[q]
                        # slice data
                        vertdata[q] = self.extract_slice_data_at_t(qsdata, "x3s")
                    f[q].create_dataset(str(t), data=vertdata[q])
                    f[q][str(t)].attrs["Time"] = qsdata.time

                # do diff variables
                if t0data is not None and q in self.qvars_diff:
                    qd = q + "-diff"
                    if qd not in f:
                        f.create_group(qd)
                    if str(t) not in f[qd]:
                        if vertdata[qd] is None:
                            if q not in vertdata.keys() or q not in mpdata.keys():
                                if qsdata is None:
                                    qsdata = self.load_q_data_raw_at_t(qs, t)
                                    # slice data
                                    vertdata[q] = self.extract_slice_data_at_t(qsdata[q], "x3s")
                            vertdata[qd] = vertdata[q] - self.extract_slice_data_at_t(t0data[q], "x3s")
                        f[qd].create_dataset(str(t), data=mpdata[qd])
                        f[qd][str(t)].attrs["Time"] = qsdata.time


    # ---------------------------------------------------
    def save_qs_slice_aAvg(self, qs=None, **kwargs):
        """ save vertical slices for variables qs, azimuthally averaged"""
        if qs is None: qs = self.variables
        force_load = kwargs.get("force_load", False)
        t0data = kwargs.get("t0data", None)

        for t in self.tinds:
            print("Loading slices from t={}".format(t))
            mpdata = {}; vertdata = {}
            with h5py.File(self.vert_slice_aAvg_path, 'a') as f:
                qsdata = None
                for q in qs:
                    if q not in f:
                        f.create_group(q)
                    if str(t) not in f[q]:
                        if qsdata is None:
                            qsdata = self.load_q_data_raw_at_t(qs, t)
                        # average data
                        vertdata[q] = np.mean(qsdata[q], axis=0)
                        self.extract_slice_data_at_t(qsdata[q], "x3s")
                        f[q].create_dataset(str(t), data=vertdata[q])

                    # do diff variables
                    if t0data is not None and q in self.qvars_diff:
                        qd = q + "-diff"
                        if qd not in f:
                            f.create_group(qd)
                        if str(t) not in f[qd]:
                            if q not in vertdata.keys():
                                if qsdata is None:
                                    qsdata = self.load_q_data_raw_at_t(qs, t)
                                vertdata[q] = self.extract_slice_data_at_t(qsdata[q], "x3s")
                            # average data
                            vertdata[qd] = np.mean(vertdata[q] - self.extract_slice_data_at_t(t0data[q], "x3s"), axis=0)
                            f[qd].create_dataset(str(t), data=mpdata[qd])


    # ---------------------------------------------------
    def reduce_data(self, qs=None, **kwargs):
        """ Reduce data, i.e. take midplane/vertical slices, azimuthally average, get extrema, etc.

        XX to do: only add diff vars
        """

        # Defaults/set-up
        if qs is None: qs = self.variables
        force_load = kwargs.get("force_load", False)
        write_header = False

        mpSlices = kwargs.get("mpSlices", True)
        vertSlices = kwargs.get("vertSlices", True)
        aAvgSlices = kwargs.get("aAvgSlices", True)
        calcVars = kwargs.get("calcVars", True)
        extremaLoad = kwargs.get("extrema", True)

        if extremaLoad:
            self.header = ["time", "variable", "vol_min", "vol_max", "x2s_min", "x2s_max", "x3s_min", "x3s_max"]
            # if file doesn't exist, export header. if it does, try to load data
            if not os.path.isfile(self.extrema_path):
                write_header = header
                extrema = pd.DataFrame()
            else:
                extrema = pd.read_csv(self.extrema_path, index_col=["time", "variable"])
                extrema.drop_duplicates(subset=None, keep='first', inplace=True)


        # t0data = self.load_q_data_raw_at_t([var.replace('-diff', '') for var in self.qvars_diff], 0)
        t0data = self.load_q_data_raw_at_t(self.variables, 0)

        # Loop through every time step
        for t in np.range(0,1): #self.tinds:
            qsdata = None
            print("------------------------------")
            print("Reducing data for time step {}".format(t))

            mpdata = {}; vertdata = {}; avgdata = {}

            # Midplane slices
            if mpSlices:
                print("Midplane slices")
                with h5py.File(self.mp_slice_path, 'a') as f:
                    for q in qs:
                        if q not in f or force_load:
                            try:
                                f.create_group(q)
                            except: pass
                        if str(t) not in f[q] or force_load:
                            if qsdata is None:
                                qsdata = self.load_q_data_raw_at_t(qs, t)

                            # slice data
                            mpdata[q] = self.extract_slice_data_at_t(qsdata[q], "x2s")
                            vertdata[q] = self.extract_slice_data_at_t(qsdata[q], "x3s")
                            try:
                                del f[q][str(t)]
                            except: pass
                            f[q].create_dataset(str(t), data=mpdata[q])
                            # f[q][str(t)].attrs["Time"] = qsdata.time

                        # do diff variables
                        if t0data is not None and q+"-diff" in self.qvars_diff:
                            qd = q + "-diff"
                            if qd not in f:
                                f.create_group(qd)
                            if str(t) not in f[qd]:
                                if q not in mpdata.keys() or q not in vertdata.keys():
                                    if qsdata is None:
                                        qsdata = self.load_q_data_raw_at_t(qs, t)
                                    mpdata[q] = self.extract_slice_data_at_t(qsdata[q], "x2s")
                                    vertdata[q] = self.extract_slice_data_at_t(qsdata[q], "x3s")
                                # slice data
                                mpdata[qd] = mpdata[q] - self.extract_slice_data_at_t(t0data[q], "x2s")
                                vertdata[qd] = vertdata[q] - self.extract_slice_data_at_t(t0data[q], "x3s")
                                f[qd].create_dataset(str(t), data=mpdata[qd])
                                f[qd][str(t)].attrs["Time"] = qsdata.time


            if vertSlices:
                print("Vertical slices")
                with h5py.File(self.vert_slice_path, 'a') as f:
                    for q in qs:
                        if q not in f:
                            f.create_group(q)
                        if str(t) not in f[q]:
                            if q not in vertdata.keys():
                                if qsdata is None:
                                    qsdata = self.load_q_data_raw_at_t(qs, t)
                                # slice data
                                vertdata[q] = self.extract_slice_data_at_t(qsdata[q], "x3s")
                            f[q].create_dataset(str(t), data=vertdata[q])
                            f[q][str(t)].attrs["Time"] = qsdata.time

                        # do diff variables
                        if t0data is not None and q+"-diff" in self.qvars_diff:
                            qd = q + "-diff"
                            if qd not in f:
                                f.create_group(qd)
                            if str(t) not in f[qd]:
                                print("getting diff data")
                                if qd not in vertdata.keys():
                                    if q not in vertdata.keys():
                                        if qsdata is None:
                                            qsdata = self.load_q_data_raw_at_t(qs, t)
                                        vertdata[q] = self.extract_slice_data_at_t(qsdata[q], "x3s")
                                    vertdata[qd] = vertdata[q] - self.extract_slice_data_at_t(t0data[q], "x3s")
                                f[qd].create_dataset(str(t), data=vertdata[qd])
                                f[qd][str(t)].attrs["Time"] = qsdata.time

            if aAvgSlices:
                print("Averaged Vertical slices")
                with h5py.File(self.vert_slice_aAvg_path, 'a') as f:
                    for q in qs:
                        if q not in f:
                            f.create_group(q)
                        if str(t) not in f[q]:
                            print("getting data for {}, t {}".format(q, t))
                            if qsdata is None:
                                qsdata = self.load_q_data_raw_at_t(qs, t)
                            # average data
                            avgdata[q] = self.extract_slice_data_at_t(qsdata[q], "x3s", avg=True)
                            f[q].create_dataset(str(t), data=avgdata[q])

                        # do diff variables
                        if t0data is not None and q+"-diff" in self.qvars_diff:
                            qd = q + "-diff"
                            if qd not in f:
                                f.create_group(qd)
                            if str(t) not in f[qd]:
                                if qsdata is None:
                                    qsdata = self.load_q_data_raw_at_t(qs, t)
                                # average data
                                avgdata[qd] = self.extract_slice_data_at_t(qsdata[q]-t0data[q], "x3s", avg=True)
                                f[qd].create_dataset(str(t), data=avgdata[qd])

            if calcVars:
                rhovrVert = None; invbetasliceVert = None;
                print("Calculating variables")
                with h5py.File(self.vert_slice_aAvg_path, 'a') as f:
                    qv = "rhovr"
                    if qv not in f:
                        f.create_group(qv)
                    if (force_load and "rhovr" in qs) or str(t) not in f[qv]:
                        print("getting data for {}, t {}".format(qv, t))
                        if qsdata is None or "rho" not in qsdata.keys() or "vel1" not in qsdata.keys():
                            qsdata = self.load_q_data_raw_at_t(["rho", "press", "Bcc1", "Bcc2", "Bcc3", "vel1", "press"], t)
                        if aAvgSlices:
                            rhovr = self.extract_slice_data_at_t(np.multiply(qsdata["rho"], qsdata["vel1"]), "x3s", avg=True)
                        if vertSlices:
                            rhovrVert = self.extract_slice_data_at_t(np.multiply(qsdata["rho"], qsdata["vel1"]), "x3s")
                        if mpSlices:
                            rhovrMp = self.extract_slice_data_at_t(np.multiply(qsdata["rho"], qsdata["vel1"]), "x2s")
                        try:
                            del f[qv][str(t)]
                        except: pass
                        f[qv].create_dataset(str(t), data=rhovr)
                    qv = "invbeta"
                    qv1 = "invbeta1"; qv2 = "invbeta2"; qv3 = "invbeta3"
                    if qv not in f:
                        f.create_group(qv)
                    if qv1 not in f:
                        f.create_group(qv1)
                    if qv2 not in f:
                        f.create_group(qv2)
                    if qv3 not in f:
                        f.create_group(qv3)
                    if (force_load and qv in qs) or str(t) not in f[qv] or str(t) not in f[qv1] or str(t) not in f[qv2] or str(t) not in f[qv3]:
                        print("getting data for {}, t {}".format(qv, t))
                        if qsdata is None or "Bcc1" not in qsdata.keys() or "Bcc2" not in qsdata.keys() or "Bcc3" not in qsdata.keys() or "press" not in qsdata.keys():
                            qsdata = self.load_q_data_raw_at_t(["Bcc1", "Bcc2", "Bcc3", "press", "rho", "vel3"], t)
                        bcc12 = np.square(qsdata["Bcc1"])
                        bcc22 = np.square(qsdata["Bcc2"])
                        bcc32 = np.square(qsdata["Bcc3"])
                        pmag = bcc12+bcc22+bcc32
                        invbeta1 = bcc12/2/qsdata["press"]
                        invbeta2 = bcc22/2/qsdata["press"]
                        invbeta3 = bcc32/2/qsdata["press"]
                        invbeta = pmag/2/qsdata["press"]
                        invbetaslice = self.extract_slice_data_at_t(invbeta, "x3s", avg=True)
                        invbeta1slice = self.extract_slice_data_at_t(invbeta1, "x3s", avg=True)
                        invbeta2slice = self.extract_slice_data_at_t(invbeta2, "x3s", avg=True)
                        invbeta3slice = self.extract_slice_data_at_t(invbeta3, "x3s", avg=True)
                        invbetasliceVert = self.extract_slice_data_at_t(invbeta, "x3s")
                        invbeta1sliceVert = self.extract_slice_data_at_t(invbeta1, "x3s")
                        invbeta2sliceVert = self.extract_slice_data_at_t(invbeta2, "x3s")
                        invbeta3sliceVert = self.extract_slice_data_at_t(invbeta3, "x3s")
                        invbetasliceMp = self.extract_slice_data_at_t(invbeta, "x2s")
                        invbeta1sliceMp = self.extract_slice_data_at_t(invbeta1, "x2s")
                        invbeta2sliceMp = self.extract_slice_data_at_t(invbeta2, "x2s")
                        invbeta3sliceMp = self.extract_slice_data_at_t(invbeta3, "x2s")
                        try:
                            del f[qv][str(t)]
                        except: pass
                        f[qv].create_dataset(str(t), data=invbetaslice)
                        try:
                            del f[qv1][str(t)]
                        except: pass
                        f[qv1].create_dataset(str(t), data=invbeta1slice)
                        try:
                            del f[qv2][str(t)]
                        except: pass
                        f[qv2].create_dataset(str(t), data=invbeta2slice)
                        try:
                            del f[qv3][str(t)]
                        except: pass
                        f[qv3].create_dataset(str(t), data=invbeta3slice)
                    qv2 = "Qtheta"; qv3 = "Qphi"
                    if self.dphi_vert is None or self.dphi_mp is None:
                        self.load_coord_ranges(True)
                    if qv2 not in f:
                        f.create_group(qv2)
                    if qv3 not in f:
                        f.create_group(qv3)
                    if (force_load and qv2 in qs) or str(t) not in f[qv2] or str(t) not in f[qv3]:
                        print("getting data for {} and {}, t {}".format(qv2, qv3, t))
                        if qsdata is None or "rho" not in qsdata.keys() or "vel3" not in qsdata.keys() or "Bcc3" not in qsdata.keys() or "Bcc2" not in qsdata.keys():
                            qsdata = self.load_q_data_raw_at_t(["rho", "vel3", "Bcc3", "Bcc2"], t)
                        phidata = np.abs(2*np.pi*np.divide(qsdata["Bcc3"], np.multiply(qsdata["vel3"], np.sqrt(qsdata["rho"]))))
                        thetadata = np.abs(2*np.pi*np.divide(qsdata["Bcc2"], np.multiply(qsdata["vel3"], np.sqrt(qsdata["rho"]))))
                        if aAvgSlices:
                            qphi = self.extract_slice_data_at_t(phidata, "x3s", avg=True)
                            qphi = np.divide(qphi, self.dphi_vert)
                            qtheta = self.extract_slice_data_at_t(thetadata, "x3s", avg=True)
                            qtheta = np.divide(qtheta, self.dtheta_vert)
                        if vertSlices:
                            qphiVert = self.extract_slice_data_at_t(phidata, "x3s")
                            qphiVert = np.divide(qphiVert, self.dphi_vert)
                            qthetaVert = self.extract_slice_data_at_t(thetadata, "x3s")
                            qthetaVert = np.divide(qthetaVert, self.dtheta_vert)
                        if mpSlices:
                            qphiMp = self.extract_slice_data_at_t(phidata, "x2s")
                            qphiMp = np.divide(qphiVert, self.dphi_mp)
                            qthetaMp = self.extract_slice_data_at_t(thetadata, "x2s")
                            qthetaMp = np.divide(qthetaVert, self.dtheta_mp)
                        try:
                            del f[qv2][str(t)]
                            del f[qv3][str(t)]
                        except: pass
                        f[qv2].create_dataset(str(t), data=qtheta)
                        f[qv3].create_dataset(str(t), data=qphi)
                    qvrr = "Rxrad"; qvmr = "Mxrad";
                    qvrp = "Rxpol"; qvmp = "Mxpol"
                    if qvrr not in f:
                        f.create_group(qvrr)
                    if qvmr not in f:
                        f.create_group(qvmr)
                    if qvrp not in f:
                        f.create_group(qvrp)
                    if qvmp not in f:
                        f.create_group(qvmp)
                    if (force_load and qvrr in qs) or str(t) not in f[qvrr] or str(t) not in f[qvmr] or str(t) not in f[qvrp] or str(t) not in f[qvmp]:
                        print("getting data for {}, t {}".format(qvrr, t))
                        if qsdata is None or "Bcc1" not in qsdata.keys() or "Bcc2" not in qsdata.keys() or "Bcc3" not in qsdata.keys() or "rho" not in qsdata.keys():
                            qsdata = self.load_q_data_raw_at_t(["Bcc1", "Bcc2", "Bcc3", "rho", "vel1", "vel2", "vel3"], t)

                        # calculate difference from mean
                        dvphi = qsdata["vel3"] - np.mean(qsdata["vel3"], axis=0)
                        # calculate stresses
                        trphi = np.multiply(qsdata["rho"], np.multiply(dvphi, qsdata["vel1"]))
                        ttphi = np.multiply(qsdata["rho"], np.multiply(dvphi, qsdata["vel2"]))
                        mrphi = -np.multiply(qsdata["Bcc1"], qsdata["Bcc3"])
                        mtphi = -np.multiply(qsdata["Bcc2"], qsdata["Bcc3"])

                        trphiAvg = self.extract_slice_data_at_t(trphi, "x3s", avg=True)
                        ttphiAvg = self.extract_slice_data_at_t(ttphi, "x3s", avg=True)
                        mrphiAvg = self.extract_slice_data_at_t(mrphi, "x3s", avg=True)
                        mtphiAvg = self.extract_slice_data_at_t(mtphi, "x3s", avg=True)
                        trphiVert = self.extract_slice_data_at_t(trphi, "x3s")
                        ttphiVert = self.extract_slice_data_at_t(ttphi, "x3s")
                        mrphiVert = self.extract_slice_data_at_t(mrphi, "x3s")
                        mtphiVert = self.extract_slice_data_at_t(mtphi, "x3s")
                        trphiMp = self.extract_slice_data_at_t(trphi, "x2s")
                        ttphiMp = self.extract_slice_data_at_t(ttphi, "x2s")
                        mrphiMp = self.extract_slice_data_at_t(mrphi, "x2s")
                        mtphiMp = self.extract_slice_data_at_t(mtphi, "x2s")
                        try:
                            del f[qvrr][str(t)]
                        except: pass
                        f[qvrr].create_dataset(str(t), data=trphiAvg)
                        try:
                            del f[qvrp][str(t)]
                        except: pass
                        f[qvrp].create_dataset(str(t), data=ttphiAvg)
                        try:
                            del f[qvmr][str(t)]
                        except: pass
                        f[qvmr].create_dataset(str(t), data=mrphiAvg)
                        try:
                            del f[qvmp][str(t)]
                        except: pass
                        f[qvmp].create_dataset(str(t), data=mtphiAvg)

                if vertSlices:
                    with h5py.File(self.vert_slice_path, 'a') as f:
                        if "rhovr" not in f:
                            f.create_group("rhovr")
                        if "invbeta" not in f:
                            f.create_group("invbeta")
                        if "invbeta1" not in f:
                            f.create_group("invbeta1")
                        if "invbeta2" not in f:
                            f.create_group("invbeta2")
                        if "invbeta3" not in f:
                            f.create_group("invbeta3")
                        if "Qphi" not in f:
                            f.create_group("Qphi")
                        if "Qtheta" not in f:
                            f.create_group("Qtheta")
                        if "Rxrad" not in f:
                            f.create_group("Rxrad")
                        if "Rxpol" not in f:
                            f.create_group("Rxpol")
                        if "Mxrad" not in f:
                            f.create_group("Mxrad")
                        if "Mxpol" not in f:
                            f.create_group("Mxpol")


                        if (force_load and "rhovr" in qs) or str(t) not in f["rhovr"]:
                            try:
                                del f["rhovr"][str(t)]
                            except: pass
                            if qsdata is None or "rho" not in qsdata.keys() or "vel1" not in qsdata.keys():
                                qsdata = self.load_q_data_raw_at_t(["rho", "press", "Bcc1", "Bcc2", "Bcc3", "vel1", "press"], t)
                            if rhovrVert is None:
                                rhovrVert = self.extract_slice_data_at_t(np.multiply(qsdata["rho"], qsdata["vel1"]), "x3s")
                            f["rhovr"].create_dataset(str(t), data=rhovrVert)
                        if (force_load and qv in qs) or str(t) not in f["invbeta"]:
                            try:
                                del f["invbeta"][str(t)]
                            except: pass
                            if qsdata is None or "Bcc1" not in qsdata.keys() or "Bcc2" not in qsdata.keys() or "Bcc3" not in qsdata.keys() or "press" not in qsdata.keys():
                                qsdata = self.load_q_data_raw_at_t(["Bcc1", "Bcc2", "Bcc3", "press", "rho", "vel3"], t)

                            if invbetasliceVert is None:
                                bcc12 = np.square(qsdata["Bcc1"])
                                bcc22 = np.square(qsdata["Bcc2"])
                                bcc32 = np.square(qsdata["Bcc3"])
                                pmag = bcc12+bcc22+bcc32
                                invbeta1 = bcc12/2/qsdata["press"]
                                invbeta2 = bcc22/2/qsdata["press"]
                                invbeta3 = bcc32/2/qsdata["press"]
                                invbeta = pmag/2/qsdata["press"]
                                invbetasliceVert = self.extract_slice_data_at_t(invbeta, "x3s")
                                invbeta1sliceVert = self.extract_slice_data_at_t(invbeta1, "x3s")
                                invbeta2sliceVert = self.extract_slice_data_at_t(invbeta2, "x3s")
                                invbeta3sliceVert = self.extract_slice_data_at_t(invbeta3, "x3s")
                                invbetasliceMp = self.extract_slice_data_at_t(invbeta, "x2s")
                                invbeta1sliceMp = self.extract_slice_data_at_t(invbeta1, "x2s")
                                invbeta2sliceMp = self.extract_slice_data_at_t(invbeta2, "x2s")
                                invbeta3sliceMp = self.extract_slice_data_at_t(invbeta3, "x2s")

                            f["invbeta"].create_dataset(str(t), data=invbetasliceVert)
                        if (force_load and qv in qs) or str(t) not in f["invbeta1"]:
                            try:
                                del f["invbeta1"][str(t)]
                            except: pass
                            f["invbeta1"].create_dataset(str(t), data=invbeta1sliceVert)
                        if (force_load and qv in qs) or str(t) not in f["invbeta2"]:
                            try:
                                del f["invbeta2"][str(t)]
                            except: pass
                            f["invbeta2"].create_dataset(str(t), data=invbeta2sliceVert)
                        if (force_load and qv in qs) or str(t) not in f["invbeta3"]:
                            try:
                                del f["invbeta3"][str(t)]
                            except: pass
                            f["invbeta3"].create_dataset(str(t), data=invbeta3sliceVert)

                        if (force_load and qv2 in qs) or str(t) not in f[qv2] or str(t) not in f[qv3]:
                            print("getting data for {} and {}, t {}".format(qv2, qv3, t))
                            if qsdata is None or "rho" not in qsdata.keys() or "vel3" not in qsdata.keys() or "Bcc3" not in qsdata.keys() or "Bcc2" not in qsdata.keys():
                                qsdata = self.load_q_data_raw_at_t(["rho", "vel3", "Bcc3", "Bcc2"], t)

                        if (force_load and qvrr in qs) or str(t) not in f["Rxrad"]:
                            try:
                                del f["Rxrad"][str(t)]
                                del f["Rxpol"][str(t)]
                                del f["Mxrad"][str(t)]
                                del f["Mxpol"][str(t)]
                            except: pass
                            # calculate difference from mean
                            dvphi = qsdata["vel3"] - np.mean(qsdata["vel3"], axis=0)
                            # calculate stresses
                            trphi = np.multiply(qsdata["rho"], np.multiply(dvphi, qsdata["vel1"]))
                            ttphi = np.multiply(qsdata["rho"], np.multiply(dvphi, qsdata["vel2"]))
                            mrphi = -np.multiply(qsdata["Bcc1"], qsdata["Bcc3"])
                            mtphi = -np.multiply(qsdata["Bcc2"], qsdata["Bcc3"])

                            trphiVert = self.extract_slice_data_at_t(trphi, "x3s")
                            ttphiVert = self.extract_slice_data_at_t(ttphi, "x3s")
                            mrphiVert = self.extract_slice_data_at_t(mrphi, "x3s")
                            mtphiVert = self.extract_slice_data_at_t(mtphi, "x3s")
                            trphiMp = self.extract_slice_data_at_t(trphi, "x2s")
                            ttphiMp = self.extract_slice_data_at_t(ttphi, "x2s")
                            mrphiMp = self.extract_slice_data_at_t(mrphi, "x2s")
                            mtphiMp = self.extract_slice_data_at_t(mtphi, "x2s")

                            f["Rxrad"].create_dataset(str(t), data=trphiVert)
                            f["Rxpol"].create_dataset(str(t), data=ttphiVert)
                            f["Mxrad"].create_dataset(str(t), data=mrphiVert)
                            f["Mxpol"].create_dataset(str(t), data=mtphiVert)


                            if (force_load and qv2 in qs) or str(t) not in f["Qtheta"]:
                                try:
                                    del f["Qtheta"][str(t)]
                                except: pass
                                thetadata = np.abs(2*np.pi*np.divide(qsdata["Bcc2"], np.multiply(qsdata["vel3"], np.sqrt(qsdata["rho"]))))
                                qthetaVert = self.extract_slice_data_at_t(thetadata, "x3s")
                                qthetaVert = np.divide(qthetaVert, self.dtheta_vert)
                                f["Qtheta"].create_dataset(str(t), data=qthetaVert)
                            if (force_load and qv2 in qs) or str(t) not in f["Qphi"]:
                                try:
                                    del f["Qphi"][str(t)]
                                except: pass
                                phidata = np.abs(2*np.pi*np.divide(qsdata["Bcc3"], np.multiply(qsdata["vel3"], np.sqrt(qsdata["rho"]))))
                                qphiVert = self.extract_slice_data_at_t(phidata, "x3s")
                                qphiVert = np.divide(qphiVert, self.dphi_vert)
                                f["Qphi"].create_dataset(str(t), data=qphiVert)
                            # print("Saved invbeta3 vert at t={}".format(t))

                if mpSlices: # NOTE: probably not right yet XX
                    with h5py.File(self.mp_slice_path, 'a') as f:
                        if "rhovr" not in f:
                            f.create_group("rhovr")
                        if "invbeta" not in f:
                            f.create_group("invbeta")
                        if "invbeta1" not in f:
                            f.create_group("invbeta1")
                        if "invbeta2" not in f:
                            f.create_group("invbeta2")
                        if "invbeta3" not in f:
                            f.create_group("invbeta3")
                        if "Rxrad" not in f:
                            f.create_group("Rxrad")
                        if "Rxpol" not in f:
                            f.create_group("Rxpol")
                        if "Mxrad" not in f:
                            f.create_group("Mxrad")
                        if "Mxpol" not in f:
                            f.create_group("Mxpol")
                        if "Qphi" not in f:
                            f.create_group("Qphi")
                        if "Qtheta" not in f:
                            f.create_group("Qtheta")

                        if force_load or str(t) not in f["rhovr"]:
                            try:
                                del f["rhovr"][str(t)]
                            except: pass
                            f["rhovr"].create_dataset(str(t), data=rhovrMp)
                            # print("Saved rhovr mp at t={}".format(t))
                        if force_load or str(t) not in f["invbeta"]:
                            try:
                                del f["invbeta"][str(t)]
                            except: pass
                            f["invbeta"].create_dataset(str(t), data=invbetasliceMp)
                            # print("Saved invbeta mp at t={}".format(t))
                        if force_load or str(t) not in f["invbeta1"]:
                            try:
                                del f["invbeta1"][str(t)]
                            except: pass
                            f["invbeta1"].create_dataset(str(t), data=invbeta1sliceMp)
                        if force_load or str(t) not in f["invbeta2"]:
                            try:
                                del f["invbeta2"][str(t)]
                            except: pass
                            f["invbeta2"].create_dataset(str(t), data=invbeta2sliceMp)
                        if force_load or str(t) not in f["invbeta3"]:
                            try:
                                del f["invbeta3"][str(t)]
                            except: pass
                            f["invbeta3"].create_dataset(str(t), data=invbeta3sliceMp)
                        if force_load or str(t) not in f["Rxrad"]:
                            try:
                                del f["Rxrad"][str(t)]
                            except: pass
                            f["Rxrad"].create_dataset(str(t), data=trphiMp)
                        if force_load or str(t) not in f["Rxpol"]:
                            try:
                                del f["Rxpol"][str(t)]
                            except: pass
                            f["Rxpol"].create_dataset(str(t), data=ttphiMp)
                        if force_load or str(t) not in f["Mxrad"]:
                            try:
                                del f["Mxrad"][str(t)]
                            except: pass
                            f["Mxrad"].create_dataset(str(t), data=mrphiMp)
                        if force_load or str(t) not in f["Mxpol"]:
                            try:
                                del f["Mxpol"][str(t)]
                            except: pass
                            f["Mxpol"].create_dataset(str(t), data=mtphiMp)
                        if force_load or str(t) not in f["Qtheta"]:
                            try:
                                del f["Qtheta"][str(t)]
                            except: pass
                            f["Qtheta"].create_dataset(str(t), data=qthetaMp)
                        if force_load or str(t) not in f["Qphi"]:
                            try:
                                del f["Qphi"][str(t)]
                            except: pass
                            f["Qphi"].create_dataset(str(t), data=qphiMp)

            if extremaLoad:
                print(extrema.loc[t].index)
                extrema = extrema.apply(pd.to_numeric)
                print(extrema.index)
                extrema = extrema.drop_duplicates(keep='first')
                print(extrema.loc[t].index)
                # if not everything is in place for this time step,
                # load the new q variables
                print(extrema.loc[t])
                condition = [q not in extrema.loc[t].index for q in qs]
                qs_to_load = np.extract(condition, qs)
                print(qs_to_load)
                new_extrema = self.extract_q_extrema_at_t(t, qs_to_load, qsdata)
                extrema.reset_index(inplace=True)
                new_extrema.reset_index(inplace=True)
                extrema = extrema.append(new_extrema, ignore_index=False)
                extrema["time"] = extrema["time"].astype('int')
                extrema.set_index(['time', 'variable'], inplace=True)

        # save to hdf5
        if extremaLoad:
            with open(self.extrema_path, 'a') as f:
                extrema.to_csv(f, header=write_header)

        # times = extrema.index.values.tolist()
            # print(times)




    # ---------------------------------------------------
    def load_q_slice_at_t(self, q, t, slicestr, sliceind=None, **kwargs):
        force_load = kwargs.get("force_load", False)
        avg = kwargs.get("avg", False)
        openfile = False

        # see if slice file exists
        if slicestr == "x2s":
            fpath = self.mp_slice_path
        elif avg:
            fpath = self.vert_slice_aAvg_path
        elif slicestr == "x3s":
            fpath = self.vert_slice_path
        else:
            print("ERROR: ")
            print("unknown slicing. Options: x2s or x3s")


        with h5py.File(fpath, 'r') as f:
            if q in f and str(t) in f[q] and not force_load:
                print("Loading slice from pre-sliced file")
                return f[q][str(t)][:]
            else:
                openfile = True
        if openfile:
            with h5py.File(fpath, 'a') as f:
                # print("Loading slice from raw data")
                # load raw data
                qdata = self.load_q_data_raw_at_t([q], t)[q]
                # slice data
                slicedata = self.extract_slice_data_at_t(qdata, slicestr, sliceind)
                # save to file
                try:
                    del f[q][str(t)]
                except:
                    pass
                finally:
                    if q not in f:
                        f.create_group(q)
                    f[q].create_dataset(str(t), data=slicedata)
                return slicedata
        return


    # ---------------------------------------------------
    def extract_slice_data_at_t(self, data, slicestr, sliceind=None, **kwargs):
        avg = kwargs.get("avg", False)
        if slicestr == "x2s": # i.e. midplane
            if sliceind is None:
                if self.nx2 % 2 == 0:
                    slice_data = np.mean(data[:, int(self.nx2 /2 - 1):int(self.nx2 /2+1), :], axis=1)
                else:
                    slice_data = data[:, int(self.nx2/2), :]
            else:
                slice_data = data[:, sliceind, :]
            # join through boundaries
            slice_data = np.vstack((slice_data[-1:, :], slice_data, slice_data[:1, :]))
        elif slicestr == "x3s": # i.e. vertical slice
            if sliceind is None:
                if avg:
                    vals_right = np.mean(data, axis=0)
                    vals_left = vals_right
                else:
                    vals_right = 0.5*(data[-1,:,:] + data[0,:,:])
                    vals_left = 0.5*(data[int((self.nx3/2)-1), :,:] + data[int(self.nx3/2), :, :])
                slice_data = np.vstack((vals_left[:1, :], vals_right,
                          vals_left[::-1, :], vals_right[:1, :]))
            else:
                slice_data = data[sliceind, :, :] # remember [phi, theta, r] not [r, theta, phi]
        else:
            print(slicestr + " is not a valid option!!")
            print("exiting...")
            return 0;

        return slice_data

    #------------------------------------------------------------------------
    def load_q_extrema_at_t(self, t, qs):
        """load all existing extrema for different types of data from extrema (reduced) file
        """
        extrema = pd.read_csv(self.extrema_path, index_col=["time", "variable"])
        # for col in extrema.columns:
            # print(extrema[col].dtypes)

        if t in extrema.index.get_level_values("time"):
            # check if variables have been loaded in
            condition = [q in extrema.loc[t].index for q in qs]
            qs_in = np.extract(condition, qs)

            # print(extrema.loc[60])
            # print(extrema.loc[(t, qs[0])])
            tuple_list = [(t, q) for q in qs_in]
            return extrema.loc[tuple_list] #also need by q
        else:
            return pd.DataFrame()

    #------------------------------------------------------------------------
    def add_q_extrema_at_t(self, t, qs=None, **kwargs):
        """add extrema for different types of data to file.

        Load extrema from file, calculate new extrema if necessary and write new extrema to csv file.
        XX To do: extend to add different types of extrema without re-calculating all of others.
        """
        self.header = ["time", "variable", "vol_min", "vol_max", "x2s_min", "x2s_max", "x3s_min", "x3s_max"]
        write_header = False
        extrema = pd.DataFrame()

        force_load = kwargs.get("force_load", False)

        if qs is None:
            qs = self.variables
        # if file doesn't exist, export header. if it does, try to load data
        if not os.path.isfile(self.extrema_path):
            write_header = True
        else:
            extrema = self.load_q_extrema_at_t(t, qs)
            # for col in extrema.columns:
                # print(extrema[col].dtypes)

        if extrema.empty or extrema.shape[0] < len(qs) or force_load:
            # will need to load some in
            if extrema.empty or force_load:
                qs_to_load = qs
            else:
                condition = [q not in extrema.loc[t].index for q in qs]
                qs_to_load = np.extract(condition, qs)
            data = self.load_q_data_raw_at_t(qs_to_load, t)
            if data is None:
                print("ERROR:")
                print("DATA DOES NOT EXIST. ")
                print("EXITING...")
                return None

            new_extrema = self.extract_q_extrema_at_t(t, qs_to_load, data)
            # print(new_extrema.dtypes)
            extrema = pd.concat([extrema, new_extrema])
            # for col in extrema.columns:
                # print(extrema[col].dtypes)

            # write to file (only new extrema since append mode)
            with open(self.extrema_path, 'a') as f:
                # print("writing extrema to file")
                new_extrema.to_csv(f, header=write_header)

        return extrema

    #------------------------------------------------------------------------
    def find_all_t(self, to_load=False):
        """Return list of all time steps found in data load path"""

        # get list of times from present raw data files
        if to_load or not os.path.isfile(self.consts_path):
            # Get max time.
            import fnmatch
            files = fnmatch.filter(os.listdir(self.data_load_fp), '*out1*.athdf')
            self.times = None

            for file in files:
                with h5py.File(self.data_load_fp + file, 'r') as f:
                    # extract code time from file
                    tval = f.attrs["Time"]
                    tind = int(file.split('.')[-2])
                    if self.times is None:
                        self.times = np.array([tind, tval])
                    else:
                        self.times = np.vstack((self.times, [tind, tval]))
            self.times = np.sort(self.times, axis=0)
            self.tinds = self.times[:, 0].astype(int)
            self.tvals = self.times[:, 1]

            # save to file
            with h5py.File(self.consts_path, 'a') as f:
                try:
                    # del f["Time"]
                    # del f["time"]
                    del f["Times"]
                    del f["Time Values"]
                    del f["Time Indices"]
                except: pass
                finally:
                    f.create_dataset("Times", data = self.times)
                    f.create_dataset("Time Values", data = self.tvals)
                    f.create_dataset("Time Indices", data = self.tinds)

        # load from time files
        else:
            with h5py.File(self.consts_path, 'r') as f:
                self.times = f["Times"][:]
                self.tinds = f["Time Indices"][:]
                self.tvals = f["Time Values"][:]

        return self.times

    #------------------------------------------------------------------------
    def get_q_extrema_trange(self, ts=None, qs=None, force_load=False):
        """Get extrema of qs quantities over time range.
        """
        print("Finding extrema...")
        extrema = pd.DataFrame()

        if qs is None:
            qs = self.variables

        if ts is None:
            ts = self.find_all_t()

        # load all available files
        for t in ts:
            new_extrema = self.add_q_extrema_at_t(t, qs, force_load=force_load)
            extrema = pd.concat([extrema, new_extrema])
        print("Done.")
        return extrema

    #------------------------------------------------------------------------
    def extract_qt_extrema(self, **kwargs):
        """ Return maximum value of q over time range ts

        options for idstr:
        - "vol": volume-averaged
        - "x1s": slice middle of x1
        - "x2s": slice middle of x2
        - "x3s": slice middle of x3

        TO DO XXX:
        - make safe when ts doesn't have raw data
        """

        qs = kwargs.get("qs", self.variables)
        ts = kwargs.get("ts", self.find_all_t())
        to_save = kwargs.get("to_save", True)
        force_load = kwargs.get("force_load", False)
        idstr = kwargs.get("idstr", "vol")
        if idstr not in self.known_eidstr:
            print("ERROR")
            print(idstr + " NOT KNOWN")
            print("EXITING")
            return None

        ext = {}
        extrema = self.get_q_extrema_trange(ts, qs, force_load)

        # discard unneeded columns
        for col in extrema.columns:
            if idstr not in col:
                extrema.drop(columns=[col], inplace=True)

        # find max/min for each variable
        for q in qs:
            # Note: need variables to be in level 1!!
            slice = extrema.xs(q, level=1)
            minval = slice[idstr+"_min"].min()
            maxval = slice[idstr+"_max"].max()
            ext[q] = [minval, maxval]

        time_extrema = pd.DataFrame(ext)

        # save in file
        if to_save:
            time_extrema.to_hdf(self.consts_path, idstr, mode='a')

        return time_extrema

    #------------------------------------------------------------------------
    def extract_q_shell_data_at_t(self, q, t, center=0, nbins=None, overwrite=False):
        # separate into quadrants
        # NOTE: this takes planet y and z position to be zero!!
        data = self.load_q_data_raw_at_t([q], t)[q]

        x1p = self.x1v[self.x1v - center > 0]
        x1m = self.x1v[self.x1v - center < 0]

        # cut off points that have (radius) |x - xp| > |L - xp|,
        # which would favor the x < xp values for xp>0 and vice versa
        max_x = np.min([np.max(np.abs(x1p - center)), np.max(np.abs(x1m - center))])
        max_y = np.max(self.x2v)
        max_z = np.max(self.x3v)
        max_r = np.min([max_x, max_y, max_z])

        x_cond = np.abs(self.x1v - center) <= max_r
        y_cond = np.abs(self.x2v) <= max_r
        z_cond = np.abs(self.x3v) <= max_r

        x1 = self.x1v[x_cond]
        x2 = self.x2v[y_cond]
        x3 = self.x3v[z_cond]

        # make grid
        xg, yg, zg = np.meshgrid(x1, x2, x3, indexing='ij')
        dist = np.sqrt((xg - center)**2 + yg**2 + zg**2)

        planetind = (np.abs(self.x1v - self.xp)).argmin()
        nbins_max = self.x1v.size - planetind
        nbins_min = int(1 + 3.322*np.log10(np.unique(dist).size)) # sturge's rule

        if nbins is None:
            nbins = int(np.mean([nbins_max, nbins_min]))
        elif nbins > nbins_max:
            print("CAUTION: there are more bins than cells in the x-direction.")
        elif nbins < nbins_min:
            print("CAUTION: there are fewer bins than recommended.")
        sphere_data = data[z_cond, :, :]
        sphere_data = sphere_data[:, y_cond, :]
        sphere_data = sphere_data[:, :, x_cond]

        # bin ranges: first bin will contain values of cells that have
        # radius between bin[0] and bin[1]; bin[1] contains values of
        # cells that have radius between bin[1] and bin[2], etc.
        binfirst = np.min(dist)
        binlast = max_r # all cells with dist between max_r and np.max(dist) will be in last bin, which should be thrown away (corners of cube in which sphere is inscribed)
        bins = np.linspace(binfirst, binlast, nbins)
        q_binned = np.zeros(bins.shape)
        bin_ncells = np.zeros(bins.shape)

        # bin the data
        for i in np.arange(dist.shape[0]):
            for j in np.arange(dist.shape[1]):
                for k in np.arange(dist.shape[2]):
                    # find bin that q should be added to
                    d = dist[i,j,k]
                    binind = np.where(d >= bins)[0][-1] # put into lower bin
                    q_binned[binind] += sphere_data[k, j, i] # remember the backwards indexing!
                    # print(data[k, j, i])
                    bin_ncells[binind] += 1
        return [bins[:-1], bin_ncells[:-1], q_binned[:-1]];

    #------------------------------------------------------------------------

    def extract_q_yashell_data_at_t(self, q, t, center=0, nbins=None, overwrite=False):
        data = self.load_q_data_raw_at_t([q], t)[q]

        # separate into quadrants
        x1p = self.x1v[self.x1v - center > 0]
        x1m = self.x1v[self.x1v - center < 0]

        # cut off points that have (radius) |x - xp| > |L - xp|,
        # which would favor the x < xp values for xp>0 and vice versa
        max_x = np.min([np.max(np.abs(x1p - center)), np.max(np.abs(x1m - center))])
        max_y = np.max(self.x2v)
        max_z = np.max(self.x3v)
        max_r = np.min([max_x, max_y, max_z])

        x_cond = np.abs(self.x1v - center) <= max_r
        y_cond = np.abs(self.x2v) <= max_r
        z_cond = np.abs(self.x3v) <= max_r

        y_condp = y_cond & (self.x2v > 0)
        y_condm = y_cond & (self.x2v < 0)

        x1 = self.x1v[x_cond]
        x2p = self.x2v[y_condp]
        x2m = self.x2v[y_condm]
        x3 = self.x3v[z_cond]


        # make grid
        xgp, ygp, zgp = np.meshgrid(x1, x2p, x3, indexing='ij')
        xgm, ygm, zgm = np.meshgrid(x1, x2m, x3, indexing='ij')
        distp = np.sqrt((xgp - center)**2 + ygp**2 + zgp**2)
        distm = np.sqrt((xgm - center)**2 + ygm**2 + zgm**2)

        planetind = (np.abs(self.x1v - self.xp)).argmin()
        nbins_max = self.x1v.size - planetind
        nbins_min = int(1 + 3.322*np.log10(np.unique(np.append(distm, distp)).size)) # sturge's rule

        if nbins is None:
            nbins = int(np.mean([nbins_max, nbins_min]))
        elif nbins > nbins_max:
            print("CAUTION: there are more bins than cells in the x-direction.")
        elif nbins < nbins_min:
            print("CAUTION: there are fewer bins than recommended.")

        sphere_datam = data[z_cond, :, :]
        sphere_datam = sphere_datam[:, y_condm, :]
        sphere_datam = sphere_datam[:, :, x_cond]
        sphere_datap = data[z_cond, :, :]
        sphere_datap = sphere_datap[:, y_condp, :]
        sphere_datap = sphere_datap[:, :, x_cond]

        # bin ranges: first bin will contain values of cells that have
        # radius between bin[0] and bin[1]; bin[1] contains values of
        # cells that have radius between bin[1] and bin[2], etc.
        binfirst = np.min([np.min(distp), np.min(distm)])
        binlast = max_r # all cells with dist between max_r and np.max(dist) will be in last bin, which should be thrown away (corners of cube in which sphere is inscribed)
        bins = np.linspace(binfirst, binlast, nbins)
        q_binnedp = np.zeros(bins.shape); q_binnedm = np.zeros(bins.shape)
        bin_ncellsp = np.zeros(bins.shape); bin_ncellsm = np.zeros(bins.shape)

        # bin the data
        [nzp, nyp, nxp] = sphere_datap.shape
        [nzm, nym, nxm] = sphere_datam.shape
        for i in np.arange(nxp):
            for j in np.arange(nyp):
                for k in np.arange(nzp):
                    # find bin that q should be added to
                    d = distp[i,j,k]
                    binind = np.where(d >= bins)[0][-1] # put into lower bin
                    q_binnedp[binind] += sphere_datap[k, j, i] # remember the backwards indexing!
                    # print(data[k, j, i])
                    bin_ncellsp[binind] += 1
        for i in np.arange(nxm):
            for j in np.arange(nym):
                for k in np.arange(nzm):
                    # find bin that q should be added to
                    d = distm[i,j,k]
                    binind = np.where(d >= bins)[0][-1] # put into lower bin
                    q_binnedm[binind] += sphere_datam[k, j, i] # remember the backwards indexing!
                    # print(data[k, j, i])
                    bin_ncellsm[binind] += 1
        return [bins[:-1], bin_ncellsp[:-1], bin_ncellsm[:-1], q_binnedp[:-1], q_binnedm[:-1]];
    #------------------------------------------------------------------------

    def extract_yagforce_data_at_t(self, t, center=0, nbins=None, overwrite=False, **kwargs):
        data = self.load_q_data_raw_at_t(["rho"], t)
        timesim = data["Time"]
        data = data["rho"]
        planetind = (np.abs(self.x1v - self.xp)).argmin()
        # print("planet ind: " + str(planetind))

        # find largest distance to planet in x-direction
        x1p = np.abs(np.max(self.x1v) - center)
        x1m = np.abs(np.min(self.x1v) - center)

        # cut off points that have (radius) |x - xp| > |L - xp|,
        # which would favor the x < xp values for xp>0 and vice versa
        max_x = np.min([x1p, x1m])
        max_y = np.max(self.x2v)
        max_z = np.max(self.x3v)
        max_r = np.min([max_x, max_y, max_z])

        x_cond = np.abs(self.x1v - center) <= max_r
        y_cond = np.abs(self.x2v) <= max_r
        z_cond = np.abs(self.x3v) <= max_r

        y_condp = y_cond & (self.x2v > 0)
        y_condm = y_cond & (self.x2v < 0)

        x1 = self.x1v[x_cond]
        x2p = self.x2v[y_condp]
        x2m = self.x2v[y_condm]
        x3 = self.x3v[z_cond]

        # make grid
        xgp, ygp, zgp = np.meshgrid(x1, x2p, x3, indexing='ij')
        xgm, ygm, zgm = np.meshgrid(x1, x2m, x3, indexing='ij')
        distp = np.sqrt((xgp - center)**2 + ygp**2 + zgp**2)
        distm = np.sqrt((xgm - center)**2 + ygm**2 + zgm**2)

        planetind = (np.abs(self.x1v - self.xp)).argmin()
        nbins_max = self.x1v.size - planetind
        nbins_min = int(1 + 3.322*np.log10(np.unique(np.append(distm, distp)).size)) # sturge's rule

        if nbins is None:
            nbins = int(np.mean([nbins_max, nbins_min]))
        elif nbins > nbins_max:
            print("CAUTION: there are more bins than cells in the x-direction.")
        elif nbins < nbins_min:
            print("CAUTION: there are fewer bins than recommended.")

        sphere_datam = data[z_cond, :, :]
        sphere_datam = sphere_datam[:, y_condm, :]
        sphere_datam = sphere_datam[:, :, x_cond]
        sphere_datap = data[z_cond, :, :]
        sphere_datap = sphere_datap[:, y_condp, :]
        sphere_datap = sphere_datap[:, :, x_cond]

        # bin ranges: first bin will contain values of cells that have
        # radius between bin[0] and bin[1]; bin[1] contains values of
        # cells that have radius between bin[1] and bin[2], etc.
        binfirst = np.min([np.min(distp), np.min(distm)])
        binlast = max_r # all cells with dist between max_r and np.max(dist) will be in last bin, which should be thrown away (corners of cube in which sphere is inscribed)
        bins = np.linspace(binfirst, binlast, nbins)
        forcep = np.zeros(bins.shape); forcem = np.zeros(bins.shape)
        bin_ncellsp = np.zeros(bins.shape); bin_ncellsm = np.zeros(bins.shape)

        # bin the data
        [nzp, nyp, nxp] = sphere_datap.shape
        [nzm, nym, nxm] = sphere_datam.shape

        for i in np.arange(nxp):
            for j in np.arange(nyp):
                for k in np.arange(nzp):
                    # find bin that q should be added to
                    d = distp[i,j,k]
                    binind = np.where(d >= bins)[0][-1] # put into lower bin
                    cell_mass = self.cell_vol * sphere_datap[k, j, i]

                    force_mag = cell_mass/d**2
                    sintheta = ygp[i, j, k]/d
                    forcep[binind] += force_mag*sintheta
                    bin_ncellsp[binind] += 1
        for i in np.arange(nxm):
            for j in np.arange(nym):
                for k in np.arange(nzm):
                    # find bin that q should be added to
                    d = distm[i,j,k]
                    binind = np.where(d >= bins)[0][-1] # put into lower bin
                    cell_mass = self.cell_vol * sphere_datam[k, j, i]
                    force_mag = cell_mass/d**2
                    rad2 = zgm[i, j, k]**2 + (xgm[i, j, k] - center)**2
                    sintheta = ygm[i, j, k]/np.sqrt(rad2 + ygm[i, j, k]**2)
                    forcem[binind] += force_mag*sintheta
                    bin_ncellsm[binind] += 1

        # save to file
        header = "Index, time, y<0 sum, y>0 sum, radial sum"
        # try loading from file:
        if os.path.isfile(self.ygforce_path):
            ygforce = np.loadtxt(self.ygforce_path, skiprows=1, delimiter=",")
            if ygforce.ndim == 1:
                tin = (ygforce[0] == t)
            else:
                ts = ygforce[:, 0]
                tin = (ts == t).any()
            if not tin:
                newrow = np.array([t, timesim, np.sum(forcem), np.sum(forcep), np.sum(forcem+forcep)])
                ygforce = np.vstack((ygforce, newrow))
                np.savetxt(self.ygforce_path, ygforce, header=header, delimiter=",")
        else:
            newrow = np.transpose(np.array([t, timesim, np.sum(forcem), np.sum(forcep), np.sum(forcem+forcep)]))
            np.savetxt(self.ygforce_path, newrow[None], header=header, delimiter=",")

        return [bins[:-1], bin_ncellsp[:-1], bin_ncellsm[:-1], forcep[:-1], forcem[:-1]];

    # ------------------------------------------------------------------------
    # -------- BROKEN BELOW -------------------------------------------------
    # ------------------------------------------------------------------------
#
#
#
    # def get_q_data_slice_at_t(self, q, t, **kwargs):
        # """Slice data and save as dictionary so can access specific sections only
#
        # Keyword arguments:
        # q: the variable to be queried. Should be in self.variables.
        # t: the timestep (i.e. number in filename, not index) to be queried.
        # kwargs:
           # midplane: boolean. If true, slice at midplane. Otherwise, vertical slice (default).
#
        # Returns: requested data.
        # """
#
        # to_write = True
        # ts = str(t)
#
        # data_slice_fname = "data_slice_"
        # grid_fname = self.data_save_fp + "grid_"
#
        # if 'midplane' in kwargs and kwargs['midplane']:
            # data_slice_fname += "midplane.hdf5"
            # grid_fname += "midplane.p"
        # else:
            # data_slice_fname += "vertical.hdf5"
            # grid_fname += "vertical.p"
            # kwargs['midplane'] = False
        # full_fname = self.data_save_fp + data_slice_fname
#
        # if os.path.isfile(full_fname):
            # with h5py.File(full_fname, 'r') as f:
                # if ts in f.keys() and q in f[ts].keys():
                    # print("loading {}, {} from .hdf5 file".format(ts, q))
                    # data = np.array(f.get(ts+'/'+q))
                    # to_write = False
#
        # if to_write:
            # if q in self.variables:
                # data = self.load_q_data_slice_at_t(q, t, **kwargs)
                # # XX safeguard for data not existing!!!
            # elif q in self.calcvars:
                # data = self.calculate_q_slice_at_t(q, t, **kwargs)
                # to_write = False
            # else:
                # print(q + " is an invalid variable.")
#
        # if to_write:
            # with h5py.File(full_fname, 'a') as f:
                # f.create_dataset(ts+'/'+q, data=data)
#
        # if not os.path.isfile(grid_fname):
            # r = self.r
            # theta = self.theta
            # phi = self.phi
            # nx3 = len(phi)
#
            # if kwargs['midplane']:
                # # phi doesn't go all the way to 2pi, so extend it
                # phi_extended = \
                    # np.concatenate((phi[-1:] - 2.0 * np.pi, phi, phi[:1] + 2.0 * np.pi))
                # phi_extended -= 0.5 * 2.0 * np.pi / nx3
                # r_grid, phi_grid = np.meshgrid(r, phi_extended)
                # x_grid = r_grid * np.cos(phi_grid)
                # y_grid = r_grid * np.sin(phi_grid)
            # else:
                # theta_extended = np.concatenate((-theta[0:1], theta, 2.0 * np.pi - theta[::-1],
                                                # 2.0 * np.pi + theta[0:1]))
                # theta_extended_corrected = theta_extended - 0.5 * (theta[1] - theta[0])
                # r_grid, theta_grid = np.meshgrid(r, theta_extended_corrected)
                # x_grid = r_grid * np.sin(theta_grid)
                # y_grid = r_grid * np.cos(theta_grid)
#
            # # ------- save scalar grid ----
            # grid = {"x_grid": x_grid, "y_grid": y_grid}
            # pickle.dump(grid, open(grid_fname, "wb"))
#
        # return data
#
#
#
    # def calculate_q_slice_at_t(self, q, t, **kwargs):
        # """Slice data and save as dictionary so can access specific sections only
#
        # Keyword arguments:
        # q: the variable to be queried. Should be either "beta" or "rhovr"
        # t: the timestep (i.e. number in filename, not index) to be queried.
        # kwargs:
           # midplane: boolean. If true, slice at midplane. Otherwise, vertical slice (default).
#
        # Returns: requested data.
        # """
        # to_write = True
        # ts = str(t)
#
        # if 'midplane' in kwargs and kwargs['midplane']:
            # data_slice_fname = "data_slice_midplane.hdf5"
        # else:
            # data_slice_fname = "data_slice_vertical.hdf5"
            # kwargs['midplane'] = False
        # full_fname = self.data_save_fp + data_slice_fname
#
#
        # if os.path.isfile(full_fname):
            # with h5py.File(full_fname, 'r') as f:
                # if ts in f.keys() and q in f[ts].keys():
                    # #print("loading {}, {} from .hdf5 file".format(ts, q))
                    # data = np.array(f.get(ts+'/'+q))
                    # to_write = False
        # if to_write:
            # if q == "beta":
                # press = self.get_q_data_slice_at_t("press", t, midplane=kwargs['midplane'])
                # bcc1 = self.get_q_data_slice_at_t("Bcc1", t, midplane=kwargs['midplane'])
                # bcc2 = self.get_q_data_slice_at_t("Bcc2", t, midplane=kwargs['midplane'])
                # bcc3 = self.get_q_data_slice_at_t("Bcc3", t, midplane=kwargs['midplane'])
#
                # b2 = bcc1**2 + bcc2**2 + bcc3**2
                # b2min = 10*sys.float_info.epsilon
                # b2[b2 <= b2min] = b2min
#
                # data = 2*press/b2
                # fmax = 10**9
                # data[data > fmax] = fmax
            # elif q == "rhovr":
                # rho = self.get_q_data_slice_at_t("rho", t, midplane=kwargs['midplane'])
                # vr = self.get_q_data_slice_at_t("vel1", t, midplane=kwargs['midplane'])
                # data = rho*vr
            # else:
                # print("error: " + q + " is not a valid option. Options are 'beta' and 'rhovr' ")
                # return 0
            # with h5py.File(full_fname, 'a') as f:
                # f.create_dataset(ts+'/'+q, data=data)
        # return data
#
#
    # #------------------------------------------------------------------------
#
    # def plot_q_radial_profile_at_t(self, q, t, to_avg=True, ri_range=[0, -1], fig_num=1, to_save=False, overwrite=False, theta=None, phi=0):
        # """Plot variable q at specific time, location.
#
        # Can specify theta at which data will be taken.
        # If to_avg, will average over azimuth (phi). Otherwise, input theta and phi (defaults are midplane for theta and 0 for phi)
#
        # Keyword Arguments:
        # q: the variable to be queried. Should be in self.variables.
        # t: the timestep (i.e. number in filename, not index) to be queried.
        # to_avg: if true, average over phi.
        # ri_range: Index range of radii over which plot should be made. Default is full range.
        # fig_num: figure on which to plot. Default is 1.
        # to_save: if true, save the figure.
        # overwrite: if true, overwrite files with the same name. If false, create file with "vx", where x is high enough to not overwrite any existing files.
        # theta: index, not actual value. Default is midplane
        # phi: index, not actual value. Default is 0. Will only be used if to_avg is false.
#
        # Improvements:
        # -plot multiple variables (qs) on same plot.
        # """
        # lab = q
        # if to_avg:
            # lab += " averaged over phi at "
            # if q not in self.data_avgPhi or t not in self.data_avgPhi[q] or theta not in self.data_avgPhi[q][t]:
                # self.get_q_data_avgPhi_at_t(q, t, theta)
            # if theta is None:
                # theta = int(self.theta.size/2)
                # lab += "midplane"
            # else:
                # lab += "theta index {:.2f}".format(theta)
            # data = self.data_avgPhi[q][t][theta]
        # else:
            # if q not in self.data_giveThetaPhi or t not in self.data_giveThetaPhi[q]:
                # self.get_q_data_giveThetaPhi_at_t(q, t, theta, phi)
            # if theta is None:
                # theta = int(self.theta.size/2)
                # lab += " at midplane, phi = {:.2f}".format(self.phi[phi])
            # else:
                # lab += " (theta, phi) = ({:.2f},{:.2f})".format(self.theta[theta], self.phi[phi])
            # data = self.data_giveThetaPhi[q][t][theta]
        # ri_min = ri_range[0]
        # ri_max = ri_range[-1]
        # ymax = np.max(data[ri_min:ri_max])
        # ymin = np.min(data[ri_min:ri_max])
        # yrange = ymax-ymin
#
        # if yrange <= 0:
            # print("Data is identically zero. Exiting.")
            # yrange = sys.float_info.epsilon
            # return 0
#
        # plt.figure(fig_num)
        # plt.plot(self.r[ri_min:ri_max], data[ri_min:ri_max], label=lab, marker='o')
        # plt.xlabel("r")
        # plt.yscale("symlog",linthreshy=yrange/50)
        # #plt.ylim([-.25, .25])
        # plt.gca().tick_params(direction='in', right=True, top=True)
        # plt.legend()
        # plt.title("t = {}".format(t))
        # plt.tight_layout()
#
        # descript = "radial-profile"
        # pdir = self.fig_save_fp + descript + "/"+self.config+"_"+self.specs+"/"
        # if to_avg:
            # pdir += "avg-phi/plots/"
            # descript = descript+"_avg-phi"
        # else:
            # pdir += "give-phi/plots/"
            # if theta is None:
                # descript = descript+"_phi{}mp".format(phi)
            # else:
                # descript = descript+"_phi{}theta{}".format(phi, theta)
        # tstr = "_t{}".format(t)
        # rstr = "_r{:.1f}-{:.1f}".format(self.r[ri_min], self.r[ri_max])
        # pname = self.create_fname("plot", q, descript+tstr+rstr, ".png")
        # if to_save:
            # pvars = self.figsave(pdir, pname, overwrite)
        # else:
            # pvars = [pdir, pname]
        # return pvars
#
#
    # #------------------------------------------------------------------------
    # def movie_q_radial_profile(self, q, t_range, to_avg=True, ri_range=[0, -1], fig_num=1, to_save=False, overwrite=False, theta=None, phi=0):
        # """Make movie of variable q evolution over specific time, location.
#
        # Can specify theta at which data will be taken.
        # If to_avg, will average over azimuth (phi). Otherwise, input theta and phi (defaults are midplane for theta and 0 for phi)
#
        # Keyword Arguments:
        # q: the variable to be queried. Should be in self.variables.
        # t: the timesteps (i.e. number in filename, not index) to be included in movie.
        # to_avg: if true, average over phi.
        # ri_range: Index range of radii over which plot should be made. Default is full range.
        # fig_num: figure on which to plot. Default is 1.
        # to_save: if true, save the figure.
        # overwrite: if true, overwrite files with the same name. If false, create file with "vx", where x is high enough to not overwrite any existing files.
        # theta: index, not actual value. Default is midplane
        # phi: index, not actual value. Default is 0. Will only be used if to_avg is false.
#
        # Improvements:
        # -plot multiple variables (qs) on same plot.
        # """
        # lab = q
        # mp = False
#
        # # --- Get max/min values ----
        # self.get_qs_extrema([q])
        # if theta is None:
            # theta = int(self.theta.size/2)
        # if theta == int(self.theta.size/2):
            # mp = True
#
        # if to_avg:
            # [minval, maxval] = self.extrema_avg[q][theta]
        # else:
            # [minval, maxval] = self.extrema_navg[q][theta][phi]
#
        # # ---- Set labels -------
        # if to_avg:
            # lab += " averaged over phi at "
            # if q not in self.data_avgPhi:
                # self.data_avgPhi[q] = {}
            # if mp:
                # lab += "midplane"
            # else:
                # lab += "theta = {:.2f}".format(self.theta[theta])
        # else:
            # if q not in self.data_giveThetaPhi:
                # self.data_giveThetaPhi[q] = {}
            # if mp:
                # lab += " at midplane, phi = {:.2f}".format(self.phi[phi])
            # else:
                # lab += " (theta, phi) = ({:.2f},{:.2f})".format(self.theta[theta], self.phi[phi])
#
        # # ---- Plot ranges/labels ----
        # ri_min = ri_range[0]
        # ri_max = ri_range[-1]
        # r_range = self.r[ri_min:ri_max]
        # fig = plt.figure()
        # yrange = maxval - minval
        # #params = [ylim, [r_range[0], r_range[-1]], "r", q]
        # l, = plt.plot([], [], marker="o")
#
        # plt.xlim(r_range[0], r_range[-1])
        # plt.ylim(minval-np.abs(minval)*.1, maxval+np.abs(maxval)*.1)
        # plt.xlabel("r")
        # plt.legend([lab], loc="upper right")
        # plt.title("Time: ")
        # plt.yscale("symlog", linthreshy=yrange/50)
        # plt.gca().tick_params(direction='in', right=True, top=True)
        # plt.tight_layout()
#
        # if to_save:
            # print("making movie...")
            # # ------- Filename ----
            # tstr = "_t{}-{}".format(t_range[0], t_range[-1])
            # rstr = "_r{:.1f}-{:.1f}".format(r_range[0], r_range[-1])
            # descript = "radial-profile"
            # pdir = self.fig_save_fp + descript + "/"+self.config+"_"+self.specs+"/"
            # if to_avg:
                # pdir += "avg-phi/movies/"
                # if mp:
                    # descript += "_avg-phi-mp"
                # else:
                    # descript += "_avg-phi-theta{}".format(theta)
            # else:
                # pdir += "give-phi/movies/"
                # if mp:
                    # descript += "_phi{}mp".format(phi)
                # else:
                    # descript += "_phi{}theta{}".format(phi, theta)
            # if not os.path.exists(pdir):
                # os.makedirs(pdir)
            # fname = self.create_fname("movie", q, descript+tstr+rstr, "_v0")
            # if not overwrite:
                # count = 0
                # while os.path.isfile(pdir+fname+".mp4"):
                    # count += 1
                    # fname = self.rreplace(fname, count)
                # print("new name to avoid overwrite: "+ fname)
                # full_descript = descript+"_"+self.config+"_"+self.specs
            # metadata = dict(title=full_descript, artist='LH', comment=full_descript)
            # writer = FFMpegWriter(fps=5, metadata=metadata)
#
            # # ---- Actually make the movie! --------
            # with writer.saving(fig, pdir+fname+".mp4", 200):
                # for t in t_range:
                    # if to_avg:
                        # t_data = self.get_q_data_avgPhi_at_t(q, t, theta)
                        # self.data_avgPhi[q][t][theta] = {}
                        # self.data_raw = {}
                    # else:
                        # t_data = self.get_q_data_giveThetaPhi_at_t(q, t, theta, phi)
                        # self.data_giveThetaPhi[q][t][theta] = {}
                        # self.data_raw = {}
#
                    # l.set_data(r_range, t_data[ri_min:ri_max])
                    # plt.title("Time: {}".format(t))
                    # writer.grab_frame()
            # print("movie saved")
#
    # #------------------------------------------------------------------------
#
    # def plot_q_slice_at_t(self, q, t, **kwargs):
        # """Plot either vertical or midplane slice of variable q at time t.
#
        # Keyword Arguments:
        # q: the variable to be queried. Should be in self.variables.
        # t: the timestep (i.e. number in filename, not index) to be queried.
        # to_avg: if true, average over phi.
        # ri_range: Index range of radii over which plot should be made. Default is full range.
        # fig_num: figure on which to plot. Default is 1.: figure on which to plot. Default is 1.
        # to_save: if true, save the figure.
        # overwrite: if true, overwrite files with the same name. If false, create file with "vx", where x is high enough to not overwrite any existing files.
#
        # """
#
        # # Deal with kwargs
        # if 'ri_max' not in kwargs:
            # ri_max = -1
        # else:
            # ri_max = kwargs['ri_max']
        # if 'plot_lims' in kwargs:
            # plot_lims = kwargs['plot_lims']
#
        # if 'fig_num' not in kwargs:
            # fig_num = 1
        # else:
            # fig_num = kwargs['fig_num']
        # if 'to_avg' not in kwargs:
            # to_avg = False
        # else:
            # to_avg = kwargs['to_avg']
        # if 'overwrite' not in kwargs:
            # overwrite = False
        # else:
            # overwrite = kwargs['overwrite']
        # if 'to_save' not in kwargs:
            # to_save = False
        # else:
            # to_save = kwargs['to_save']
        # if 'logc' not in kwargs:
            # kwargs['logc'] = True
        # if 'vmin' in kwargs:
            # vmin = kwargs['vmin']
        # else: vmin = None
        # if 'vmax' in kwargs:
            # vmax = kwargs['vmax']
        # else: vmax = None
#
        # # --- get data, set up ------
#
        # if 'midplane' not in kwargs or not kwargs['midplane']:
            # tit_str = "Vertical slice of "
            # kwargs['midplane'] = False
        # else:
            # tit_str = "Midplane slice of "
        # tit_str += r"{}".format(self.varnames[q]) + " at t = {}".format(t)
#
        # grid_fname = self.data_save_fp + "grid_"
        # if kwargs['midplane']:
            # grid_fname += "midplane.p"
        # else:
            # grid_fname += "vertical.p"
#
        # vals = self.get_q_data_slice_at_t(q, t, midplane=kwargs['midplane'])
        # grid = pickle.load(open(grid_fname, "rb"))
        # x_grid = grid["x_grid"]
        # y_grid = grid["y_grid"]
#
        # if self.r is None:
            # self.r = pickle.load(open(self.data_save_fp + "r_steps.p", "rb"))
        # r_max = self.r[ri_max]
        # ri_min = 0
#
        # # ---- clean data ----
        # # Remove inf
        # fmax = 10**9
        # vals[vals > fmax] = fmax
        # # ---- Determine colormapping properties -----
        # ymin = np.min(vals)
        # ymax = np.max(vals)
        # yrange = ymax-ymin
#
        # if yrange <= 0:
            # print("Data is identically zero. Exiting")
            # yrange = sys.float_info.epsilon
            # return 0
#
        # if ymin < 0:
            # if 'cmap' not in kwargs:
                # kwargs['cmap'] = 'RdBu'
            # # Make cmap symmetric about zero
            # if vmin is None:
                # vmin = -np.max([np.abs(ymin), np.abs(ymax)])
            # if vmax is None:
                # vmax = -vmin
        # else:
            # if q == "beta":
                # kwargs['cmap'] = 'Blues_r'
            # else:
                # kwargs['cmap'] = 'Blues'
            # if vmin is None:
                # vmin = ymin
            # if vmax is None:
                # vmax = ymax
#
        # cmap = plt.get_cmap(kwargs['cmap'])
#
        # if kwargs['logc']:
            # if ymin <= 0:
                # norm = colors.SymLogNorm(linthresh=yrange/50.)
            # else:
                # norm = colors.LogNorm()
        # else:
            # norm = colors.Normalize()
#
#
        # # ------- Figure stuff------------
        # plt.figure(fig_num)
        # im = plt.pcolormesh(x_grid, y_grid, vals, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
        # if 'plot_lims' not in kwargs:
            # plot_lims = [[-r_max, r_max], [-r_max, r_max]]
#
        # plt.xlim(plot_lims[0])
        # plt.ylim(plot_lims[1])
#
        # plt.gca().set_aspect('equal')
        # #plt.xlim((-r_max, r_max))
        # #plt.ylim((-r_max, r_max))
        # #plt.xlim((0, r_max))
        # #plt.ylim((-r_max/2, r_max/2))
#
        # plt.gca().tick_params(direction='in', right=True, top=True)
        # plt.title(tit_str)
        # plt.colorbar(im)
        # #r_string = r'\log_{10}(r)\ ' if kwargs['logr'] else r'r\ '
        # r_string = r'r '
        # angle_string_x = r'\cos(\phi)' if kwargs['midplane'] else r'\sin(\theta)'
        # angle_string_y = r'\sin(\phi)' if kwargs['midplane'] else r'\cos(\theta)'
        # plt.xlabel('$' + r_string + angle_string_x + '$')
        # plt.ylabel('$' + r_string + angle_string_y + '$')
        # plt.tight_layout()
#
#
        # # ---------- naming and saving ----------
        # descript = "slice"
#
        # pdir = self.fig_save_fp + descript + "/"+self.config+"_"+self.specs+"/"
        # if kwargs['midplane']:
            # pdir += "midplane_slice/plots/"
            # descript += "_midplane"
        # else:
            # pdir += "vertical_slice/plots/"
            # descript += "_vertical"
        # if kwargs['logc']:
            # descript += "_logc"
        # if to_avg:
            # descript += "_avg-phi"
        # tstr = "_t{}".format(t)
        # rstr = "_r{:.1f}-{:.1f}".format(self.r[ri_min], self.r[ri_max])
        # pname = self.create_fname("plot", q, descript+tstr+rstr, ".png")
        # if to_save:
            # self.figsave(pdir, pname, overwrite)
#
#
    # #------------------------------------------------------------------------
#
    # def movie_q_slice(self, q, t_range, **kwargs):
        # """Plot either vertical or midplane slice of variable q at time t.
#
        # Keyword Arguments:
        # q: the variable to be queried. Should be in self.variables.
        # t: the timestep (i.e. number in filename, not index) to be queried.
        # to_avg: if true, average over phi.
        # ri_range: Index range of radii over which plot should be made. Default is full range.
        # fig_num: figure on which to plot. Default is 1.
        # overwrite: if true, overwrite files with the same name. If false, create file with "vx", where x is high enough to not overwrite any existing files.
#
        # """
#
        # # Deal with kwargs
        # if 'ri_max' not in kwargs:
            # ri_max = -1
        # else:
            # ri_max = kwargs['ri_max']
        # if 'plot_lims' in kwargs:
            # plot_lims = kwargs['plot_lims']
        # if 'fig_num' not in kwargs:
            # fig_num = 1
        # else:
            # fig_num = kwargs['fig_num']
        # if 'to_avg' not in kwargs:
            # to_avg = False
        # else:
            # to_avg = kwargs['to_avg']
        # if 'overwrite' not in kwargs:
            # overwrite = False
        # else:
            # overwrite = kwargs['overwrite']
        # if 'logc' not in kwargs:
            # kwargs['logc'] = True
#
        # # --- Get min/max values --------
        # [minval, maxval] = self.get_q_slice_extrema(q)
        # fmax = 10**9
        # if maxval > fmax:
            # maxval = fmax
        # valrange = maxval - minval
        # if self.r is None:
            # self.r = pickle.load(open(self.data_save_fp + "r_steps.p", "rb"))
        # r = self.r
        # r_max = r[ri_max]
#
        # if 'midplane' not in kwargs or not kwargs['midplane']:
            # tit_str = "Vertical slice of "
            # kwargs['midplane'] = False
        # else:
            # tit_str = "Midplane slice of "
        # tit_str += r"{}".format(self.varnames[q])
#
#
        # # -------- Load scalar grid, data --------
        # grid_fname = self.data_save_fp + "grid_"
        # if kwargs['midplane']:
            # grid_fname += "midplane.p"
        # else:
            # grid_fname += "vertical.p"
#
        # #vals = self.get_q_data_slice_at_t(q, t, midplane=kwargs['midplane'])
        # grid = pickle.load(open(grid_fname, "rb"))
        # x_grid = grid["x_grid"]
        # y_grid = grid["y_grid"]
#
#
        # # ---- colorbar mapping -----
#
        # if valrange <= 0:
            # print("Data is identically zero. Exiting")
            # return 0
#
        # if minval < 0:
            # if 'cmap' not in kwargs:
                # kwargs['cmap'] = 'RdBu'
            # # Make cmap symmetric about zero
            # vmin = -np.max([np.abs(minval), np.abs(maxval)])
            # vmax = -vmin
        # else:
            # if q == "beta":
                # kwargs['cmap'] = 'Blues_r'
            # else:
                # kwargs['cmap'] = 'Blues'
            # vmin = minval
            # vmax = maxval
#
        # fig = plt.figure()
        # cmap = plt.get_cmap(kwargs['cmap'])
        # if kwargs['logc']:
            # if minval <= 0:
                # norm = colors.SymLogNorm(linthresh=(maxval-minval)/50.)
            # else:
                # norm = colors.LogNorm()
        # else:
            # norm = colors.Normalize()
#
#
        # # ---- plot limits ----
        # plt.gca().set_aspect('equal')
        # if 'plot_lims' not in kwargs:
            # plot_lims = [[-r_max, r_max], [-r_max, r_max]]
        # plt.xlim(plot_lims[0])
        # plt.ylim(plot_lims[1])
#
        # plt.gca().tick_params(direction='in', right=True, top=True)
        # r_string = r'r '
        # angle_string_x = r'\cos(\phi)' if kwargs['midplane'] else r'\sin(\theta)'
        # angle_string_y = r'\sin(\phi)' if kwargs['midplane'] else r'\cos(\theta)'
        # plt.xlabel('$' + r_string + angle_string_x + '$')
        # plt.ylabel('$' + r_string + angle_string_y + '$')
        # cbar = None
#
        # print("making movie...")
        # # ------- Filename ----
        # tstr = "_t{}-{}".format(t_range[0], t_range[-1])
        # rstr = "_r{:.1f}-{:.1f}".format(self.r[0], r_max)
        # descript = "slice"
        # pdir = self.fig_save_fp + descript + "/"+self.config+"_"+self.specs+"/"
        # if kwargs['midplane']:
            # pdir += "midplane_slice/movies/"
            # descript += "_midplane"
        # else:
            # pdir += "vertical_slice/movies/"
            # descript += "_vertical"
        # if kwargs['logc']:
            # descript += "_logc"
        # if to_avg:
            # descript += "_avg-phi"
        # if not os.path.exists(pdir):
            # os.makedirs(pdir)
        # fname = self.create_fname("movie", q, descript+tstr+rstr, "_v0")
        # if not overwrite:
            # count = 0
            # while os.path.isfile(pdir+fname+".mp4"):
                # count += 1
                # fname = self.rreplace(fname, count)
            # print("new name to avoid overwrite: "+ fname)
            # full_descript = descript+"_"+self.config+"_"+self.specs
        # metadata = dict(title=full_descript, artist='LH', comment=full_descript)
        # writer = FFMpegWriter(fps=5, metadata=metadata)
#
        # # ---- Actually make the movie! --------
        # print(pdir+fname)
        # with writer.saving(fig, pdir+fname+".mp4", 200):
            # for t in t_range:
                # t_str = tit_str + " at t = {}".format(t)
                # plt.title(t_str)
                # data = self.get_q_data_slice_at_t(q, t, midplane=kwargs['midplane'])
                # print("Time: {}".format(t))
#
                # # ----- plot data -----
                # im = plt.pcolormesh(x_grid, y_grid, data, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm)
                # if cbar is None:
                    # cbar = plt.colorbar(im)
                # writer.grab_frame()
                # del data
                # self.data_raw = {}
        # print("movie saved")
#
#
    # #------------------------------------------------------------------------
    # def figsave(self, idir, iname, overwrite=False):
        # """Expanded function of savefig to make directory if need be, not overwrite.
#
        # Keyword Arguments:
        # idir: Directory to save figure in. Must have / at the end.
        # iname: Name to use for file. Will be added on to if file already exists (and overwrite is false).
        # overwrite: if true, write file with iname no matter what. Otherwise, append "vx" where x is a number sufficiently high so that no file with that name exists.
        # """
        # if not os.path.exists(idir):
            # os.makedirs(idir)
        # count = 0
        # if not overwrite:
            # iname = iname.split(".")
            # iname[-2] += "_v0"
            # iname = ".".join(iname)
            # print(idir+iname)
#
            # while os.path.isfile(idir+iname):
                # count += 1
                # iname = self.rreplace(iname, count)
            # print("new name to avoid overwrite: "+ iname)
        # plt.savefig(idir+iname, bbox_inches='tight')
        # return [idir, iname]
#
#
    # #------------------------------------------------------------------------
    # def rreplace(self, string, count):
        # """
        # Improvements:
        # -can't handle double digit count"""
        # rstring = string[::-1]
        # newrstring = rstring.replace(str(count-1), str(count), 1)
        # return newrstring[::-1]
#
    # #------------------------------------------------------------------------
    # def create_fname(self, ftype, q, descript, ext):
        # return ftype + "_" + q + "_" + self.config + "_" + self.specs + "_" + descript + ext
#
    # #------------------------------------------------------------------------
    # def deep_getsizeof(self, o, ids):
        # """Find the memory footprint of a Python object
#
        # This is a recursive function that drills down a Python object graph
        # like a dictionary holding nested dictionaries with lists of lists
        # and tuples and sets.
#
        # The sys.getsizeof function does a shallow size of only. It counts each
        # object inside a container as pointer only regardless of how big it
        # really is.
#
        # :param o: the object
        # :param ids:
        # :return:
        # """
        # if id(o) in ids:
            # return 0
#
        # r = sys.getsizeof(o)
        # ids.add(id(o))
#
        # if isinstance(o, str) or isinstance(0, str):
            # return r
#
        # if isinstance(o, Mapping):
            # return r + sum(self.deep_getsizeof(k, ids) + self.deep_getsizeof(v, ids) for k, v in o.items())
#
        # if isinstance(o, Container):
            # return r + sum(self.deep_getsizeof(x, ids) for x in o)
#
        # return r
