import numpy as np
import h5py
import os
import athena_read
import pickle

def save_constants(fpath, data, overwrite):
    """Save data under qname->t in fpath .hdf5 file"""
    # first, see if the file exists

    (nx1, nx2, nx3) = np.shape(data["rho"])

    with h5py.File(fpath, 'a') as f:
        if "x1v" not in f.keys():
            f.create_dataset("x1v", data=data["x1v"])
        elif overwrite:
            f["x1v"][...] = data["x1v"]
        if "x2v" not in f.keys():
            f.create_dataset("x2v", data=data["x2v"])
        elif overwrite:
            f["x2v"][...] = data["x2v"]
        if "x3v" not in f.keys():
            f.create_dataset("x3v", data=data["x3v"])
        elif overwrite:
            f["x3v"][...] = data["x3v"]
        if "mp_grid" not in f.keys() or overwrite:
            x3v = data["x3v"]
            phi_extended = \
                np.concatenate((x3v[-1:] - 2.0 * np.pi, x3v, x3v[:1] + 2.0 * np.pi))
            phi_extended -= 0.5 * 2.0 * np.pi / nx3

            r_grid, phi_grid = np.meshgrid(data["x1v"], phi_extended)

            mp_grid = [r_grid*np.cos(phi_grid), r_grid*np.sin(phi_grid)]
            try:
                f.create_dataset("mp_grid", data=mp_grid)
            except:
                del f["mp_grid"]
                f.create_dataset("mp_grid", data=mp_grid)
        if "vert_grid" not in f.keys() or overwrite:
            x2v = data["x2v"]
            theta_extended = np.concatenate((-x2v[0:1], x2v, 2.0 * np.pi - x2v[::-1],
                                             2.0 * np.pi + x2v[0:1]))

            theta_extended_corrected = theta_extended - 0.5 * (x2v[1] - x2v[0])

            r_grid, theta_grid = np.meshgrid(data["x1v"], theta_extended_corrected)
            vert_grid = [r_grid*np.sin(theta_grid), r_grid*np.cos(theta_grid)]
            try:
                f.create_dataset("vert_grid", data=vert_grid)
            except:
                del f["vert_grid"]
                f.create_dataset("vert_grid", data=vert_grid)


        # XX make grids
        # theta
        # self.x2vcorr = theta_extended_corrected
        # theta_extended2 = np.concatenate((-self.x2f[0:1], self.x2f, 2.0 * np.pi - self.x2f[::-1],
                                        # 2.0 * np.pi + self.x2f[0:1]))
        # theta_extended_corrected2 = theta_extended2 - 0.5 * (self.x2f[1] - self.x2f[0])
        # self.x2fcorr = theta_extended_corrected2
        # phi_extended = \
            # np.concatenate((self.x3v[-1:] - 2.0 * np.pi, self.x3v, self.x3v[:1] + 2.0 * np.pi))
        # phi_extended -= 0.5 * 2.0 * np.pi / self.nx3

        # --- calculate spherical slice grids ---
        # theta-phi plane (weird)
        # self.x1s_grid = np.meshgrid(theta_extended_corrected, phi_extended)
        # r-phi plane (midplane slice)
        # r_grid, phi_grid = np.meshgrid(self.x1v, phi_extended)
        # r_grid, phi_grid = np.meshgrid(x1v, x3v)
# 
        # x2s_grid = [r_grid*np.cos(phi_grid), r_grid*np.sin(phi_grid)]
        # r-theta plane (vertical slice)
        # r_grid, theta_grid = np.meshgrid(self.x1v, theta_extended_corrected)
        # r_grid, theta_grid = np.meshgrid(x1v, x2v)
        # x3s_grid = [r_grid*np.sin(theta_grid), r_grid*np.cos(theta_grid)]


    return


def save_to_hdf5(fpath, qname, t, data, overwrite):
    """Save data under qname->t in fpath .hdf5 file"""
    # first, see if the file exists
    with h5py.File(fpath, 'a') as f:
        if qname not in f:
            f.create_group(qname)
        if str(t) not in f[qname]:
            print("Saving " + qname + " at time {} to ".format(t) + fpath)
            f[qname].create_dataset(str(t), data=data)
        elif overwrite:
            print("Saving " + qname + " at time {} to ".format(t) + fpath)
            del f[qname][str(t)]
            f[qname].create_dataset(str(t), data=data)
        else:
            print(qname + " already saved for time {}".format(t))
    return

def calculate_qdata(q, data):
    variables = ['press', 'Bcc1', 'Bcc2', 'Bcc3', 'rho', 'vel1', 'vel2', 'vel3']
    calc_vars = ['invbeta', 'qtheta', 'qphi']
    if q in variables:
        return data[q]
    if q not in calc_vars:
        print("Error: variable name " + q+ " is not a valid option")
    elif q == "invbeta":
        return calculate_invbeta(data)
    elif q == "qtheta":
        return calculate_qtheta(data)
    elif q == "qphi":
        return calculate_qphi(data)


def calculate_qtheta(data):
    """Calculates qtheta at time t
    data should hold Bcc2, vel3, and rho 3d data in a dictionary"""
    numerator = 2*np.pi*data["Bcc2"]
    denominator = np.multiply(data["vel3"], np.sqrt(data["rho"]))
    qtheta = np.divide(numerator, denominator)
    return qtheta

def calculate_qphi(data):
    """Calculates qphi at time t
    data should hold Bcc2, vel3, and rho 3d data in a dictionary"""
    numerator = 2*np.pi*data["Bcc3"]
    denominator = np.multiply(data["vel3"], np.sqrt(data["rho"]))
    qtheta = np.divide(numerator, denominator)
    return qtheta

def calculate_invbeta(data):
    """Calculates invbeta at time t
    data should hold Bcc1, Bcc2, Bcc3, and press 3d data in a dictionary"""
    mag_pressure = (data["Bcc1"]**2+data["Bcc2"]**2+data["Bcc3"]**2)/2.0 # note that 4PI = 1
    invbeta = mag_pressure/data["press"]
    return invbeta

def get_slice(data, midplane=False, average=False):
    (nx1, nx2, nx3) = np.shape(data)
    if midplane:
        if nx2 % 2 == 0:
            slice_data = np.mean(data[:, int(nx2 /2 - 1):int(nx2 /2+1), :], axis=1)
        else:
            slice_data = data[:, int(nx2/2), :]
        # join through boundaries
        slice_data = np.vstack((slice_data[-1:, :], slice_data, slice_data[:1, :]))
    else:
        if average:
            vals_right = np.mean(data, axis=0)
            vals_left = vals_right
        else:
            vals_right = 0.5*(data[-1,:,:] + data[0,:,:])
            vals_left = 0.5*(data[int((nx3/2)-1), :,:] + data[int(nx3/2), :, :])
        slice_data = np.vstack((vals_left[:1, :], vals_right,
                                vals_left[::-1, :], vals_right[:1, :]))
    return slice_data

def calculate_extrema(data):
    return

def save_extrema(fpath, qname, t, extrema_data, overwrite=False):
    """add extrema for different types of data to file.

    Load extrema from file, calculate new extrema if necessary and write new extrema to csv file.
    XX To do: extend to add different types of extrema without re-calculating all of others.
    """
    header = ["vol_min", "vol_max", "mp_min", "mp_max", "vert_min", "vert_max"]

    # if file doesn't exist, export header. if it does, try to load data
    if not os.path.isfile(fpath):
        extrema = {}
    else:
        with open(fpath, 'rb') as f:
            extrema = pickle.load(f)

    # qname, then time, then string
    if qname in extrema.keys():
        if str(t) in extrema[qname].keys():
            if not overwrite:
                print(qname + " extrema already saved at time {}".format(t))
            else:
                print("Overwriting " + qname + " at time {} to ".format(t) + fpath)
                extrema[qname][str(t)] = extrema_data
        else:
            print("Saving " + qname + " at time {} to ".format(t) + fpath)
            extrema[qname][str(t)] = extrema_data
    else:
        extrema[qname] = {}
        print("Saving " + qname + " at time {} to ".format(t) + fpath)
        extrema[qname][str(t)] = extrema_data

    # Now get the max/min over all time
    if "all" not in extrema[qname]:
        all_time = np.zeros((8, 1))
    else:
        all_time = extrema[qname]["all"]

    if extrema_data[0] < all_time[0]: all_time[0] = extrema_data[0]
    if extrema_data[1] > all_time[1]: all_time[1] = extrema_data[1]
    if extrema_data[2] < all_time[2]: all_time[2] = extrema_data[2]
    if extrema_data[3] > all_time[3]: all_time[3] = extrema_data[3]
    if extrema_data[4] < all_time[4]: all_time[4] = extrema_data[4]
    if extrema_data[5] > all_time[5]: all_time[5] = extrema_data[5]
    if extrema_data[6] < all_time[6]: all_time[6] = extrema_data[6]
    if extrema_data[7] > all_time[7]: all_time[7] = extrema_data[7]

    extrema[qname]["all"] = all_time

    # Save to file!
    with open(fpath, 'wb') as f:
        pickle.dump(extrema, f)
    return extrema
