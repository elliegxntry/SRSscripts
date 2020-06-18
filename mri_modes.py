"""
This script extracts the magnetic modes from each direction of the magnetic field by doing an fft.
Also includes plots of the first mode with specs and a plot of every mode calculated on the same plot
INPUTS:
    - Radius
    - Quantity (which direction of mag field)
    - distribution
    - time steps to run
    - which modes to calculate
OUTPUTS:
    - A .txt file with the frequency each mode appears for every timestep
    - A plot of the m1 frequencies over time
    - A plot of every mode calculated over time as lines
"""

# Import packages
import numpy as np
import os
import sys
sys.path.append("C:/Users/Ellie/Downloads/nerd/scripts/modules")
import new_athena_read as read
import matplotlib.pyplot as plt

# Specifications
radius = 10.0
theta = np.pi/2
quantity = 'Bcc1'
dist = 'B'
time_steps = np.arange(0, 671)
wavelengths = None
mstart = 1
mend = 3

# Paths to load and save data
config = "1.1.1-torus2_b-gz2_a0beta500tor" + dist + "_br32x32x64rl2x2"
datapath = "C:/Users/Ellie/Downloads/nerd/SRSData/"
data_load_path = datapath + config + "/"
data_save_path = "C:/Users/Ellie/Downloads/nerd/SRSData/Reduced/Constant" + dist + "/Fourier/"
if not os.path.isdir(data_save_path):
    os.mkdir(data_save_path)
figsavepath = "C:/Users/Ellie/Downloads/nerd/SRSPlots/Profiles/time_quantity_profiles/timemodes/" + dist + "/"
if not os.path.isdir(figsavepath):
    os.mkdir(figsavepath)
amp_path = data_save_path + "fourier-amplitude-" + quantity + "_data_r{}.txt".format(radius)
header = ""

# Calculate modes for each timestep
for time in time_steps:
    # load data
    fname = data_load_path + config + ".prim.{:05d}.athdf".format(time)
    #print("Loading time {}".format(time))

    # read in data only as function of phi
    raw_data = read.athdf(fname,quantities=[quantity],x1_min=radius, x1_max=radius, x2_min=theta, x2_max=theta+0.01) 
    phi_values = raw_data['x3v']
    N_phi = len(phi_values)
    phi_profile = raw_data[quantity].reshape((N_phi,))
    dphi = phi_values[1] - phi_values[0]

    # Do a spatial fourier transform of the values as a function of phi
    data_fft = np.fft.fft(phi_profile)/N_phi # N_phi normalizes output
    # create list of wavenumbers and wavelengths that match the fft
    if wavelengths is None:
        wavenumbers = np.linspace(0, 1.0/dphi, N_phi)
        wavelengths = 1.0/wavenumbers

    # # The below commented code plots the entire fourier spectrum. Do NOT plot for every time step!
    # plt.plot(wavelengths[:N_phi//2]*N_phi/(2*np.pi), np.abs(data_fft)[:N_phi//2], marker="o")
    # # # plt.plot(0.0, np.abs(data_fft)[:N_phi//2][0], marker="o")
    # plt.xlabel("Wavelength in cells")
    # plt.show()

    # Extract relevant mode amplitudes
    const_amplitude = np.abs((data_fft[:N_phi//2])[0]) # Magnetic field's constant value
    m1_amplitude = np.abs((data_fft[:N_phi//2])[1]) # The m=1 mode (i.e. sin(phi/2))
    m2_amplitude = np.abs((data_fft[:N_phi//2])[2]) # The m=2 mode (i.e. sin(phi))

    # Save to file
    code_time = np.array(raw_data['Time']).reshape(1,)
    output_data = np.hstack([code_time, np.abs(data_fft[:N_phi//2])]).reshape((1, N_phi//2 + 1))
    if not os.path.isfile(amp_path):
        header = "time, const amplitude"
        for i in np.arange(1, N_phi//2):
            header += ", m{} amplitude".format(i)
    with open(amp_path, "a") as f:
        np.savetxt(f, output_data, header=header)

# clean up mdot by removing duplicates and sorting
with open(amp_path, "r") as f:
    ampdata = np.loadtxt(f, skiprows=1)
unique_times, inds = np.unique(ampdata[:, 0], return_index=True)
sorted_ampdata = ampdata[inds]
# new_ampdata = np.transpose(np.stack([unique_times, ampdata[:,1][inds], ampdata[:,2][inds], ampdata[:,3][inds]]))
# print(new_ampdata.shape)
# print(sorted_ampdata.shape)
header = "time, const amplitude"
N_phi = 128
for i in np.arange(1, N_phi//2):
    header += ", m{} amplitude".format(i)
with open(amp_path, "w") as f:
    np.savetxt(f, sorted_ampdata, header=header)

# Plot m1 over time
plt.figure()
plt.plot(unique_times, sorted_ampdata[:, 2], '+')
plt.xlabel("Time [GM/c^3]")
plt.ylabel("M1")
plt.title(quantity + "\nConstant " + dist + "\n radius={} [GM/c^2]".format(radius) + config)
plt.tight_layout()

# Plot several modes over time
plt.figure()
# if N_phi is None: N_phi = 128
# mend = N_phi//2 - 1
num_plots = mend - mstart + 1
colormap = plt.cm.viridis
plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.viridis(np.linspace(0, 1, num_plots))))
labels=[]
for i in np.arange(mstart, mend + 1, 1):
    plt.plot(unique_times, sorted_ampdata[:, i+1])
    labels.append(r"$m = {}$".format(i))
plt.xlabel("Time [GM/c^3]")
plt.legend(labels, ncol=8)
plt.ylabel("Mode amplitude")
plt.title(quantity + "\nConstant " + dist + "\n radius={} [GM/c^2]".format(radius))
plt.tight_layout()
plt.savefig(figsavepath + quantity + "_r{}.png".format(radius))
plt.show()
