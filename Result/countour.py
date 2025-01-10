import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import subprocess
import threading
from root2mpl import *
import ROOT
import concurrent.futures


import numpy as np
import corner
import matplotlib.pyplot as plt

# Load the chain and log_prob files
chain = np.load("chain_thin_1p2_4par_40w_binned.npy")  # Shape: (nsteps, nwalkers, ndim)
log_prob = np.load("log_prob_thin_1p2_4par_40w_binned.npy")  # Shape: (nsteps, nwalkers)

# Flatten the chain and log_prob arrays
chain_flat = chain.reshape(-1, chain.shape[-1])  # Shape: (nsteps * nwalkers, ndim)
log_prob_flat = log_prob.flatten()  # Shape: (nsteps * nwalkers,)

# Convert log_prob to chi2
chi2 = -2 * log_prob_flat

# Find the minimum chi2
chi2_min = np.min(chi2)
# print correspondant parameter values
min_index = np.argmin(chi2)
min_params = chain_flat[min_index]
print(f"Minimum chi2: {min_params}")

# Filter parameter sets satisfying chi2 <= chi2_min + 1
condition = chi2 >= (chi2_min + 1 - 0.1)
filtered_chain1 = chain_flat[condition]
filtered_chi2 = chi2[condition]
condition = filtered_chi2 <= (chi2_min + 1 + 0.1)
filtered_chain = filtered_chain1[condition]

min_values = np.min(filtered_chain, axis=0)
max_values = np.max(filtered_chain, axis=0)
bin_values = [5, 5, 1, 0.002]

custom_bins = [
    28, 
    44,
    100,
    100
]

for i in range(4):
    custom_bins[i] = round((max_values[i]-min_values[i])/bin_values[i])

print(f"Original chain shape: {chain_flat.shape}")
print(f"Filtered chain shape: {filtered_chain.shape}")
print(f"Chi2 min: {chi2_min}")

# Parameter labels for the corner plot
labels = [f"param{i+1}" for i in range(chain.shape[-1])]  # Example: param1, param2, ...

# Generate and save the corner plot
fig = corner.corner(filtered_chain, labels=labels, bins = custom_bins)
fig.savefig("filtered_corner_chi2min+1.png")

print("Corner plot saved as 'filtered_corner_chi2min+1.png'.")
