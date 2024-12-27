import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import subprocess
import threading


# Define the parameters
par_start = [
    -10, 1.82,
    40, 3.57,
    50, 400, 
    50, 6000
]

par_stop = [
    50, 2.12,
    100, 3.87, 
    200, 600,
    200, 6200
]

par_step = [
    2, 0.01, 
    2, 0.01,
    5, 5, 
    5, 5
]

par_optimal = [
    85.2619/4, 1.9796,
    288.154/4, 3.72882,
    85, 525,
    110, 6100
]

# Initialize the correlation matrices
H_Correlation = [[None for _ in range(8)] for _ in range(8)]
TOTAL = 0

# Create the correlation matrices
for first in range(8):
    for second in range(8):
        x_bins = int(abs(par_start[first] - par_stop[first]) / par_step[first])
        y_bins = int(abs(par_start[second] - par_stop[second]) / par_step[second])
        H_Correlation[first][second] = np.zeros((x_bins, y_bins))
        TOTAL += x_bins * y_bins

print(f"TOTAL: {TOTAL}")

# Create a figure with subplots
fig, axes = plt.subplots(8, 8, figsize=(16, 16))

start_time = time.time()
counter_progress = 0

# Function to minimize (dummy implementation)
def run_root_macro(params, thread):
    param_str = " ".join(map(str, params))

    command = f"./RBS_All_forMCMC {param_str} {thread}"

    print(command)

    result = subprocess.run(
        command,
        shell=True,  
        # stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, 
        # universal_newlines=True,
    )
    try:
        with open(f"tmp/{thread}.txt", "r") as f:
            lines = f.readlines()
            chi2 = float(lines[-1])

        # reading cout from command line
        # chi2 = float(result.stdout.split("\n")[-2])
        
    except (IndexError, ValueError):
        raise RuntimeError(f"Failed to parse chi2 from ROOT output: {result.stdout}")
    
    return chi2

# Fill the correlation matrices using multithreading
def process_parameters(first, second, i, par1, j, par2):
    par_it = par_optimal.copy()
    par_it[first] = par1
    par_it[second] = par2

    par = [
        par_it[4], par_it[5],
        par_it[6], par_it[7],
        par_it[0], par_it[1],
        par_it[2], par_it[3]
        
    ]

    Chi2 = run_root_macro(par,threading.get_ident())
    return first, second, i, j, Chi2


for first in range(8):
        for second in range(8):
            if first <= second:
                continue
            
            x_bins = int(abs(par_start[first] - par_stop[first]) / par_step[first])
            y_bins = int(abs(par_start[second] - par_stop[second]) / par_step[second])
            
            for i in range(x_bins):
                par1 = par_start[first] + i * par_step[first]   
                for j in range(y_bins):
                    par2 = par_start[second] + j * par_step[second]
                    chi2 = process_parameters(first, second, i, par1, j, par2)
                    H_Correlation[first][second][i, j] = chi2

                    counter_progress += 1
                    progress = (counter_progress / TOTAL) * 100
                    elapsed_time = time.time() - start_time
                    estimated_total_time = elapsed_time / (counter_progress / TOTAL)
                    remaining_time = estimated_total_time - elapsed_time
                    print(f"Progress: {progress:.2f}% - Estimated Remaining Time: {remaining_time:.2f}s", end='\r')

# Plot the correlation matrices
for first in range(8):
    for second in range(8):
        if first <= second:
            continue
        ax = axes[first, second]
        cax = ax.imshow(H_Correlation[first][second], norm=LogNorm(), aspect='auto', origin='lower')
        fig.colorbar(cax, ax=ax)

plt.tight_layout()
plt.show()