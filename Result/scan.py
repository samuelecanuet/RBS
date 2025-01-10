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


# Function to minimize (dummy implementation)
def run_root_macro(params):
    param_str = " ".join(map(str, params))

    command = f"./RBS_2SCAN {param_str}"

    result = subprocess.run(
        command,
        shell=True,  
        # stdout=subprocess.PIPE,
        # stderr=subprocess.PIPE, 
        # universal_newlines=True,
    )


def process_pair(i, j):
    run_root_macro([i, j])

# Use ThreadPoolExecutor to multithread the loop
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as executor:
    futures = []
    for i in range(8):
        for j in range(8):
            if i <= j:
                continue
            else:
                futures.append(executor.submit(process_pair, i, j))

    # Optionally, wait for all futures to complete
    for future in concurrent.futures.as_completed(futures):
        future.result()  # This will raise any exceptions caught during execution


## create 8x8 grid of plots
fig, axs = plt.subplots(8, 8, figsize=(20, 20))
for i in range(8):
    for j in range(8):
        if i <= j:
            axs[i, j].axis("off")
            continue
        try:
            f = ROOT.TFile(f"saved_scan/Correlation_Saved_{i}_{j}.root")
            h = f.Get(f"Correlation_{i}_{j}")
            DisplayTH2D(h, ax=axs[i, j])
        except:
            axs[i, j].axis("off")
            continue
        

plt.show()
