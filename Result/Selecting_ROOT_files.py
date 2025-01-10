from ROOT import TFile
from root2mpl import *
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import subprocess
import multiprocessing

def int_to_color(integer, colormap_name='viridis', vmin=None, vmax=None):
    cmap = cm.get_cmap(colormap_name)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    normalized_value = norm(integer)
    return cmap(normalized_value)

def process_file(i, var, step, n, energy, Al, Mylar, face, offset):
    mod = int((var - step * n / 2) + i * step)

    # Construct the filename based on the variable
    if var == Mylar:
        filename = f"../../../../../../mnt/hgfs/shared-2/ROOT_files/{energy}_{Al}_{mod}_{Al}_{face}_{offset}.root"
    else:
        filename = f"../../../../../../mnt/hgfs/shared-2/ROOT_files/{energy}_{mod}_{Mylar}_{mod}_{face}_{offset}.root"

    # Try opening the file
    try:
        file = TFile(filename, "READ")
    except:
        return

    # Try getting the histogram
    try:
        hist = file.Get("RBS_1_6")
    except:
        return

    # Check conditions and remove the file if necessary
    if hist.GetMaximum() / hist.Integral() > 0.0036:
        os.remove(filename)

        if var == Mylar:
            command = f"./RBS_All_forMCMC 85 525 {Al} {mod} 21 1.9697 71 3.6 0"
        else:
            command = f"./RBS_All_forMCMC 85 525 {mod} {Mylar} 21 1.9697 71 3.6 0"

        subprocess.run(command, shell=True, stderr=subprocess.PIPE)

def main():
    # Input parameters
    energy = 1.2
    Al = 100
    Mylar = 6145
    face = "B"
    offset = 5.25
    step = 5
    n = 100
    var = Mylar


    # Create a pool of workers
    with multiprocessing.Pool(processes=6) as pool:
        # Distribute the work among the processes
        pool.starmap(
            process_file,
            [(i, var, step, n, energy, Al, Mylar, face, offset) for i in range(n)]
        )

if __name__ == "__main__":
    main()
