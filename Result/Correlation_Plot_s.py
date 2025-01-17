from ROOT import TFile
from root2mpl import *
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import subprocess
import multiprocessing

###correlation plot for 8 paramaters

ndim = 8

minimum_value = [85, 525, 110, 6100, 22.0887, 1.9681, 72.0385, 3.72882]

range_value = [
    [50, 150],
    [400, 600],
    [50, 150],
    [6000, 6100],
    [0, 40],
    [1.95, 2.0],
    [40, 100],
    [3.6, 3.8]
]

step_value = [
    5,
    5,
    5,
    5,
    1,
    0.005,
    1,
    0.005
]

step_value = [ step_value[i] for i in range(ndim) ]

unit_value = [
    " [nm]", 
    " [nm]",
    " [nm]",
    " [nm]",
    " [keV]",
    " [keV/CH]",
    " [keV]",
    " [keV/CH]"
]

labels = [
    r"$\epsilon_{Al}$",
    r"$\epsilon_{Mylar}$",
    r"$\epsilon_{Al}$",
    r"$\epsilon_{Mylar}$",
    r"$b_{1.2}$",
    r"$a_{1.2}$",
    r"$b_{3.0}$",
    r"$a_{3.0}$"
]

fig, ax = plt.subplots(ndim, ndim, figsize = (10, 10))
# figg, axx = plt.subplots()
# fig.tight_layout(pad=1.0)
for i, irange in enumerate(range_value):
    for j, jrange in enumerate(range_value):
        if i <= j:
            ax[i, j].axis('off')
            continue


## 8x8 list 
data_x = [
    [[] for _ in range(ndim)] for j in range(ndim)
]

data_y = [
    [[] for _ in range(ndim)] for j in range(ndim)
]

data_chi2 = [
    [[] for _ in range(ndim)] for j in range(ndim)
]

chi2_min = [
    3.36, 3.51, 2.14, 3.28, 2.22, 2.82, 1.25, 1.67
]

thread = 1
import numpy as np
import os
from ROOT import TFile
from root2mpl import *
import matplotlib.pyplot as plt
import subprocess
from scipy.ndimage import gaussian_filter
# File to save/load data
save_file = "data_save1.npy"

# Load existing data or initialize
if os.path.exists(save_file):
    saved_data = np.load(save_file, allow_pickle=True).tolist()
    saved_dict = {tuple(entry[0]): entry[1] for entry in saved_data}  # Convert to dict for fast lookup
else:
    saved_data = np.array([], dtype=[('parameters', object), ('chi2', float)])
    saved_dict = {}

# Initialize the rest of the arrays and settings as before
data_x = [[[] for _ in range(ndim)] for j in range(ndim)]
data_y = [[[] for _ in range(ndim)] for j in range(ndim)]
data_chi2 = [[[] for _ in range(ndim)] for j in range(ndim)]

for i, irange in enumerate(range_value):
    for j, jrange in enumerate(range_value):
        if i <= j:
            ax[i, j].axis("off")
            continue

        print(f"Processing {labels[i]} vs {labels[j]}")

        for ix, x in enumerate(np.arange(irange[0], irange[1] + step_value[i], step_value[i])):
            for iy, y in enumerate(np.arange(jrange[0], jrange[1] + step_value[j], step_value[j])):
                x = round(x, 4)
                print(f"{x}{unit_value[i]} and {y}{unit_value[j]}")

                # Create parameter set
                parameter = minimum_value.copy()
                parameter[i] = x
                parameter[j] = y

                # Check if the parameter set is already saved
                param_tuple = tuple(parameter)
                if param_tuple in saved_dict:
                    chi2 = saved_dict[param_tuple]
                    # print(f"Found saved value for {param_tuple}: {chi2}")
                else:
                    # Compute chi2 as before
                    command = f"./RBS_All_forMCMC {parameter[0]} {parameter[1]} {parameter[2]} {parameter[3]} {parameter[4]} {parameter[5]} {parameter[6]} {parameter[7]} {thread}"
                    subprocess.run(command, shell=True, stderr=subprocess.PIPE)

                    try:
                        with open(f"tmp/{thread}.txt", "r") as f:
                            lines = f.readlines()
                            chi2 = np.sum([float(line) for line in lines])
                    except (IndexError, ValueError):
                        chi2 = 0

                    # Save new value
                    saved_data.append([param_tuple, chi2])
                    saved_dict[param_tuple] = chi2

                # Update data arrays
                data_x[i][j].append(x)
                data_y[i][j].append(y)
                data_chi2[i][j].append(chi2)

        # Save progress after each subplot
        np.save(save_file, saved_data)

        # Generate and save histogram
        hist, xedges, yedges = np.histogram2d(
            data_y[i][j],
            data_x[i][j],
            bins=(
                np.arange(jrange[0], jrange[1], step_value[j]),
                np.arange(irange[0], irange[1], step_value[i])
            ),
            weights=data_chi2[i][j]
        )

        smoothed_hist = gaussian_filter(hist, sigma=1)
        img = ax[i, j].imshow(smoothed_hist.T, extent=[jrange[0], jrange[1], irange[0], irange[1]],
                        origin='lower', aspect='auto', vmax=50, vmin=10, cmap='ocean')
        
        ax[i, j].set_ylim(irange[0]+step_value[i], irange[1]-step_value[i])
        ax[i, j].set_xlim(jrange[0]+step_value[j], jrange[1]-step_value[j])
        
        ## line x and y of truth value
        ax[i, j].plot([jrange[0], jrange[1]], [minimum_value[i], minimum_value[i]], color='red')
        ax[i, j].plot([minimum_value[j], minimum_value[j]], [irange[0], irange[1]], color='red')
        ax[i, j].scatter(minimum_value[j], minimum_value[i], color='red', s=10)

        ## plotting countour line
        # ax[i, j].contour(smoothed_hist.T, levels=[chi2_min[i] + 2.3, chi2_min[i] + 6.18, chi2_min[i] + 11.83], extent=[jrange[0], jrange[1], irange[0], irange[1]], colors='red')

        if j != 0:
            ax[i, j].set_yticklabels([])
        if i != ndim - 1:
            ax[i, j].set_xticklabels([])
        if j == 0:
            ax[i, j].set_ylabel(f"{labels[i]}" + unit_value[i], fontsize=10)
        if i == ndim - 1:
            ax[i, j].set_xlabel(f"{labels[j]}" + unit_value[j], fontsize=10)

        plt.savefig("correlation_plot1.png")
        print(np.min(data_chi2[i][j]), np.max(data_chi2[i][j]))
        print(data_x[i][j][np.argmin(data_chi2[i][j])], data_y[i][j][np.argmin(data_chi2[i][j])])
        print(labels[i], labels[j])
        