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
    [1.90, 2.0],
    [40, 100],
    [3.6, 3.8]
]

step_value = [
    5,
    5,
    5,
    5,
    1,
    0.01,
    1,
    0.01
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
    r"$a_{1.2}$",
    r"$b_{1.2}$",
    r"$a_{3.0}$",
    r"$b_{3.0}$"
]

fig, ax = plt.subplots(ndim, ndim, figsize = (10, 10))
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

thread = 0

for i, irange in enumerate(range_value):
    for j, jrange in enumerate(range_value):
        if i <= j:
            ax[i, j].axis('off')
            continue

        print(f"Processing {labels[i]} vs {labels[j]}")
        
        # loop on x and y for 2D plot
        for ix, x in enumerate(range(irange[0], irange[1]+step_value[i], step_value[i])):
            for iy, y in enumerate(range(jrange[0], jrange[1]+step_value[j], step_value[j])):
                
                print(f"{x}{unit_value[i]} and {y}{unit_value[j]}")
                
                parameter = minimum_value.copy()
                parameter[i] = x
                parameter[j] = y


                data_y[i][j].append(y)
                data_x[i][j].append(x)

                filename1 = f"../../../../../../mnt/hgfs/shared-2/ROOT_files/1.2_{parameter[0]}_{parameter[1]}_{parameter[0]}_A_-2.root"
                filename2 = f"../../../../../../mnt/hgfs/shared-2/ROOT_files/1.2_{parameter[0]}_{parameter[1]}_{parameter[0]}_B_5.root"
                filename3 = f"../../../../../../mnt/hgfs/shared-2/ROOT_files/3.0_{parameter[0]}_{parameter[1]}_{parameter[0]}_A_-2.root"
                filename4 = f"../../../../../../mnt/hgfs/shared-2/ROOT_files/3.0_{parameter[0]}_{parameter[1]}_{parameter[0]}_B_5.root"
                filename5 = f"../../../../../../mnt/hgfs/shared-2/ROOT_files/1.2_{parameter[2]}_{parameter[3]}_{parameter[2]}_A_-2.root"
                filename6 = f"../../../../../../mnt/hgfs/shared-2/ROOT_files/1.2_{parameter[2]}_{parameter[3]}_{parameter[2]}_B_5.root"
                filename7 = f"../../../../../../mnt/hgfs/shared-2/ROOT_files/3.0_{parameter[2]}_{parameter[3]}_{parameter[2]}_A_-2.root"
                filename8 = f"../../../../../../mnt/hgfs/shared-2/ROOT_files/3.0_{parameter[2]}_{parameter[3]}_{parameter[2]}_B_5.root"

                try:
                    file1 = TFile(filename1, "READ")
                    file2 = TFile(filename2, "READ")
                    file3 = TFile(filename3, "READ")
                    file4 = TFile(filename4, "READ")
                    file5 = TFile(filename5, "READ")
                    file6 = TFile(filename6, "READ")
                    file7 = TFile(filename7, "READ")
                    file8 = TFile(filename8, "READ")
                except:
                    data_chi2[i][j].append(0)
                    continue
 
                command = f"./RBS_All_forMCMC {parameter[0]} {parameter[1]} {parameter[2]} {parameter[3]} {parameter[4]} {parameter[5]} {parameter[6]} {parameter[7]} {thread}"
                subprocess.run(command, shell=True, stderr=subprocess.PIPE)

                try:
                    with open(f"tmp/{thread}.txt", "r") as f:
                        lines = f.readlines()
                        chi2 = np.sum([float(line) for line in lines])
                    data_chi2[i][j].append(chi2)
                except (IndexError, ValueError):
                    data_chi2[i][j].append(0)
                    continue
                    raise RuntimeError(f"Failed to parse chi2 from ROOT output")
                
        hist, xedges, yedges = np.histogram2d(data_x[i][j], data_y[i][j], bins=(range(range_value[i][0], range_value[i][1], step_value[i]), range(range_value[j][0], range_value[j][1], step_value[j])), weights=data_chi2[i][j])
        ax[i, j].imshow(hist.T, extent=[range_value[i][0], range_value[i][1], range_value[j][0], range_value[j][1]], origin='lower', aspect='auto')

        if j != 0: 
            ## deleting label on axis, not the ticks
            ax[i, j].set_yticklabels([])
        if i != ndim-1:
            ax[i, j].set_xticklabels([])
        if j == 0:
            ax[i, j].set_ylabel(f"{labels[i]}" + unit_value[i], fontsize=10)
        if i == ndim-1:
            ax[i, j].set_xlabel(f"{labels[j]}" + unit_value[j], fontsize=10)


        np.save("data_save1.npy", data)
        plt.savefig("correlation_plot.png")
plt.show()