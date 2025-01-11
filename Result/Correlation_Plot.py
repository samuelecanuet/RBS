from ROOT import TFile
from root2mpl import *
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import subprocess
import multiprocessing

###correlation plot for 8 paramaters

minimum_value = [85, 525, 110, 6100, 22.0887, 1.9681, 72.0385, 3.72882]

range_value = [
    [50, 150],
    [400, 600],
    [50, 150],
    [5900, 6200],
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
    r"$\epsilon$",
]

fig, ax = plt.subplots(8, 8, figsize=(20, 20))
fig.tight_layout(pad=2.0)

thread = 1
for i, irange in enumerate(range_value):
    for j, jrange in enumerate(range_value):
        if i >= j:
            continue
        
        # loop on x and y for 2D plot
        for x in range(irange[0], irange[1], step_value[i]):
            for y in range(jrange[0], jrange[1], step_value[j]):
                
                parameter = minimum_value.copy()
                parameter[i] = x
                parameter[j] = y
 
                command = f"./RBS_All_forMCMC {parameter[0]} {parameter[1]} {parameter[2]} {parameter[3]} {parameter[4]} {parameter[5]} {parameter[6]} {parameter[7]} {thread}"
                # subprocess.run(command, shell=True, stderr=subprocess.PIPE)

                try:
                    with open(f"tmp/{thread}.txt", "r") as f:
                        lines = f.readlines()
                        chi2 = float(lines[-1])
        
                except (IndexError, ValueError):
                    raise RuntimeError(f"Failed to parse chi2 from ROOT output: {result.stdout}")




Al = [5, 150]
Mylar = [400, 600]

al_data = []
mylar_data = []
weight_data = []

fig, ax = plt.subplots()

for al in range(Al[0], Al[1], 5):
    for mylar in range(Mylar[0], Mylar[1], 5):
        Problem = False
        al_data.append(al)
        mylar_data.append(mylar)
        filename = f"../../../../../../mnt/hgfs/shared-2/ROOT_files/3.0_{al}_{mylar}_{al}_B_5.root"
        try:
            file = TFile(filename, "READ")
        except:
            weight_data.append(0)
            command = f"./RBS_All_forMCMC -1111 -1111 {int(al)} {int(mylar)} 21 1.9697 71 3.6 1"
            # subprocess.run(command, shell=True, stderr=subprocess.PIPE)
            continue

        try:
            hist = file.Get("RBS_1_6")
            if hist.GetMaximum() / hist.Integral() > 0.007:
                DisplayTH1D(hist, ax = ax, normalized=True, color='red')
            else:
                DisplayTH1D(hist, ax = ax, normalized=True, color='blue')
        except:
            weight_data.append(0)
            continue

        if hist.GetMaximum() / hist.Integral() > 0.007:
        #     # os.remove(filename)
            weight_data.append(2)
        else:
            weight_data.append(1)


plt.show()

plt.close("all")



hist, xedges, yedges = np.histogram2d(al_data, mylar_data, bins=(range(Al[0], Al[1], 5), range(Mylar[0], Mylar[1], 5)), weights=weight_data)
plt.imshow(hist.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], origin='lower', aspect='auto')
plt.colorbar()
plt.show()
