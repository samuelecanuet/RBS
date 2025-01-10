from ROOT import TFile
from root2mpl import *
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import subprocess
import multiprocessing


Al = [5, 150]
Mylar = [400, 600]

al_data = []
mylar_data = []
weight_data = []

for al in range(Al[0], Al[1], 5):
    for mylar in range(Mylar[0], Mylar[1], 5):
        Problem = False
        al_data.append(al)
        mylar_data.append(mylar)
        filename = f"../../../../../../mnt/hgfs/shared-2/ROOT_files/1.2_{al}_{mylar}_{al}_A_-2.root"
        try:
            file = TFile(filename, "READ")
        except:
            weight_data.append(0)
            command = f"./RBS_All_forMCMC -1111 -1111 {int(al)} {int(mylar)} 21 1.9697 71 3.6 1"
            # subprocess.run(command, shell=True, stderr=subprocess.PIPE)
            continue

        try:
            hist = file.Get("RBS_1_6")
        except:
            weight_data.append(0)
            continue

        # if hist.GetMaximum() / hist.Integral() > 0.0036:
        #     # os.remove(filename)
        #     weight_data.append(2)
        # else:
        weight_data.append(1)

hist, xedges, yedges = np.histogram2d(al_data, mylar_data, bins=(range(Al[0], Al[1], 5), range(Mylar[0], Mylar[1], 5)), weights=weight_data)

plt.imshow(hist.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], origin='lower', aspect='auto')
plt.colorbar()
plt.show()
