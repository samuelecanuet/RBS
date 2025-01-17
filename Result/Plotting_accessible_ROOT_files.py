from ROOT import TFile
from root2mpl import *
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import subprocess
import multiprocessing

#1.2B 0.005
#3.0B 0.009

#3.0A 0.007
#1.2A 0.0045

Al = [5, 150]
Mylar = [5900, 6200]

al_data = []
mylar_data = []
weight_data = []

fig, ax = plt.subplots()

for al in range(Al[0], Al[1], 5):
    for mylar in range(Mylar[0], Mylar[1], 5):
        Problem = False
        al_data.append(al)
        mylar_data.append(mylar)
        filename = f"../../../../../../mnt/hgfs/shared-2/ROOT_files/3.0_{al}_{mylar}_{al}_A_-2.root"
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
            # os.remove(filename)
            weight_data.append(2)
        else:
            weight_data.append(1)


plt.show()

plt.close("all")



hist, xedges, yedges = np.histogram2d(al_data, mylar_data, bins=(range(Al[0], Al[1], 5), range(Mylar[0], Mylar[1], 5)), weights=weight_data)
plt.imshow(hist.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], origin='lower', aspect='auto')
plt.colorbar()
plt.show()
