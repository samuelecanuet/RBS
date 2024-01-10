import numpy as np
from matplotlib import pyplot as plt
import csv
from ROOT import *
import subprocess
import argparse
import ctypes
import seaborn as sns
from scipy.optimize import *


def DisplayTH1D(Hist, ax, color=None, label=None, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, xtick=None, ytick=None, ylog=None, xlog=None, rebin=None, normalized=None):
        if rebin   != None: Hist.Rebin(rebin)
        if normalized == True: integral = Hist.Integral()
        else : integral = 1.
        nbins_x = Hist.GetNbinsX()
        hist_data = np.zeros(nbins_x)
        bin_centers_x = np.zeros(nbins_x)

        for i in range(1, nbins_x + 1):
            hist_data[i - 1] = Hist.GetBinContent(i)/integral
            bin_centers_x[i - 1] = Hist.GetXaxis().GetBinCenter(i)

        

        if color    == None: color = "black"
        if label    == None: label = Hist.GetTitle()               
        if title    == None: title = Hist.GetTitle()
        if xlabel   == None: xlabel = Hist.GetXaxis().GetTitle()
        if ylabel   == None: ylabel = Hist.GetYaxis().GetTitle()
        if normalized==True: ylabel = "Normalized" + ylabel
        if xlim     == None: xlim = ( bin_centers_x.min(), bin_centers_x.max() )
        if ylim     == None and hist_data.max()*1.1 > ax.get_ylim()[1] : ylim = ( 0, hist_data.max()*1.1 )
        if xtick    != None: ax.set_xticks(np.linspace(xlim[0], xlim[1], xtick))
        if ytick    != None: ax.set_yticks(np.linspace(ylim[0], ylim[1], ytick))
        if xlog     != None: ax.set_xscale('log')
        if ylog     != None: 
            ax.set_yscale('log')
            ylim = ( 1, hist_data.max()*1.1 )
        
        # ax.bar(bin_centers_x, hist_data, label = label, color=color)
        ax.bar(bin_centers_x, hist_data, label="SRIM", edgecolor='black', color='white', linewidth=0.5, width=bin_centers_x[0]-bin_centers_x[1])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        return [bin_centers_x, hist_data]

def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-(x - mean)**2 / (2 * stddev**2))



# fig, ax = plt.subplots()
# rootfile = TFile("test.root")
# histo = DisplayTH1D(rootfile.Get("RBS"), ax = ax)

# dic={}

# def which_mat(valeur, dictionnaire):
#     for cle, intervalle in dictionnaire.items():
#         if intervalle[0] <= valeur <= intervalle[1]:
#             return cle
#     return None

# for name in rootfile.GetListOfKeys():
#     if "G4" in name.GetTitle():
#         h1 = rootfile.Get(name.GetTitle())
#         total_content = 0
#         bin_list = []
#         for bin in range(1, h1.GetNbinsX() + 1):
#             bin_content = h1.GetBinContent(bin)
            
#             if (bin_content != 0.0):
#                 total_content += bin_content
#                 bin_list.append(bin)
#         dic[total_content / len(bin_list)] = [h1.GetBinLowEdge(min(bin_list)), h1.GetBin(max(bin_list))]
        

# print(dic)
# value=[]
# tree = rootfile.Get("Catcher")
# histo = TH1D("ok", "ok", 1200, 0, 1200)
# for event in tree:
#     if (event.z_vertex == -9999999.999999998):
#         continue
#     if 45 <= event.z_vertex <= 695:
#         coef=193.9753
#     else:
#         coef = 359.46
#     histo.Fill(event.RBS_res, 1/event.crosssection_rbs*coef)
#     #print(1/event.crosssection_rbs*coef)
    
# DisplayTH1D(histo, ax = ax)

##############calib

fig, ax = plt.subplots()
rootfile = TFile("vide.root")
#histo = DisplayTH1D(rootfile.Get("STIM_Ekin"), ax = ax)


tree = rootfile.Get("Detectors")
histo1 = TH1D("STIM_Ekin", "STIM_Ekin", 1300, 0, 1300)
for event in tree:
    if (event.STIM_Ekin != 0):
        histo1.Fill(event.STIM_Ekin + np.random.normal(0, 6.7))
    
histo = DisplayTH1D(histo1, ax = ax)

data=[]
condition = False
with open("2023W43/STIM_trou_bis.mpa", "r") as file:
    for line in file:
        if  "[DATA14,1024 ]\n" in str(line) :
            condition = True
            continue
        if str(line) == "[DATA15,1024 ]\n":
            break
        if condition == True:
            data.append(int(line))

bin = data.index(max(data))
data = [max(histo[1])/max(data)*i for i in data]

x = np.linspace(1, 1024, 1024)
initial_guess = [8000, 550, 20]
popt, pcov = curve_fit(gaussian, x, data, p0=initial_guess)

popt1, pcov1 = curve_fit(gaussian, histo[0], histo[1], p0=[5000, 1200, 10])
#print(popt1)
#plt.plot(histo[0], gaussian(histo[0], *popt1))
x_calib = np.linspace(1, 1.0005*popt1[1]/popt[1]*1024, 1024)

resolution = 1.0005*popt1[1]/popt[1] * popt[2]
print("RESOLUTION : ", 1.0005*popt1[1]/popt[1] * popt[2])
plt.plot(x_calib, data)
plt.xlim(1100, 1300)
#plt.show()

#####CALIB STIM 
print("CALIBRATION STIM : ", 1.0005*popt1[1]/popt[1], "keV/channel")


##############data
# fig, ax = plt.subplots()
# rootfile = TFile("AuSi.root")
# histo = DisplayTH1D(rootfile.Get("RBS"), ax = ax)


# data=[]
# condition = False
# with open("2023W43/RBS_STDAuSi.mpa", "r") as file:
#     for line in file:
#         if  "[DATA15,1024 ]\n" in str(line) :
#             condition = True
#             continue
#         if str(line) == "[DATA16,1024 ]\n":
#             break
#         if condition == True:
#             data.append(int(line))

# data = [max(histo[1])/max(data)*i for i in data]

# x = np.linspace(1, 1.99*1024, 1024)
# plt.plot(x, data, color='red')
# plt.xlim(0, 1200)
# #plt.show()

# print("CALIBRATION RBS : ", 1.994, "keV/channel")


##############data
fig, ax = plt.subplots()
rootfile = TFile("../build/test.root")
# x = []
# y=[]
# tree = rootfile.Get("Detectors")
# histo1 = TH1D("RBS", "RBS", 1200, 0, 1200)
# for event in tree:
#     if (event.RBS_Ekin != 0):
#         x.append(event.RBS_Ekin)
#         y.append(event.RBS_weight)
# plt.scatter(x, y)
histo = DisplayTH1D(rootfile.Get("RBS_Ekin_weight"), ax = ax)


data=[]
condition = False
with open("2023W43/RBS_Pos4.mpa", "r") as file:
    for line in file:
        if  "[DATA15,1024 ]\n" in str(line) :
            condition = True
            continue
        if str(line) == "[DATA16,1024 ]\n":
            break
        if condition == True:
            data.append(int(line))

data = [max(histo[1])/max(data)*i for i in data]

x = np.linspace(1, 1.995*1024, 1024)
plt.plot(x, data, color='red')
plt.xlim(800, 1100)
plt.show()

