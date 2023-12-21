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



fig, ax = plt.subplots()
rootfile = TFile("test.root")
histo = DisplayTH1D(rootfile.Get("RBS2"), ax = ax)

##############calib
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

x = np.linspace(1, 1024, 1024)
initial_guess = [8000, 550, 20]
popt, pcov = curve_fit(gaussian, x, data, p0=initial_guess)
#925 for 1Mev
ratio = 1112/popt[1]
print(popt[2]*ratio)
print(ratio)

##############data
data=[]
condition = False
with open("2023W43/RBS_Pos3.mpa", "r") as file:
    for line in file:
        if  "[DATA15,1024 ]\n" in str(line) :
            condition = True
            continue
        if str(line) == "[DATA16,1024 ]\n":
            break
        if condition == True:
            data.append(2.26*int(line))

x = np.linspace(1, 1024*ratio, 1024)
plt.plot(x, data, color='red')
plt.xlim(800, 1100)
plt.show()

