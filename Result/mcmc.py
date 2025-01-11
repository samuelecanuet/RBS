import emcee
import numpy as np
import subprocess
import corner
import sys
from multiprocessing import Pool
import threading
import matplotlib.pyplot as plt
import os
import StretchMoveInteger
from root2mpl import *
from ConfidenceLevel_Dim import *
from ROOT import *
from tqdm import tqdm
from colorama import Fore, Style

def run_root_macro(params, thread):
    counter= 0
    for i in range(len(parameters_full)):
        if parameters_full[i] != -1111:
            parameters_full[i] = params[counter]
            counter+=1


    param_str = " ".join(map(str, parameters_full))

    command = f"./RBS_All_forMCMC {param_str} {thread}"
    result = subprocess.run(
        command,
        shell=True,  
        stderr=subprocess.PIPE, 
    )
    try:
        with open(f"tmp/{thread}.txt", "r") as f:
            lines = f.readlines()
            chi2 = float(lines[-1])
        
    except (IndexError, ValueError):
        raise RuntimeError(f"Failed to parse chi2 from ROOT output: {result.stdout}")
    
    return chi2

def log_prob(params):
    
    #define limits
    if np.any((params > parlimits[:, 1]) | (params < parlimits[:, 0])):
        return -np.inf
    
    process_id = os.getpid()

    chi2 = run_root_macro(params, process_id)

    if (chi2 == 0):
        return -np.inf  
    
    return -0.5 * chi2

def save_progress(sampler, step, save_interval=100):
    """Save the chain and log-probability during MCMC every save_interval steps."""
    if step % save_interval == 0:
        np.save(f"chain_step{step}.npy", sampler.get_chain())
        np.save(f"log_prob_step{step}.npy", sampler.get_log_prob())
        if step-save_interval > 0:
            os.remove(f"chain_step{step-save_interval}.npy")
            os.remove(f"log_prob_step{step-save_interval}.npy")

########################## PARAMETERS ##########################
type = ["THIN"]
energy = [1.2]

nsteps = 100
burnin = 100

ndim = len(type) * 2 + len(energy) * 2
nwalkers = 10*ndim

Al_thickness_thin = 85
Al_thickness_thin_lim = [5, 150]
Mylar_thikness_thin = 525
Mylar_thikness_thin_lim = [400, 650]

Al_thickness_thick = 110
Al_thickness_thick_lim = [50, 200]
Mylar_thikness_thick = 6100
Mylar_thikness_thick_lim = [5950, 6150]

offset_calib_1 = 21.3155
offset_calib_1_lim = [-25, 75]
coefficients_calib_1 = 1.9697
coefficients_calib_1_lim = [1.85, 2.08]

offset_calib_3 = 72.0385
offset_calib_3_lim = [30, 100]
coefficients_calib_3 = 3.72882
coefficients_calib_3_lim = [3.5, 4]

thickness = 20
offsetcalib = 20
coefficientscalib = 0.06

########## SELECTED PARAMETERS ##########
parameters = []
parameters_full = []
parlimits = []
labels = []
units = []
integer_params=[]
bin_width = []
bins =[0]*ndim
initial_position = [[] for _ in range(nwalkers)]
if "THIN" in type:
    parameters.append(Al_thickness_thin)
    parameters_full.append(Al_thickness_thin)
    parlimits.append(Al_thickness_thin_lim)
    bin_width.append(5)
    labels.append(r"$\epsilon_{Al}$")
    units.append(" [nm]")
    integer_params.append(0)
    for i in range(nwalkers): initial_position[i].append(round((Al_thickness_thin + np.random.normal(0, thickness))/5)*5)
    parameters.append(Mylar_thikness_thin)
    parameters_full.append(Mylar_thikness_thin)
    parlimits.append(Mylar_thikness_thin_lim)
    bin_width.append(5)
    labels.append(r"$\epsilon_{Mylar}$")
    units.append(" [nm]")
    integer_params.append(1)
    for i in range(nwalkers): initial_position[i].append(round((Mylar_thikness_thin + np.random.normal(0, thickness))/5)*5)
else:
    parameters_full.append(-1111)
    parameters_full.append(-1111)
if "THICK" in type:
    parameters.append(Al_thickness_thick)
    parameters_full.append(Al_thickness_thick)
    parlimits.append(Al_thickness_thick_lim)
    bin_width.append(5)
    labels.append(r"$\epsilon_{Al}$")
    units.append(" [nm]")
    if len(type) == 1: integer_params.append(0)
    else: integer_params.append(2)
    for i in range(nwalkers): initial_position[i].append(round((Al_thickness_thick + np.random.normal(0, thickness))/5)*5)
    parameters.append(Mylar_thikness_thick)
    parameters_full.append(Mylar_thikness_thick)
    parlimits.append(Mylar_thikness_thick_lim)
    bin_width.append(5)
    labels.append(r"$\epsilon_{Mylar}$")
    units.append(" [nm]")
    if len(type) == 1: integer_params.append(1)
    else: integer_params.append(3)
    for i in range(nwalkers): initial_position[i].append(round((Mylar_thikness_thick + np.random.normal(0, thickness))/5)*5)
else:
    parameters_full.append(-1111)
    parameters_full.append(-1111)
for en in energy:
    if en == 1.2:
        parameters.append(offset_calib_1)
        parameters_full.append(offset_calib_1)
        parlimits.append(offset_calib_1_lim)
        bin_width.append(1)
        labels.append(r"$b_{1.2}$")
        units.append(" [keV]")
        for i in range(nwalkers): initial_position[i].append(offset_calib_1 + np.random.normal(0, offsetcalib))
        parameters.append(coefficients_calib_1)
        parameters_full.append(coefficients_calib_1)
        parlimits.append(coefficients_calib_1_lim)
        bin_width.append(0.002)
        labels.append(r"$a_{1.2}$")
        units.append(" [keV/CH]")
        if len(energy) == 1: 
            parameters_full.append(-1111)
            parameters_full.append(-1111)
        for i in range(nwalkers): initial_position[i].append(coefficients_calib_1 + np.random.normal(0, coefficientscalib))
    if en == 3.0:
        if len(energy) == 1: 
            parameters_full.append(-1111)
            parameters_full.append(-1111)
        parameters.append(offset_calib_3)
        parameters_full.append(offset_calib_3)
        parlimits.append(offset_calib_3_lim)
        bin_width.append(1)
        labels.append(r"$b_{3.0}$")
        units.append(" [keV]")
        for i in range(nwalkers): initial_position[i].append(offset_calib_3 + np.random.normal(0, offsetcalib))
        parameters.append(coefficients_calib_3)
        parameters_full.append(coefficients_calib_3)
        parlimits.append(coefficients_calib_3_lim)
        bin_width.append(0.002)
        labels.append(r"$a_{3.0}$")
        units.append(" [keV/CH]")
        for i in range(nwalkers): initial_position[i].append(coefficients_calib_3 + np.random.normal(0, coefficientscalib))

parameters_full=np.array(parameters_full)
#########################################
################################################################


########################## STARTING ##########################
if (len(sys.argv) > 1):
    if (sys.argv[1] == "run"):
        stretch_move = StretchMoveInteger.StretchMoveInteger(a=2.0, integer_params=integer_params)
        with Pool(processes=8) as pool:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool, moves=[stretch_move])

            # Run MCMC with periodic saving
            for step, _ in enumerate(sampler.sample(initial_position, iterations=nsteps, progress=True)):
                save_progress(sampler, step + 1)
    
        np.save("chain.npy", sampler.get_chain(flat=True))
        np.save("log_prob.npy", sampler.get_log_prob())

try:
    chain = np.load("chain_thin_1p2_4par_40w_binned_wxs.npy")
    # chain = np.load("chain_all_8par_80w.npy")
    # chain = chain.reshape(-1, nwalkers, ndim)
    lnprob = np.load("log_prob_thin_1p2_4par_40w_binned_wxs.npy")
    # lnprob = np.load("log_prob_all_8par_80w.npy")
except FileNotFoundError:
    print("Error: chain.npy or lnprob.npy not found. Ensure they are in the current directory.")
    exit(1)

##################################################################

########################## PLOTTING ##############################
# Analyze the results
# chain = chain[burnin:, :]
# 
# lnprob = lnprob.reshape(-1, nwalkers)
lnprob = lnprob.flatten()
# lnprob = lnprob[burnin:]
flat_chain = chain.reshape(-1, ndim)

## BINING
min_values = np.min(flat_chain, axis=0)
max_values = np.max(flat_chain, axis=0)
for i in range(ndim): bins[i] = round((max_values[i]-min_values[i])/bin_width[i])


print("---- MEAN parameters:")
param_means = np.mean(flat_chain, axis=0)
param_stds = np.std(flat_chain, axis=0)
for i, (mean, std) in enumerate(zip(param_means, param_stds)):
    print(f"Parameter {i}: {mean:.4f} ± {std:.4f}")

########################### WALKERS PLOT ##########################
fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
for i in range(ndim):
    ax = axes[i]
    ax.plot(chain[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(chain))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number")
fig.savefig("Walkers.png")
##################################################################

########################### CORNER PLOT ##########################
minimum = flat_chain[np.argmin(-2 * lnprob)]

fig = corner.corner(flat_chain, labels=labels, quantiles=[0.16, 0.5, 0.84], 
                    show_titles=True, title_fmt=".4f", levels=[0.68, 0.95], smooth=2.0,
                    plot_datapoints=False,
                    truths=minimum, bins=bins
)
fig.suptitle("Corner Plot", fontsize=16)
axes = np.array(fig.axes).reshape((ndim, ndim))
for i in range(ndim):
    # Update x-axis label
    axes[-1, i].set_xlabel(labels[i] + units[i])
    # Update y-axis label
    if i > 0:
        axes[i, 0].set_ylabel(labels[i] + units[i])

fig.savefig("Corner.png")


print("---- BEST parameters - GLOBAL ERROR : CHI2+9.30: ", np.min(-2*lnprob))
indices_chi2p1 = -2 * lnprob.flatten() <= np.min(-2 * lnprob) + find_delta(sigma=1, nu=ndim)
params_chi2p1 = flat_chain[indices_chi2p1]

params_error = [
    [
abs(np.min(params_chi2p1[:, i]))-minimum[i], abs(np.max(params_chi2p1[:, i])-minimum[i])] for i in range(ndim)
]
for i, min, err, u in zip(range(ndim), minimum, params_error, units):
    print(f"Parameter {i}: {min:.4f} ± {max(err):.4f} {u}")

##################################################################

########################### RESTRICTED CORNER PLOT ##########################
## restricting paramaters
def cleaning(chain_, lnprob_, restricting):
    filtered_chain_ = []
    filtered_ln_prob_ = []

    # Iterate through the rest of the steps
    for step in range(nwalkers, len(lnprob_)-nwalkers):
        taking = True
        for par, ranges in enumerate(restricting):
            if (ranges[0] == 0) and (ranges[1] == 0):
                continue
            elif (ranges[0] < chain[step][par] < ranges[1]):
                continue
            else:
                taking = False
        if (chain_[step-nwalkers] == chain_[step]).all() or (chain_[step+nwalkers] == chain_[step]).all():
            #replace by nan array
            chain_[step] = np.full(ndim, np.nan)
            lnprob_[step] = np.nan
        if taking:
            filtered_chain_.append(chain_[step])
            filtered_ln_prob_.append(lnprob_[step])
    filtered_chain_ = np.array(filtered_chain_)
    filtered_ln_prob_ = np.array(filtered_ln_prob_)

    return filtered_chain_, filtered_ln_prob_


restricting = [
]

filtered_chain, filtered_ln_prob = cleaning(flat_chain, lnprob, restricting)

## BINING
#delete nan values from filtered chain
filtered_chain.reshape(-1, ndim)
min_values = np.nanmin(filtered_chain, axis=0)
max_values = np.nanmax(filtered_chain, axis=0)

for i in range(ndim): bins[i] = round((max_values[i]-min_values[i])/bin_width[i])

##
filtered_chain_wonan = filtered_chain[~np.isnan(filtered_chain).any(axis=1)]
filtered_ln_prob_wonan = filtered_ln_prob[~np.isnan(filtered_chain).any(axis=1)]

minimum = filtered_chain_wonan[np.argmin(-2 * filtered_ln_prob_wonan)]

if (ndim == 8):
    minimum = np.array([85, 525, 110, 6100, 22.0887, 1.9681, 72.0385, 3.72882])
else:
    minimum = np.array([85, 525, 22.0887, 1.9681])

fig1 = corner.corner(filtered_chain_wonan, labels=labels, quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_fmt=".4f", levels=[0.68, 0.95], smooth=2.0,
                    plot_datapoints=False,
                    truths=minimum, 
                    bins=bins
)

fig1.suptitle("Corner Plot", fontsize=16)
axes = np.array(fig1.axes).reshape((ndim, ndim))
for i in range(ndim):
    # Update x-axis label
    axes[-1, i].set_xlabel(labels[i] + units[i])
    # Update y-axis label
    if i > 0:
        axes[i, 0].set_ylabel(labels[i] + units[i])

fig1.savefig("Corner_restricted.png")

print("---- BEST parameters - LOCAL ERROR ( CHI2 + 4.30 for each histograms ): ", np.min(-2*filtered_ln_prob_wonan))
### Chi² + 1
ntype = len(integer_params)/2
nenergy = ( len(parameters) - len(integer_params) ) /2
if ntype == 2 and nenergy == 2:
    n=4
else:
    n = max(ntype, nenergy)


indices_chi2p1 = -2 * filtered_ln_prob_wonan <= np.min(-2 * filtered_ln_prob_wonan) + find_delta(sigma=0.25, nu=ndim)
params_chi2p1 = filtered_chain_wonan[indices_chi2p1]

params_error = [
    [
abs(np.min(params_chi2p1[:, i])-minimum[i]), abs(np.max(params_chi2p1[:, i])-minimum[i])] for i in range(ndim)
]
for i, min, err, u in zip(range(ndim), minimum, params_error, units):
    print(f"Parameter {i}: {min:.4f} ± {max(err):.4f} {u}")

filtered_chain = filtered_chain.reshape(-1, nwalkers, ndim)
fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
for i in range(ndim):
    ax = axes[i]
    ax.plot(filtered_chain[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(filtered_chain))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number")
fig.savefig("Walkers_restricted.png")

####################################################################
########################### SPECTRUM PLOT ##########################
## PLOTTING MINIUM FROM MANUAL OR MCMC
#######
## PLOTTING A FILL BETWEEN WITH CHI2+1 RESPECTED FOR EACH HISTOGRAMS (NOT GLOBAL)
#######


## Prepare number of paramater and so one for the plot
ntype = len(integer_params)/2
nenergy = ( len(parameters) - len(integer_params) ) /2

if ntype == 1 and nenergy == 1:
    ny=1
    if (parameters_full[4] == -1111):
        Energies=[3.0]
    else:
        Energies=[1.2]
    if (parameters_full[0] == -1111):
        Types=["THICK"]
    else:
        Types=["THIN"]
if ntype == 2 and nenergy == 1:
    ny=2
    Types=["THIN", "THICK"]
    if (parameters_full[0] == -1111):
        Energies=[3.0]
if ntype == 1 and nenergy == 2:
    ny=2
    Energies=[1.2, 3.0]
    if (parameters_full[0] == -1111):
        Types=["THICK"]
if ntype == 2 and nenergy == 2:
    ny=4
    Types=["THIN", "THICK"]
    Energies=[1.2, 3.0]

#windows for plot
Exp_hist = []
Sim_hist = []
windows={}
windows["THIN"] = {"1.2": (800, 1100), "3.0": (2100, 2800)}
windows["THICK"] = {"1.2": (200, 1200), "3.0": (1900, 2700)}


# For best set of parameters
counter= 0
for i in range(len(parameters_full)):
    if parameters_full[i] != -1111:
        parameters_full[i] = minimum[counter]
        counter+=1
command = f"./RBS_All_forMCMC {' '.join(map(str, parameters_full))} 0"
parameters_min = parameters_full.copy()
result = subprocess.run(command, shell=True, stderr=subprocess.PIPE)
## read file for chi2
try:
    with open(f"tmp/0.txt", "r") as f:
        lines = f.readlines()
        chi2_mins = np.array([float(line) for line in lines])
except (IndexError, ValueError):
    raise RuntimeError(f"Failed to parse chi2 from ROOT output: {result.stdout}")
file = TFile("RBS_Results.root", "READ")

## do not open all the previous fig
plt.close("all")


# Starting plotting the spectrum with minimum
fig, ax = plt.subplots(ny, 2, figsize=(12, 3*ny))
fig.tight_layout(pad=2.0)
ax = ax.reshape(ny, -1)
c=[]
bin_centers = {}
bin_centers["1.2"] = []
bin_centers["3.0"] = []
counter = -1
for type in Types:
    for energy in Energies:
        counter+=1
        for i, face in enumerate(["A", "B"]):
            if i == 1 : ylabel = "" 
            else: ylabel = "Counts/keV"
            if counter==ny-1: xlabel="Energy [keV]"
            else: xlabel= ""

            ##Experiment
            Exp_hist.append(file.Get(f"Exp_Hist_{type}_{energy}_{face}"))
            c.append(DisplayTH1D(Exp_hist[-1], ax=ax[counter, i], xlim = windows[type][str(energy)], lw=1, titlesize=10, xlabel=xlabel, ylabel=ylabel, labelsize=15, ticksize=10, title = f"{type} {energy} MeV {face}", label="Experiment"))
            
            ##Simulated, 3.0
            Sim_hist.append(file.Get(f"Sim_Hist_conv_{type}_{energy}_{face}_4keV"))
            c.append(DisplayTH1D(Sim_hist[-1], ax=ax[counter, i], xlim = windows[type][str(energy)], lw=1, titlesize=10, xlabel=xlabel, ylabel=ylabel, labelsize=15, ticksize=10, color="red", title = f"{type} {energy} MeV {face}", label = "Simulated"))

            print(chi2_mins[counter*2+i])
            ax[counter, i].text(0.11, 0.9, r"$\chi^{2}_{\nu}$" + f" = {chi2_mins[counter*2+i]:.2f}", horizontalalignment='center', verticalalignment='center', color='red', transform=ax[counter, i].transAxes)

            bin_centers[str(energy)]=[]
            for i in range(Exp_hist[-1].GetNbinsX()):
                bin_centers[str(energy)].append(Exp_hist[-1].GetBinCenter(i+1))


### Lower and Upper error bar for the best parameters
## Chi² + ( 1 * N_hist) [ GLOBAL chi2 + 1 ]
# just to get reduced set of parameter more probably in the chi2+1 per histogram
credible_region = params_chi2p1

files=[]

hist_data = {}
hist_data["THIN"] = {}
hist_data["THICK"] = {}
hist_data["THIN"]["1.2"] = [[], []]
hist_data["THIN"]["3.0"] = [[], []]
hist_data["THICK"]["1.2"] = [[], []]
hist_data["THICK"]["3.0"] = [[], []]

hist_sim = {}
hist_sim["THIN"] = {}
hist_sim["THICK"] = {}
hist_sim["THIN"]["1.2"] = [[], []]
hist_sim["THIN"]["3.0"] = [[], []]
hist_sim["THICK"]["1.2"] = [[], []]
hist_sim["THICK"]["3.0"] = [[], []]
d=[]
paramaters_error = np.zeros(ndim)
print(len(credible_region))
for n, error in tqdm(enumerate(credible_region), desc="Loading"):
    if n%100 != 0:
        continue

    d.append(error[0])
    counter= 0
    for j in range(len(parameters_full)):
        if parameters_full[j] != -1111:
            parameters_full[j] = error[counter]
            counter+=1
    command = f"./RBS_All_forMCMC {' '.join(map(str, parameters_full))} 0"
    result = subprocess.run(command, shell=True, stderr=subprocess.PIPE)

    ## read file for chi2
    try:
        with open(f"tmp/0.txt", "r") as f:
            lines = f.readlines()
            chi2 = np.array([float(line) for line in lines])
        
    except (IndexError, ValueError):
        raise RuntimeError(f"Failed to parse chi2 from ROOT output: {result.stdout}")
    
    # selection chi2+1 per histogram
    if (chi2_mins + find_delta(1, 4) < chi2 ).any():
        continue

    ## saving error bar 
    for i in range(ndim):
        paramaters_error[i] = max(paramaters_error[i], abs(parameters_full[i]-minimum[i]))    

    files.append(TFile("RBS_Results.root", "READ"))
    counter = -1
    for type in Types:
        for energy in Energies:
            counter+=1
            for i, face in enumerate(["A", "B"]):
                ##Experiment
                Exp_hist.append(files[-1].Get(f"Exp_Hist_{type}_{energy}_{face}"))
                hist_data[type][str(energy)][i].append(DisplayTH1D(Exp_hist[-1], ax=ax[counter, i], color="gray", visible = False, xlim = windows[type][str(energy)], titlesize=10, xlabel="", ylabel="", labelsize=15, ticksize=10, title = f"{type} {energy} MeV {face}"))
            
                ##Simulated
                Sim_hist.append(files[-1].Get(f"Sim_Hist_conv_{type}_{energy}_{face}_4keV"))
                hist_sim[type][str(energy)][i].append(DisplayTH1D(Sim_hist[-1], ax=ax[counter, i], color="blue", visible=False, xlim = windows[type][str(energy)], titlesize=10, xlabel="", ylabel="", labelsize=15, ticksize=10, title = f"{type} {energy} MeV {face}"))

def moving_average(data, window_size=10):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')


# plotting fill between error bar
counter=-1
for type in Types:
    for energy in Energies:
        counter+=1
        hist_data[type][str(energy)][0] = [line.get_ydata() for line in hist_data[type][str(energy)][0]]
        hist_data[type][str(energy)][1] = [line.get_ydata() for line in hist_data[type][str(energy)][1]]
        hist_sim[type][str(energy)][0] = [line.get_ydata() for line in hist_sim[type][str(energy)][0]]
        hist_sim[type][str(energy)][1] = [line.get_ydata() for line in hist_sim[type][str(energy)][1]]
        hist_data_array = [np.array(d) for d in hist_data[type][str(energy)]]
        hist_sim_array = [np.array(d) for d in hist_sim[type][str(energy)]]
        hist_data_min = [np.min(d, axis=0) for d in hist_data_array]
        hist_data_max = [np.max(d, axis=0) for d in hist_data_array]
        hist_sim_min = [np.min(d, axis=0) for d in hist_sim_array]
        hist_sim_max = [np.max(d, axis=0) for d in hist_sim_array]

        alpha=0.2
        move = 4
        ax[counter, 0].fill_between(bin_centers[str(energy)][int(move/2-1):-int(move/2)], moving_average(hist_data_min[0], move), moving_average(hist_data_max[0], move), color="black", alpha=alpha, zorder=1, label=r"Experiment $\chi^{2}_{\nu}+1$")
        ax[counter, 0].fill_between(bin_centers[str(energy)][int(move/2-1):-int(move/2)], moving_average(hist_sim_min[0], move), moving_average(hist_sim_max[0], move), color="red", alpha=alpha, zorder=1, label=r"Simulated $\chi^{2}_{\nu}+1$")
        ax[counter, 1].fill_between(bin_centers[str(energy)][int(move/2-1):-int(move/2)], moving_average(hist_data_min[1], move), moving_average(hist_data_max[1], move), color="black", alpha=alpha, zorder=1, label=r"Experiment $\chi^{2}_{\nu}+1$")
        ax[counter, 1].fill_between(bin_centers[str(energy)][int(move/2-1):-int(move/2)], moving_average(hist_sim_min[1], move), moving_average(hist_sim_max[1], move), color="red", alpha=alpha, zorder=1, label=r"Simulated $\chi^{2}_{\nu}+1$")

ax[0, 0].legend(loc="upper right", fontsize = 10)
s = r"$\chi^{2}_{\nu}$ = " + "{:.2f}".format(np.min(-2*filtered_ln_prob_wonan)/2/ny)
print("---- BEST ANALYSIS parameters: ", "{:.2f}".format(np.sum(chi2_mins)))
for i, min, err, u in zip(range(ndim), minimum, paramaters_error, units):
    print(f"Parameter {i}: {min:.4f} ± {err:.4f} {u}")

fig.savefig("Spectrum.png", dpi=600)
plt.show()

