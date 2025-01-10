import numpy as np
import corner
import matplotlib.pyplot as plt
from ROOT import *
from root2mpl import *
import subprocess


# Define labels for parameters
labels = [r"$\epsilon_{\mathrm{Al}}$", r"$\epsilon_{mylar}$", r"$\epsilon_{\mathrm{Al}}$", r"$\epsilon_{mylar}$", r"$b$", r"$a$", r"$c$", r"$d$"]
labels = [r"$\epsilon_{\mathrm{Al}}$", r"$\epsilon_{mylar}$", r"$b$", r"$a$"]


true_values = [
    85,       # Al_thickness_thin
    525,      # Mylar_thikness_thin
    # 110,      # Al_thickness_thick
    # 6100,     # Mylar_thikness_thick
    21.3155, # offset_calib_1
    1.9697,   # coefficients_calib_1
    # 72.0385,  # offset_calib_3
    # 3.72882   # coefficients_calib_3
]

# Load saved MCMC data
try:
    chain = np.load("mcmc_chain_step100.npy")
    lnprob = np.load("mcmc_log_prob_step100.npy")
except FileNotFoundError:
    print("Error: chain.npy or lnprob.npy not found. Ensure they are in the current directory.")
    exit(1)

################ CLEANING ################
def cleaning(chain):
    """
    Remove identical consecutive steps in the MCMC chain for each walker.

    Parameters:
        chain (ndarray): The MCMC chain of shape (nsteps, ndim) or (nsteps, nwalkers, ndim).

    Returns:
        ndarray: Filtered chain with identical consecutive steps removed.
    """
    # nsteps, ndim = chain.shape  # nsteps: number of steps, ndim: number of parameters

    # Initialize a list to store the filtered chain
    filtered_chain = [chain[0]]  # Always keep the first step
    log = []

    # Iterate through the rest of the steps
    for i in range(1, len(chain)):
            if np.all(chain[i][3] > 1.9):
                if np.all(chain[i][0] > 5) and np.all(chain[i][0] < 150):
                    filtered_chain.append(chain[i])

        # if not np.all(i % 2):
        #     filtered_chain.append(chain[i])
    # Convert the list back to a numpy array
    filtered_chain = np.array(filtered_chain)

    return filtered_chain

################################################################
# chain = chain[:: , :]
# lnprob = lnprob[200::, :]
# Determine dimensions
# nsteps, ndim = chain.shape
# print(f"Loaded chain with shape: {chain.shape}")
ndim=4
# Discard burn-in and flatten the chain
burnin = 10000
flat_samples = chain.reshape(-1, ndim)

# Plot parameter traces
fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
for i in range(ndim):
    ax = axes[i]
    ax.plot(chain[:, i], "k", alpha=0.3)
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
axes[-1].set_xlabel("step number")
fig.suptitle("Parameter Traces", fontsize=16)
fig.savefig("mcmc_trace.png")
plt.close(fig)

min_values = np.min(flat_samples, axis=0)
max_values = np.max(flat_samples, axis=0)
bin_values = [5, 5, 5, 5, 1, 0.002, 1, 0.002]
bin_values = [5, 5, 1, 0.002]

custom_bins = [
    28, 
    44,
    # 100,
    # 100,
    100,
    100,
    # 100,
    # 100
]

for i in range(ndim):
    custom_bins[i] = round((max_values[i]-min_values[i])/bin_values[i])


flat_samples_cleaned = cleaning(flat_samples)
param_means = np.mean(flat_samples_cleaned, axis=0)
param_stds = np.std(flat_samples_cleaned, axis=0)
# true_values = chain.reshape(-1, chain.shape[-1])[np.argmin(-2 * lnprob.flatten())]


unit = [
    " [nm]",
    " [nm]"
    # " [nm]",
    # " [nm]",
    " [keV]",
    " [keV/CH]"
    # " [keV]",
    # " [keV/CH]"
]
    
scatter_kws = {"s": 50, "alpha": 0.3, "color": "blue"}  # Spot size, transparency, and color


# Generate corner plot
fig = corner.corner(flat_samples_cleaned, labels=labels, quantiles=[0.16, 0.5, 0.84], 
                    show_titles=True, title_fmt=".3f", bins=custom_bins, levels=[0.68, 0.95], smooth=2.0,   
                    plot_datapoints=False,
                    truths=true_values,
)
fig.suptitle("Corner Plot", fontsize=16)
axes = np.array(fig.axes).reshape((ndim, ndim))
# for i in range(ndim):
    # Update x-axis label
    # axes[-1, i].set_xlabel(labels[i] + unit[i])
    # # Update y-axis label
    # if i > 0:
    #     axes[i, 0].set_ylabel(labels[i] + unit[i])
fig.savefig("corner_plot.png")
plt.show()

#Generate spectrum plot
#get quanils or each paramater
# quantiles = np.percentile(flat_samples, [16, 84], axis=0)
# print(quantiles)
# # # Generate a spectrum plot for optimum values
# fig, ax = plt.subplots(2, 1)

# command = f"./RBS_All_forMCMC {true_values[0]} {true_values[1]} 100 6100 {true_values[2]} {true_values[3]} 70 3.76 0"

# result = subprocess.run(
#         command,
#         shell=True,  
#         # stdout=subprocess.PIPE,
#         stderr=subprocess.PIPE, 
#         # universal_newlines=True,
#     )

# file = TFile(f"RBS_Results.root", "READ")
# hist=[]
# # EXPERIMENTAL 
# hist.append(file.Get("Exp_Hist_calib_A"))
# DisplayTH1D(hist[-1], ax=ax[0], color="black", label="Experimental")
# hist.append(file.Get("Exp_Hist_calib_B"))
# DisplayTH1D(hist[-1], ax=ax[1], color="black", label="Experimental")

# # SIMULATED OPTIMUM
# hist.append(file.Get(f"Sim_Hist_conv_1.2_{true_values[0]}.0_{true_values[1]}.0_A_-2_4keV"))
# DisplayTH1D(hist[-1], ax=ax[0], color="red", label="Simulated")
# hist.append(file.Get(f"Sim_Hist_conv_1.2_{true_values[0]}.0_{true_values[1]}.0_B_5.25_4keV"))
# DisplayTH1D(hist[-1], ax=ax[1], color="red", label="Simulated")


# bin_centers = []
# for i in range(hist[-1].GetNbinsX()):
#     bin_centers.append(hist[0].GetBinCenter(i+1))
# file.Close()

# for i in range(len(true_values)):
#     for j in range(len(true_values)):
#         if i == j:
#             continue
#         for add in quantiles:
#             current_values = true_values.copy()
#             current_values[j] = add[j]
#             command = f"./RBS_All_forMCMC {current_values[0]} {current_values[1]} 100 6100 {current_values[2]} {current_values[3]} 70 3.76 0"
#             subprocess.run(
#                 command,
#                 shell=True,  
#                 stderr=subprocess.PIPE, 
#             )
#             file = TFile(f"RBS_Results.root", "READ")
#             for key in file.GetListOfKeys():
#                     print(key.GetName())
#             if j < 2:
#                 hist.append(file.Get(("Sim_Hist_conv_1.2_{:.1f}_{:.1f}_A_-2_4keV").format(current_values[0], current_values[1])))
#                 DisplayTH1D(hist[-1], ax=ax[0], color="red")
#                 hist.append(file.Get(("Sim_Hist_conv_1.2_{:.1f}_{:.1f}_B_5.25_4keV").format(current_values[0], current_values[1])))
#                 DisplayTH1D(hist[-1], ax=ax[1], color="red")
#             else:
#                 hist.append(file.Get("Exp_Hist_calib_A"))
#                 DisplayTH1D(hist[-1], ax=ax[0], color="gray")
#                 hist.append(file.Get("Exp_Hist_calib_B"))
#                 DisplayTH1D(hist[-1], ax=ax[1], color="gray")
#             file.Close()


flat_chain = chain.reshape(-1, chain.shape[-1])
flat_log_prob = lnprob.flatten()

# Compute chi-squared values
chi2 = -2 * flat_log_prob

# Find the minimum chi2 value and the threshold for the 68% credible region
chi2_min = np.min(chi2)
print(f"Minimum chi2: {chi2_min} or {chi2_min/len(flat_chain[0])} per degree of freedom")
print(f"Minimum chi2 parameters: {flat_chain[np.argmin(-2 * flat_log_prob)]}")
chi2_threshold = chi2_min + 1  # 68% corresponds to Δχ² ≤ 1

# Filter the chain for parameter sets within the credible region
credible_region_indices = chi2 <= chi2_threshold
credible_region_params = flat_chain[credible_region_indices]

list_th=[]
list_calib=[]
counter=0
hist_data = [[], []]
hist_sim = [[], []]
for params in credible_region_params:
    th=False
    calib=False
    if [params[0], params[1]] not in list_th:
        th= True
        list_th.append([params[0], params[1]])
    if [params[2], params[3]] not in list_calib:
        calib=True
        list_calib.append([params[2], params[3]])

    if not th and not calib:
        continue

    command = f"./RBS_All_forMCMC {params[0]} {params[1]} 100 6100 {params[2]} {params[3]} 70 3.76 0"
    subprocess.run(
        command,
        shell=True,  
        stderr=subprocess.PIPE, 
    )

    file = TFile(f"RBS_Results.root", "READ")

    if calib:
        hist.append(file.Get("Exp_Hist_calib_A"))
        hist_data[0].append(DisplayTH1D(hist[-1], ax=ax[0], color="gray", visible=False))
        hist.append(file.Get("Exp_Hist_calib_B"))
        hist_data[1].append(DisplayTH1D(hist[-1], ax=ax[1], color="gray", visible=False))
    if th:
        hist.append(file.Get("Sim_Hist_conv_1.2_{:.1f}_{:.1f}_A_-2_4keV".format(params[0], params[1])))
        hist_sim[0].append(DisplayTH1D(hist[-1], ax=ax[0], color="red", visible=False))
        hist.append(file.Get("Sim_Hist_conv_1.2_{:.1f}_{:.1f}_B_5.25_4keV".format(params[0], params[1])))
        hist_sim[1].append(DisplayTH1D(hist[-1], ax=ax[1], color="red", visible=False))

    file.Close()
    counter+=1

hist_data[0] = [line.get_ydata() for line in hist_data[0]]
hist_data[1] = [line.get_ydata() for line in hist_data[1]]
hist_sim[0] = [line.get_ydata() for line in hist_sim[0]]
hist_sim[1] = [line.get_ydata() for line in hist_sim[1]]
hist_data_array = [np.array(d) for d in hist_data]
hist_sim_array = [np.array(d) for d in hist_sim]
hist_data_min = [np.min(d, axis=0) for d in hist_data_array]
hist_data_max = [np.max(d, axis=0) for d in hist_data_array]
hist_sim_min = [np.min(d, axis=0) for d in hist_sim_array]
hist_sim_max = [np.max(d, axis=0) for d in hist_sim_array]

alpha=0.2
ax[0].fill_between(bin_centers, hist_data_min[0], hist_data_max[0], color="black", alpha=alpha)
ax[0].fill_between(bin_centers, hist_sim_min[0], hist_sim_max[0], color="red", alpha=alpha)
ax[1].fill_between(bin_centers, hist_data_min[1], hist_data_max[1], color="black", alpha=alpha)
ax[1].fill_between(bin_centers, hist_sim_min[1], hist_sim_max[1], color="red", alpha=alpha)



plt.show()
plt.savefig("spectrum_plot.png")

# Print estimated parameter values
param_means = np.mean(flat_samples, axis=0)
param_stds = np.std(flat_samples, axis=0)
print("Estimated parameters:")
for i, (mean, std) in enumerate(zip(param_means, param_stds)):
    print(f"Parameter {labels[i]}: {mean:.4f} ± {std:.4f}")


