import emcee
import numpy as np
import subprocess
import corner
import sys
from multiprocessing import Pool
import threading
import matplotlib.pyplot as plt
import os


def run_root_macro(params, thread):
    params[0] = round((params[0]+2.5)/5)*5
    params[1] = round((params[1]+2.5)/5)*5
    # params[2] = round((params[2]+2.5)/5)*5
    # params[3] = round((params[3]+2.5)/5)*5

    ####
    # print(params)
    #add value to numpy array
    params = np.append(params, 0)
    params = np.append(params, 0)
    params = np.append(params, 0)
    params = np.append(params, 0)
    params = np.append(params, offset_calib_3)
    params = np.append(params, coefficients_calib_3)
    # params[2] = 100
    # params[3] = 6000
    ####
    param_str = " ".join(map(str, params))

    command = f"./RBS_All_forMCMC {param_str} {thread}"

    result = subprocess.run(
        command,
        shell=True,  
        # stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, 
        # universal_newlines=True,
    )
    try:
        with open(f"tmp/{thread}.txt", "r") as f:
            lines = f.readlines()
            chi2 = float(lines[-1])

        # reading cout from command line
        # chi2 = float(result.stdout.split("\n")[-2])
        
    except (IndexError, ValueError):
        raise RuntimeError(f"Failed to parse chi2 from ROOT output: {result.stdout}")
    
    return chi2

def log_prob(params):
    # Define parameter bounds 
    if (params[0] > 150 or params[0] < 50):
        return -np.inf
    if (params[1] > 650 or params[1] < 400):
        return -np.inf
    # if (params[2] > 200 or params[2] < 50):
    #     return -np.inf
    # if (params[3] > 6500 or params[3] < 5500):
    #     return -np.inf
    # if (params[2] < -200 or params[2] > 200):
    #     return -np.inf
    # if (params[3] < 1.7 or params[3] > 2.2):
    #     return -np.inf
    
    process_id = os.getpid()

    chi2 = run_root_macro(params, process_id)

    if (chi2 == 0):
        return -np.inf  
    
    return -0.5 * chi2

ndim = 2  # Number of parameters
nwalkers = 10*ndim  # Number of walkers

Al_thickness_thin = 85
Mylar_thikness_thin = 525

Al_thickness_thick = 110
Mylar_thikness_thick = 6100

offset_calib_1 = 21.3155
coefficients_calib_1 = 1.9697

offset_calib_3 = 72.0385
coefficients_calib_3 = 3.72882


thickness = 5
offsetcalib = 2
coefficientscalib = 0.01

initial_position = [
    [   Al_thickness_thin + np.random.normal(-thickness, thickness),
        Mylar_thikness_thin + np.random.normal(-thickness, thickness)
        # Al_thickness_thick + np.random.normal(-thickness, thickness),
        # Mylar_thikness_thick + np.random.normal(-thickness, thickness),
        # offset_calib_1 + np.random.normal(-offsetcalib, offsetcalib),
        # coefficients_calib_1 + np.random.normal(-coefficientscalib, coefficientscalib),
        # offset_calib_3 + np.random.normal(-offsetcalib, offsetcalib),
        # coefficients_calib_3 + np.random.normal(-coefficientscalib, coefficientscalib)
        ]
    for _ in range(nwalkers)
     ]
    


# Run the MCMC sampling
nsteps = 10000
burnin = 200

with Pool(processes=8) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool)
    sampler.run_mcmc(initial_position, nsteps, progress=True)

# sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob)
# sampler.run_mcmc(initial_position, nsteps, progress=True)

# Analyze the results
samples = sampler.get_chain(discard=burnin, flat=True)  # Discard burn-in and thin samples
param_means = np.mean(samples, axis=0)
param_stds = np.std(samples, axis=0)

## saving
np.save("chain.npy", sampler.get_chain(flat=True))
np.save("log_prob.npy", sampler.get_log_prob())


print("Estimated parameters:")
for i, (mean, std) in enumerate(zip(param_means, param_stds)):
    print(f"Parameter {i}: {mean:.4f} ± {std:.4f}")

fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["al_thin", "mylar_thin", "off1", "coeff1"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number")
fig.savefig("mcmc.png")

fig = corner.corner(
    samples, labels=labels, truths=[Al_thickness_thin, Mylar_thikness_thin, Al_thickness_thick, Mylar_thikness_thick, offset_calib_1, coefficients_calib_1, offset_calib_3, coefficients_calib_3]
)

fig.savefig("corner.png")

