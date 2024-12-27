import numpy as np
import corner
import matplotlib.pyplot as plt

# Define labels for parameters
labels = ["al_thin", "mylar_thin", "off1", "coeff1"]

# Define true parameter values for reference in the corner plot
true_values = [
    73,       # Al_thickness_thin
    538.7,      # Mylar_thikness_thin
    # 110,      # Al_thickness_thick
    # 6100,     # Mylar_thikness_thick
    # 22.25, # offset_calib_1
    # 1.9609,   # coefficients_calib_1
    # 72.0385,  # offset_calib_3
    # 3.72882   # coefficients_calib_3
]

# Load saved MCMC data
try:
    chain = np.load("chain.npy")
    lnprob = np.load("log_prob.npy")
except FileNotFoundError:
    print("Error: chain.npy or lnprob.npy not found. Ensure they are in the current directory.")
    exit(1)

chain = chain[::, :]
# Determine dimensions
nsteps, ndim = chain.shape
print(f"Loaded chain with shape: {chain.shape}")

# Discard burn-in and flatten the chain
burnin = 200
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

# Generate corner plot
fig = corner.corner(flat_samples, labels=labels, truths=true_values, bins=40)
fig.suptitle("Corner Plot", fontsize=16)
fig.savefig("corner_plot.png")
plt.close(fig)

# Print estimated parameter values
param_means = np.mean(flat_samples, axis=0)
param_stds = np.std(flat_samples, axis=0)
print("Estimated parameters:")
for i, (mean, std) in enumerate(zip(param_means, param_stds)):
    print(f"Parameter {labels[i]}: {mean:.4f} ± {std:.4f}")
