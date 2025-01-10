import numpy as np
import corner
import matplotlib.pyplot as plt


def remove_identical_neighbors(chain):
    """
    Remove identical consecutive steps in the MCMC chain for each walker.

    Parameters:
        chain (ndarray): The MCMC chain of shape (nsteps, ndim) or (nsteps, nwalkers, ndim).

    Returns:
        ndarray: Filtered chain with identical consecutive steps removed.
    """
    nsteps, ndim = chain.shape  # nsteps: number of steps, ndim: number of parameters

    # Initialize a list to store the filtered chain
    filtered_chain = [chain[0]]  # Always keep the first step

    # Iterate through the rest of the steps
    for i in range(1, nsteps):
        # Check if the current step is different from the previous one
        # if not np.all(chain[i][2] < -20) and np.all(chain[i][2] < 70):
        #     if not np.all(chain[i][3] < 1.84):
        #         filtered_chain.append(chain[i])

        if not np.all(i % 2):
            filtered_chain.append(chain[i])
    # Convert the list back to a numpy array
    filtered_chain = np.array(filtered_chain)

    return filtered_chain

# Example usage:
# Load the MCMC chains
chain = np.load("chain.npy")  # Assuming shape is (nsteps, ndim)
print(f"Loaded chain shape: {chain.shape}")
chain = chain[10000::, :]
# Remove identical consecutive steps
filtered_chain = remove_identical_neighbors(chain)

# Print filtered chain shape
print(f"Filtered chain shape: {filtered_chain.shape}")


# Define labels for parameters
labels = ["al_thin", "mylar_thin", "off1", "coeff1"]

# Define true parameter values for reference in the corner plot
#corner plot
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

ndim=4
flat_samples = filtered_chain.reshape(-1, ndim)

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

#mean
mean_values = np.mean(flat_samples, axis=0)

#get min and max of each parameters

min_values = np.min(flat_samples, axis=0)
max_values = np.max(flat_samples, axis=0)
bin_values = [5, 5, 1, 0.002]


custom_bins = [28, 44, 100, 100]
for i in range(ndim):
    custom_bins[i] = round((max_values[i]-min_values[i])/bin_values[i])

# Generate corner plot
fig = corner.corner(flat_samples, labels=labels, truths=mean_values, quantiles=[0.16, 0.5, 0.84], show_titles=True, bins=custom_bins)
fig.suptitle("Corner Plot", fontsize=16)
fig.savefig("corner_plot.png")
plt.close(fig)

# Print estimated parameter values
param_means = np.mean(flat_samples, axis=0)
param_stds = np.std(flat_samples, axis=0)
print("Estimated parameters:")
for i, (mean, std) in enumerate(zip(param_means, param_stds)):
    print(f"Parameter {labels[i]}: {mean:.4f} ± {std:.4f}")