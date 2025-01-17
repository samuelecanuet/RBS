from ROOT import TFile
from root2mpl import *
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import subprocess
import multiprocessing
from ConfidenceLevel_Dim import *
from scipy.optimize import curve_fit
from scipy.optimize import minimize

###correlation plot for 8 paramaters

ndim = 4

minimum_value = [85, 500, 22.0887, 1.9681]

range_value = [
    [50, 120],
    [400, 600],
    [18, 25],
    [1.96, 1.98],
]

step_value = [
    5,
    5,
    0.4,
    0.001,
]

step_value = [ step_value[i] for i in range(ndim) ]

unit_value = [
    " [nm]", 
    " [nm]",
    " [keV]",
    " [keV/CH]",
]

labels = [
    r"$\epsilon_{\mathrm{Al}}$",
    r"$\epsilon_{Mylar}$",
    r"$b_{1.2}$",
    r"$a_{1.2}$",
]

cov_matrix = []

fig, ax = plt.subplots(ndim, ndim, figsize = (10, 10))
# figg, axx = plt.subplots()
# fig.tight_layout(pad=1.0)
for i, irange in enumerate(range_value):
    for j, jrange in enumerate(range_value):
        if i < j:
            ax[i, j].axis('off')
            continue
        
        elif i > j:
            if j != 0:
                ax[i, j].set_yticklabels([])
            if i != ndim - 1:
                ax[i, j].set_xticklabels([])
            if j == 0:
                ax[i, j].set_ylabel(f"{labels[i]}" + unit_value[i], fontsize=10)
            if i == ndim - 1:
                ax[i, j].set_xlabel(f"{labels[j]}" + unit_value[j], fontsize=10)


## 8x8 list 
data_x = [
    [[] for _ in range(ndim)] for j in range(ndim)
]

data_y = [
    [[] for _ in range(ndim)] for j in range(ndim)
]

data_chi2 = [
    [[] for _ in range(ndim)] for j in range(ndim)
]

chi2_min = [
    3.36, 3.51, 2.14, 3.28, 2.22, 2.82, 1.25, 1.67
]

def ellipse(t, cx, cy, a, b, theta):
    """
    Parametric equation of an ellipse.
    t: Parametric angle (array-like)
    cx, cy: Center of the ellipse
    a, b: Semi-major and semi-minor axes
    theta: Rotation angle of the ellipse
    """

    if (a < b ):
        a, b = b, a
    x = cx + a * np.cos(t) * np.cos(theta) - b * np.sin(t) * np.sin(theta)
    y = cy + a * np.cos(t) * np.sin(theta) + b * np.sin(t) * np.cos(theta)
    return x, y

def fit_ellipse(Y, Z):
    """
    Fit an ellipse to the given Y, Z points using scipy.optimize.

    Parameters:
        Y (list or np.array): Y-coordinates of the points.
        Z (list or np.array): Z-coordinates of the points.

    Returns:
        tuple: (a, b, center, angle) parameters of the fitted ellipse.
    """
    Y, Z = np.array(Y), np.array(Z)

    # Initial guess: centered ellipse with no rotation
    y0, z0 = np.mean(Y), np.mean(Z)
    a0, b0 = (np.max(Y) - np.min(Y)) / 2, (np.max(Z) - np.min(Z)) / 2
    angle0 = 0  # Initial rotation angle

    def ellipse_loss(params):
        y0, z0, a, b, angle = params
        cos_angle = np.cos(angle)
        sin_angle = np.sin(angle)

        if(a<b):
            a,b = b,a

        # Transform points into the rotated frame of the ellipse
        Y_rot = cos_angle * (Y - y0) + sin_angle * (Z - z0)
        Z_rot = -sin_angle * (Y - y0) + cos_angle * (Z - z0)

        # Distance from the ellipse equation
        distances = ((Y_rot / a) ** 2 + (Z_rot / b) ** 2) - 1
        return np.sum(distances ** 2)

    # Optimize parameters
    initial_params = [y0, z0, a0, b0, angle0]
    bounds = [(None, None), (None, None), (1e-3, None), (1e-3, None), (-np.pi, np.pi)]
    result = minimize(ellipse_loss, initial_params, bounds=bounds)

    # Extract results
    y0, z0, a, b, angle = result.x
    return a, b, (y0, z0), angle


def covariance_from_ellipse(a, b, theta):

    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta), np.cos(theta)]])
    
    if theta < -np.pi/4:
        a, b = b, a
    A = np.array([[a**2, 0],
                  [0, b**2]])
    return R @ A @ R.T


thread = 10
import numpy as np
import os
from ROOT import TFile
from root2mpl import *
import matplotlib.pyplot as plt
import subprocess
from scipy.ndimage import gaussian_filter
save_file = "data_save1_thin_wxs.npy"

if os.path.exists(save_file):
    saved_data = np.load(save_file, allow_pickle=True).tolist()
    saved_dict = {tuple(entry[0]): entry[1] for entry in saved_data}  # Convert to dict for fast lookup
else:
    saved_data = np.array([], dtype=[('parameters', object), ('chi2', float)])
    saved_dict = {}


min=1000
# Initialize the rest of the arrays and settings as before
data_x = [[[] for _ in range(ndim)] for j in range(ndim)]
data_y = [[[] for _ in range(ndim)] for j in range(ndim)]
data_chi2 = [[[] for _ in range(ndim)] for j in range(ndim)]

for i, irange in enumerate(range_value):
    for j, jrange in enumerate(range_value):
        print(labels[i], labels[j], i, j)
        if i < j:
            continue

        elif i == j:
            for ix, x in enumerate(np.arange(irange[0], irange[1] + step_value[i], step_value[i])):
                x = round(x, 4)
                # print(f"{x}{unit_value[i]}")

                # Create parameter set
                parameter = minimum_value.copy()
                parameter[i] = x

                param_tuple = tuple(parameter)
                if param_tuple in saved_dict:
                    chi2 = saved_dict[param_tuple]
                    # print(f"Found saved value for {param_tuple}: {chi2}")
                else:
                    # Compute chi2 as before
                    command = f"./RBS_All_forMCMC {parameter[0]} {parameter[1]} -1111 -1111 {parameter[2]} {parameter[3]} -1111 -1111 {thread}"
                    subprocess.run(command, shell=True, stderr=subprocess.PIPE)

                    try:
                        with open(f"tmp/{thread}.txt", "r") as f:
                            lines = f.readlines()
                            chi2 = np.sum([float(line) for line in lines])
                    except (IndexError, ValueError):
                        chi2 = 0

                    # print(f"Chi2 for {param_tuple}: {chi2}")

                    saved_data.append([param_tuple, chi2])
                    saved_dict[param_tuple] = chi2

                data_x[i][j].append(x)
                data_y[i][j].append(x)
                data_chi2[i][j].append(chi2)
            
            np.save(save_file, saved_data)

            ax[i, j].scatter(data_x[i][j], data_chi2[i][j], s=10)
            ax[i, j].set_title(f"{labels[i]}")
            if j != ndim - 1:
                ax[i, j].set_xticklabels([])

            ax[i, j].tick_params(axis='y', labelsize=8)
            ax[i, j].yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
            ax[i, j].set_xlim(irange[0], irange[1])
            ax[i, j].set_xticks([irange[0], (irange[0]+irange[1])/2, irange[1]])
            ax[i, j].tick_params(axis='x', labelsize=10)
            if i == 0:
                ax[i, j].set_ylabel(r"$\chi^2_{\nu}$")

            # print(f"Processing {labels[i]} vs {labels[j]}")

        else:
            for ix, x in enumerate(np.arange(irange[0], irange[1] + step_value[i], step_value[i])):
                for iy, y in enumerate(np.arange(jrange[0], jrange[1] + step_value[j], step_value[j])):
                    x = round(x, 4)
                    y = round(y, 4)
                    # print(f"{x}{unit_value[i]} and {y}{unit_value[j]}")

                    parameter = minimum_value.copy()
                    parameter[i] = x
                    parameter[j] = y

                    param_tuple = tuple(parameter)
                    if param_tuple in saved_dict:
                        chi2 = saved_dict[param_tuple]
                        # print(f"Found saved value for {param_tuple}: {chi2}")
                    else:
                        command = f"./RBS_All_forMCMC {parameter[0]} {parameter[1]} -1111 -1111 {parameter[2]} {parameter[3]} -1111 -1111 {thread}"
                        print(command)
                        subprocess.run(command, shell=True, stderr=subprocess.PIPE)

                        try:
                            with open(f"tmp/{thread}.txt", "r") as f:
                                lines = f.readlines()
                                chi2 = np.sum([float(line) for line in lines])
                        except (IndexError, ValueError):
                            chi2 = 0

                        saved_data.append([param_tuple, chi2])
                        saved_dict[param_tuple] = chi2

                    data_x[i][j].append(x)
                    data_y[i][j].append(y)
                    data_chi2[i][j].append(chi2)

            np.save(save_file, saved_data)

            hist, xedges, yedges = np.histogram2d(
                data_y[i][j],
                data_x[i][j],
                bins=(
                    np.arange(jrange[0], jrange[1], step_value[j]),
                    np.arange(irange[0], irange[1], step_value[i])
                ),
                weights=data_chi2[i][j]
            )

            # if np.min(data_chi2[i][j]) < min:
            #     min = np.min(data_chi2[i][j])
            #     print(min)
            #     print(data_x[i][j][np.argmin(data_chi2[i][j])], data_y[i][j][np.argmin(data_chi2[i][j])])
                

            smoothed_hist = gaussian_filter(hist, sigma=1)
            img = ax[i, j].imshow(smoothed_hist.T, extent=[jrange[0], jrange[1], irange[0], irange[1]],
                            origin='center', aspect='auto', vmax = 20, vmin = 4, cmap='Blues_r')
            
            ax[i, j].set_ylim(irange[0], irange[1])
            ax[i, j].set_xlim(jrange[0], jrange[1])
            
            ax[i, j].plot([jrange[0], jrange[1]], [minimum_value[i], minimum_value[i]], color='red')
            ax[i, j].plot([minimum_value[j], minimum_value[j]], [irange[0], irange[1]], color='red')
            ax[i, j].scatter(minimum_value[j], minimum_value[i], color='red', s=10)

            if j != 0:
                ax[i, j].set_yticklabels([])
            if i != ndim - 1:
                ax[i, j].set_xticklabels([])
            if j == 0:
                ax[i, j].set_ylabel(f"{labels[i]}" + unit_value[i], fontsize=10)
                ax[i, j].tick_params(axis='y', labelsize=10)
            if i == ndim - 1:
                ax[i, j].set_xlabel(f"{labels[j]}" + unit_value[j], fontsize=10)
                ax[i, j].tick_params(axis='x', labelsize=10)
            ax[i, j].set_xticks([jrange[0], (jrange[0]+jrange[1])/2, jrange[1]])
            ax[i, j].set_yticks([irange[0], irange[0]+abs((irange[0]-irange[1]))/4, (irange[0]+irange[1])/2, irange[1]-abs(irange[0]-irange[1])/4, irange[1]])



            if (i == 3 and j == 2):
                continue

            contour= ax[i, j].contour(smoothed_hist.T, levels=[np.min(data_chi2[i][j]) + find_delta(1, 2)], extent=[jrange[0], jrange[1], irange[0], irange[1]], colors='none')
            
            #################################################""
            paths = contour.collections[0].get_paths()  
            vertices = paths[0].vertices                

            vertices = vertices[(vertices[:, 0] >= jrange[0]) & (vertices[:, 0] <= jrange[1])]
            vertices = vertices[(vertices[:, 1] >= irange[0]) & (vertices[:, 1] <= irange[1])]
        
            x_contour, y_contour = vertices[:, 0], vertices[:, 1]

            par = fit_ellipse(x_contour, y_contour)
            a_fit, b_fit, cx_fit, cy_fit, theta_fit = par[0], par[1], par[2][0], par[2][1], par[3]

            t_fitted = np.linspace(0, 2 * np.pi, 500)
            x_fitted, y_fitted = ellipse(t_fitted, cx_fit, cy_fit, a_fit, b_fit, theta_fit)

            ax[i, j].plot(x_fitted, y_fitted, color='black', label='Fitted Ellipse', linewidth=1.5, ls='--') 

            cov_matrix.append(covariance_from_ellipse(a_fit, b_fit, theta_fit))
            print(theta_fit)
            print("Covariance Matrix:")
            print(cov_matrix[-1])
            sigma_x, sigma_y = np.sqrt(np.diag(cov_matrix[-1]))
            print(f"1-sigma Errors: σ_x = {sigma_x}, σ_y = {sigma_y}")

            #################################################""

            
            print(np.min(data_chi2[i][j]), np.max(data_chi2[i][j]))
            print(data_x[i][j][np.argmin(data_chi2[i][j])], data_y[i][j][np.argmin(data_chi2[i][j])])
            print(labels[i], labels[j])
    
        plt.savefig("correlation_plot1_thin.png", dpi=300)

### FINAL COVARIANCE MATRIX ###
        
Final_Cov = np.zeros((4, 4))

Final_Cov[0, 0] = np.sum([cov_matrix[0][0, 0], cov_matrix[1][0, 0], cov_matrix[3][0, 0]])
Final_Cov[1, 1] = np.sum([cov_matrix[0][1, 1], cov_matrix[2][0, 0], cov_matrix[4][0, 0]])
Final_Cov[2, 2] = np.sum([cov_matrix[1][1, 1], cov_matrix[2][1, 1]])
Final_Cov[3, 3] = np.sum([cov_matrix[3][1, 1], cov_matrix[4][1, 1]])

Final_Cov[0, 1] = cov_matrix[0][0, 1]
Final_Cov[0, 2] = cov_matrix[1][0, 1]
Final_Cov[0, 3] = cov_matrix[3][0, 1]

Final_Cov[1, 2] = cov_matrix[2][0, 1]
Final_Cov[1, 3] = cov_matrix[4][0, 1]

# Final_Cov[2, 3] = cov_matrix[5][0, 1]

#make it symetric 

for i in range(4):
    for j in range(i+1, 4):
        Final_Cov[j, i] = Final_Cov[i, j]

print("Final Covariance Matrix:")
print(Final_Cov)

for i in range(4):
    print(f"\sigma {i}: {np.sqrt(Final_Cov[i, i])}")
        
