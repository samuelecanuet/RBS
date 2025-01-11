# from normal_corner import normal_corner
import numpy as np

mean = np.array([22.0887, 1.9681])
varlabels = ['x1', 'x2']

# ## covariance_matrix
# #      |      0    |      1    |      2    |      3    |      4    |
# # ----------------------------------------------------------------------
# #    0 |      24.73      0.1225       80.49     0.05598       106.5 
# #    1 |     0.1225   0.0006597       0.408   0.0002835      0.5399 
# #    2 |      80.49       0.408       272.8       0.186       354.7 
# #    3 |    0.05598   0.0002835       0.186    0.000253      0.2466 
# #    4 |      106.5      0.5399       354.7      0.2466       470.1 
# #    5 |      116.4      0.5902       387.7      0.2696       512.8 
# #    6 |      52.23      0.2648       173.9      0.1209         230 
# #    7 |      261.2       1.324       869.7      0.6047        1150 


# #      |      5    |      6    |      7    |
# # ----------------------------------------------------------------------
# #    0 |      116.4       52.23       261.2 
# #    1 |     0.5902      0.2648       1.324 
# #    2 |      387.7       173.9       869.7 
# #    3 |     0.2696      0.1209      0.6047 
# #    4 |      512.8         230        1150 
# #    5 |      562.5       251.5        1257 
# #    6 |      251.5       113.2         564 
# #    7 |       1257         564        2830


covm = np.array([[24.73, 0.1225],
                  [0.1225, 0.0006597]])


# figure_1 = normal_corner.normal_corner(covariance_matrix,mean_vector,variable_labels)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# If running from command line, above should be ran before
# importing normal_corner
import numpy as np
from normal_corner import normal_corner

### EXAMPLE 1: plotting one covariance matrix ###

# Covariance matrix, as a numpy array
covm = np.array([[24.73,0.1225,80.49], [0.1225, 0.01 ,0.408], [80.49,0.408, 272.8]])
print(covm)

# Mean matrix, as a numpy array
mean = np.array([22, 1.96, 70])
print(mean)

# Variable labels for plotting, as a list of strings
# in LaTeX format, between $$ symbols
varlabels = ['$var_1$','$var\\frac{2}{3}$','$var3^3$']

# Make a corner plot
fig1 = normal_corner.normal_corner(covm,mean,varlabels)
# plt.savefig('example_1.png')
plt.close()
