import numpy as np
import matplotlib.pyplot as plt


data = np.load("data_save (copy).npy", allow_pickle=True)

for par, chi2 in zip(data[0], data[1]):
    print(par, chi2)