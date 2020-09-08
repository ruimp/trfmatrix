#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

n_in = 1.46
n_pol = 1.6
n_out = 1.0
n_layer = np.ones(99) * n_pol
d_list = np.arange(5, 55, 5)
d = 10
wl_list = np.linspace(400, 3000, 1000)

T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()


def f(x, a, b):
    return a * x + b

slopes = np.zeros(len(wl_list))
for i, wl in enumerate(wl_list):
    E = T.get_field(n_in, n_layer, n_out, d, i)
    y = np.absolute(E)**2
    popt, pcov = opt.curve_fit(f, np.arange(len(y)//4), y[:len(y)//4])
    slopes[i] = 1/popt[0]
plt.xlabel("Wavelenght (nm)")
plt.ylabel("Inverse of slope - 1 / m")
plt.plot(wl_list, slopes)
plt.savefig("Plots/Joule_Heating/2c.pdf", dpi=800, bbox_inches="tight")
plt.show()

wl_list2 = [1000]
T = trf.TransferMatrix(300, -0.2, wl_list2)
T.get_graphene_xi()
n_out_list = np.linspace(1, 3, 1000)
for i, n_out in enumerate(n_out_list):
    E = T.get_field(n_in, n_layer, n_out, d, 0)
    y = np.absolute(E)**2
    popt, pcov = opt.curve_fit(f, np.arange(len(y)//4), y[:len(y)//4], xtol = 1e-8)
    slopes[i] = - 1/popt[0]
plt.xlabel("Exit refractive index ")
plt.ylabel("Inverse of slope - 1 / m")
plt.plot(n_out_list, slopes)
plt.savefig("Plots/Joule_Heating/2d.pdf", dpi=800, bbox_inches="tight")
plt.show()
