#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import trf_matrix as trf
import scipy.optimize as opt
import os
os.chdir("/home/rmp/Projects/Transfer_Matrix/Joule_Heating/Data_2a")

n_in = 1.46
n_pol = 1.6
n_out_list = np.arange(1, 3, .2)
#n_out_list = [1.8, 2.8]
xi = 0.149
d = 10
wl_list = np.arange(600, 2600, 200)
#wl_list = [800, 1800]

T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_constant_xi(xi)

n_layers = [19, 39, 59, 79, 99]
#ds = 500 / np.array(list(n_layers), dtype=float)
ds = np.ones(5, dtype=int) * d

def linfunc(x, a, b):
    return a*x + b

for n_out in n_out_list:
    for i, wl in enumerate(wl_list):
        x = np.arange(n_layers[-1] + 1)
        for d, n_layer in zip(ds, n_layers):
            n_pol_list = np.ones(n_layer) * n_pol
            E = T.get_field(n_in, n_pol_list, n_out, d, i)
            y = np.absolute(E)**2
            popt, pcov = opt.curve_fit(linfunc, x[:len(y)//5], y[:len(y)//5])
            data = np.column_stack((np.arange(len(y), dtype=int), y))
            plt.plot(x[:len(y)], y, label="{} layers".format(n_layer + 1))
            residuals = y[:len(y)//5] - linfunc(x[:len(y) // 5], *popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y[:len(y)//5] -np.mean(y[:len(y)//5]))**2)
            r_squared = 1 - (ss_res / ss_tot)
            print(r_squared)
            with open("n_out:{:.1f}_wl:{}_nlayers:{}_d:{}.dat".format(n_out, wl, n_layer+1, d), "w+") as f:
                np.savetxt(f, data, delimiter="\t", header="nlayers\t|E|^2\tslope = {}\tr^2  = {}".format(popt[0], r_squared))
        plt.yscale("log")
        plt.title(r"$\lambda = ${} | $n_f$ $=$ {:.1f}".format(wl, n_out))
        plt.legend()
        plt.savefig("plot_n_out:{:.1f}_wl:{}_nlayers:{}.pdf".format(n_out, wl, n_layer+1), dpi=800, bbox_inches="tight")
        plt.close()
#        plt.show()
