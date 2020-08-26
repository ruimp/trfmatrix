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
xi = 0.149
d = 10
wl_list = np.arange(600, 2600, 200)

T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()

n_layers = [19, 39, 59, 79, 99]

def linfunc(x, a, b):
    return a*x + b

for n_out in n_out_list:
    for i, wl in enumerate(wl_list):
        x = np.arange(n_layers[-1] + 1)
        for n_layer in n_layers:
            n_pol_list = np.ones(n_layer) * n_pol
            E = T.get_field(n_in, n_pol_list, n_out, d, i)
            y = np.absolute(E)**2
            popt, pcov = opt.curve_fit(linfunc, x[:len(y) // 2], y[:len(y)//2])
            data = np.column_stack((np.arange(len(y), dtype=int), y, np.ones(len(y))*popt[0]))
            plt.plot(x[:len(y)], y, label="{} layers".format(n_layer + 1))
            with open("n_out:{:.1f}_wl:{}_nlayers:{}.dat".format(n_out, wl, n_layer+1), "w+") as f:
                np.savetxt(f, data, delimiter="\t", header="nlayers\t|E|^2\tslope")
        plt.yscale("log")
        plt.title(r"$\lambda = ${} | $n_f$ $=$ {:.1f}".format(wl, n_out))
        plt.legend()
        plt.savefig("plot_n_out:{:.1f}_wl:{}_nlayers:{}.pdf".format(n_out, wl, n_layer+1), dpi=800, bbox_inches="tight")
        plt.close()
