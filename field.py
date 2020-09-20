#!/usr/bin/env python3
#
import trf_matrix as trf
import numpy as np
import os
import matplotlib.pyplot as plt
os.chdir("/home/rmp/Projects/TransferMatrix/Field")

n_in = 1.46
n_pol = 1.6
n_out_list = [1.2, 1.4, 1.6, 1.8]
d = 10
n_layers = 30
wl_list = [650, 700, 750]
s = [-1, 0, 1]
n_pol_list = np.ones(n_layers - 1) * n_pol
d_list = np.ones(n_layers) * d
T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()
x = np.arange(n_layers)

n_out = n_out_list[0]
for n_out in n_out_list:
    for j in range(T.n_wls):
        data = np.zeros((n_layers, 2))
        E = T.get_field(n_in, n_pol_list, n_out, d_list, j)
        y = np.absolute(E)**2
        data[:, 0] = x
        data[:, 1] = y
        filename = "N{}:ni{}:np{}:nf{}:wl{}:s:{}.dat".format(n_layers, n_in, n_pol, n_out, wl_list[j], s)
        np.savetxt(filename, data, delimiter="\t")
