#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import os
import matplotlib.pyplot as plt
os.chdir("/home/rmp/Projects/TransferMatrix/Reflectance")

n_in_list = [1.46, 1.67]
n_pol = 1.6
n_out_list = np.arange(1, 1.9, .1)
d_list = np.arange(5, 46, 5)
n_layers_list = np.arange(15, 41, 5)
wl_list = np.linspace(400, 1200, 400)

T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()

for n_in in n_in_list:
    for d in d_list:
        for n_layers in n_layers_list:
            d_list = np.ones(n_layers) * d
            n_pol_list = np.ones(n_layers - 1) * n_pol
            for n_out in n_out_list:
                data = np.zeros((T.n_wls, 2))
                for j, wl in enumerate(wl_list):
                    T.trf_routine(n_in, n_pol_list, n_out, d_list, j)
                    T.get_R()
                    data[j, 0] = wl
                    data[j, 1] = T.R
                filename = "N{}:ni{:.2f}:np{:.2f}:nf{:.2f}:d{}:wl[400,1200].dat".format(n_layers, n_in, n_pol, n_out, d)
                np.savetxt(filename, data, delimiter="\t")
