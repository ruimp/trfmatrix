#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import os
os.chdir("/home/rmp/Projects/TransferMatrix/Nlayers")

n_in = 1.46
n_pol = 1.6
n_out_list = np.arange(1, 2, .2)
n_layers_list = np.arange(15, 50, 5)
ds = np.arange(5, 45, 5)
wl_list = np.arange(400, 1400, 200)
T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()

for n_layers in n_layers_list:
    n_pol_list = np.ones(n_layers - 1) * n_pol
    for d in ds:
        d_list = np.ones(n_layers) * d
        for j in range(T.n_wls):
            for i, n_out in enumerate(n_out_list):
                data = np.zeros((n_layers, 2))
                R = np.zeros(n_layers)
                T.init_trf_matrix()
                T.multiply_by_layer(n_in, n_out, j)
                T.get_R()
                data[0, 0] = 1
                data[0, 1] = T.R
                for n in range(1, n_layers):
                    T.trf_routine(n_in, n_pol_list[:n], n_out, d_list[:n+1], j)
                    T.get_R()
                    R[n] = T.R*100
                    data[n, 0] = n+1
                    data[n, 1] = T.R
                filename = "N{}:ni{:.2f}:np{:.2f}:nf{:.2f}:d{:.2f}:wl{}.dat".format(n_layers, n_in, n_pol, n_out, d, wl_list[j])
                np.savetxt(filename, data, delimiter="\t")
