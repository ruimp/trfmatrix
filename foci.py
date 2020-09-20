#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import os
os.chdir("/home/rmp/Projects/TransferMatrix/Foci")

n_in_list = [1.46, 1.77]
n_pol = 1.6
n_out1 = 1.3
n_out2 = 1.7
n_layers_list = [20, 25, 30, 35, 40, 45]
ds = np.linspace(5, 36, 50)
wl_list = np.linspace(300, 1400, 200)
T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()

w_list = [0, .1, .2, .3]

for n_in in n_in_list:
    for n_layers in n_layers_list:
        n_pol_list = np.ones(n_layers-1) * n_pol
        for w in w_list:
            foci = np.zeros((len(ds), 2))
            for i, d in enumerate(ds):
                d_list = d * np.random.normal(loc=1, scale=w, size=n_layers-1)
                R1 = np.zeros(T.n_wls)
                R2 = np.zeros(T.n_wls)
                for j in range(T.n_wls):
                    T.trf_routine(n_in, n_pol_list, n_out1, d_list, j)
                    T.get_R()
                    R1[j] = T.R*100
                    T.trf_routine(n_in, n_pol_list, n_out2, d_list, j)
                    T.get_R()
                    R2[j] = T.R*100
                R = R2 - R1
                foci[i, 0] = d
                foci[i, 1] = len(np.where(np.diff(np.sign(R)))[0])

            filename = "N{}:w{}:ni{}:np{}:wl[300,1400].dat".format(n_layers, w, n_in, n_pol)
            np.savetxt(filename, foci, delimiter="\t")
