#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import os
import sys
os.chdir("/home/rmp/Projects/TransferMatrix/TEST")

n_in = 1.444
wl_list = np.linspace(1, 10000, 10000)
n_out_list = np.linspace(1.3, 1.6, 20)
n_pol = 1.6
d = 10

T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_constant_xi(.15)

n_layers = int(sys.argv[-1])

n_pol_list = np.ones(n_layers - 1) * n_pol
d_list = np.ones(n_layers) * d
for n_out in n_out_list:
    data = np.zeros((T.n_wls, 3))
    for j in range(T.n_wls):
        T.trf_routine(n_in, n_pol_list, n_out, d_list, j)
        r = T.get_r()
        data[j, 0] = wl_list[j]
        data[j, 1] = np.absolute(r)
        data[j, 2] = np.angle(r)
    filename = "N{}:ni{:.3f}:np1.6:nf{:.3f}:d{}.dat".format(n_layers, n_in, n_out, d)
    np.savetxt(filename, data, delimiter="\t")
    print("done")
    #print(n_layers)
