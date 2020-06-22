#!/usr/bin/env python3


import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt

n_in = 1.46
n_pol = 1.6
n_out = 1.0
d_list = [10, 20, 30, 40]
wls = np.linspace(400, 2000, 200)
T = trf.TransferMatrix(300, -0.2, wls)
#ds = np.linspace(2, 100, 200)
n_20layer = np.ones(2) * n_pol

T.get_constant_xi(.149)
#T.get_graphene_xi()
#T.init_trf_matrix()
#T.multiply_by_layer(n_in, n_pol, 0)
#for i in range(18):
#    T.multiply_by_chunk(n_pol, n_pol, d_list[1], 0)
#T.multiply_by_chunk(n_pol, n_out, d_list[1], 0)
#print(T.trf_matrix)
T.plot_r_vs_wl(n_in, n_20layer, n_out, d_list)
#T.plot_r_vs_wl_2layer(n_in, n_pol, n_out, d_list)
#T.plot_r_heatmap(n_in, n_20layer, n_out, ds)
#plt.show()
