#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt

n_in = 1.46
n_pol = 1.6
n_out = 1.0
xi = 0.149
n_layer = np.ones(39) * n_pol
d = 10
d_list = [5, 10, 15, 20]
wl_list = [800, 1000, 1200]

T = trf.TransferMatrix(300, -0.2, wl_list)

T.get_constant_xi(xi)
T.plot_field_vs_d(n_in, n_layer, n_out, d_list, wl_ind = 0)
plt.show()


T.init_trf_matrix()
T.trf_routine(n_in, n_layer, n_out, d,0)
n_layers = n_layer.size
v = np.zeros((n_layers + 2, 2), dtype=complex)
v[0] = np.array([1, T.get_r()])
T.init_trf_matrix()
T.multiply_by_layer(n_in, n_layer[0], 0)
v[1] = T.trf_matrix.dot(v[0])
T.init_trf_matrix()
T.multiply_by_chunk(n_layer[0], n_layer[1], d, 0)
v[2] = T.trf_matrix.dot(v[1])
T.init_trf_matrix()
T.multiply_by_chunk(n_layer[1], n_out, d, 0)
v[3] = T.trf_matrix.dot(v[2])
E = v.sum(axis=-1)
plt.plot(np.absolute(E))
plt.show()
