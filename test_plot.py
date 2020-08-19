#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt

n_in = 1.46
n_pol = 1.6
n_out = 1.0
xi = 0.149
n_layer = np.ones(99) * n_pol
d = 10
d_list = [5, 10, 15, 20]
wl_list = [800, 1000, 1300, 1500, 1800, 2000]

T = trf.TransferMatrix(300, -0.2, wl_list)

T.get_constant_xi(xi)
T.plot_field_sum_vs_layers(n_in, n_layer, n_out, d)
plt.show()
