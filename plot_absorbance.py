#!/usr/bin/env python3

#!/usr/bin/env python3
#
import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt

n_in = 1.46
n_pol = 1.6
n_layer = np.ones(19) * n_pol
d = 10
wl_list = [800, 1000, 1200, 1400, 1600, 1800]
n_out = 1

T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()

T.plot_r_vs_layers(n_in, n_layer, n_out, d)

plt.show()
