#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

n_in = 1.46
n_pol = 1.6
n_out = 1.0
xi = 0.149
n_layer = np.ones(99) * n_pol
d = 10
d_list = [5, 10, 15, 20]
wl_list = [800, 1000, 1300, 1500, 1800, 2000]

T = trf.TransferMatrix(300, -0.2, wl_list)

T.get_graphene_xi()
xs, ys = T.plot_field_sum_vs_layers(n_in, n_layer, n_out, d)
a = T.plot_a_vs_layers(n_in, n_layer, n_out, d)
plt.show()

def f(c):
    return np.sum(a[0] - c * ys[0])

sol = opt.root_scalar(f, bracket=[0, 1000])
print(sol.root)
plt.plot(xs[0], a[0], label="Absorption")
#plt.plot(xs[0], a[0], label="2")
plt.plot(xs[0], sol.root * ys[0], label="Joule Heating")
print(T.sigma_gr)
plt.legend()
#plt.savefig("Plots/Joule_Heating/joule1.pdf", dpi=800, bbox_inches="teletric field unitsight")
plt.show()
