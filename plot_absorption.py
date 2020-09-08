#!/usr/bin/env python3
#
import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt

n_in = 1.46
n_pol = 1.6
n_layer = np.ones(19) * n_pol
#d_list = np.arange(5, 55, 5)
d = 10
wl_list = [800, 1000, 1200, 1400, 1600, 1800]
n_out_list = np.linspace(1, 3, 1000)

T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()

for i in range(T.n_wls):
    y = np.zeros(len(n_out_list))
    for j, n_out in enumerate(n_out_list):
        T.trf_routine(n_in, n_layer, n_out, d, i)
        T.get_coeffs(n_in, n_out)
        y[j] = T.A
    plt.plot(n_out_list, y, label=r"$\lambda$ = {}".format(wl_list[i]))
plt.legend()
plt.savefig("Plots/Absorption/absorption_vs_n_out.pdf", dpi=800, bbox_inches="tight")
plt.show()
