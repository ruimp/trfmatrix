#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

n_in = 1.46
n_pol = 1.6
n_out = 1.0
n_layer = np.ones(99) * n_pol
d_list = np.arange(5, 55, 5)
wl_list = np.arange(600, 2600, 200)

T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()

coeffs = np.zeros((len(d_list) * T.n_wls, 5))

def f(c, i):
    return np.sum(a[i] - c * ys[i])

c = 3.0e8
e0 = 8.854*10**(-12)

for i, d in enumerate(d_list):
    xs, ys = T.plot_field_sum_vs_layers(n_in, n_layer, n_out, d)
    a = T.get_a_list(n_in, n_layer, n_out, d)
    for j in range(T.n_wls):
        sol = opt.root_scalar(lambda c: f(c, j), bracket=[0, 10000])
#        plt.plot(xs[j], a[j], label="Absorption")
#        plt.plot(xs[j], sol.root * ys[j], label="Joule Heating")
#
        b = T.sigma_gr.real[j] / n_pol / e0 / c
        coeffs[10*i + j] = (d, wl_list[j], sol.root, b, sol.root / b)
#        plt.legend()
#        plt.savefig("Plots/Joule_Heating/d{}_wl{}.pdf".format(d, wl_list[i]), dpi=800, bbox_inches="tight")
#        plt.show()

np.savetxt("Joule_Heating/coeffs.dat", coeffs, delimiter = "\t", header="d\twl\tcalculatedcalculated  coeff\ttheoretical value\tratio".format(d_list, wl_list))
