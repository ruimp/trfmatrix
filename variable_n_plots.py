#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt

n_in = 1.46
n_pol_list = np.arange(1.2, 2, .2)
n_out = 1.0
xi = 0.149
d = 10
d_list = [5, 10, 15, 20]
wls = np.linspace(500, 1800, 1000)
n_layers = 20
T = trf.TransferMatrix(300, -0.2, wls)
T.get_constant_xi(xi)


fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, gridspec_kw={'hspace':0, 'wspace':0.03})
for i, d in enumerate(d_list):
    for n_pol in n_pol_list:
        R = np.zeros(T.n_wls)
        for j in range(T.n_wls):
            T.init_trf_matrix()
            T.multiply_by_layer(n_in, n_pol, j)
            for p in range(n_layers - 2):
                T.multiply_by_chunk(n_pol, n_pol, d, j)
            T.multiply_by_chunk(n_pol, n_out, d, j)
            T.get_coeffs(n_in, n_out)
            R[j] =T.R*100
        axs[i//2, i%2].plot(T.wl_list, R, label="n = {:.2f} | D = {}nm".format(n_pol, d))
    axs[i//2, i%2].legend()
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Reflectance", labelpad=10)
plt.savefig("Plots/Plot_vs_Npol/R_vs_WL_vs_npol_20layers_Xic.pdf", dpi=1200, bbox_inches="tight")
plt.show()
