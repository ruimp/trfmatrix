#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt

n_in = 1.8
n_in_list = [1.2, 1.46, 1.8]
n_pol = 1.6
n_out_list = [1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]
n_layers_list = [15, 20, 25, 30, 35, 40]
n_layers = 40

d = 10
wl_list = np.linspace(500, 1600, 500)

T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()

fig, axs = plt.subplots(3, 2, sharex=True, sharey=True, gridspec_kw={'hspace':0.16, 'wspace':0.03})
#for i, n_layers in enumerate(n_layers_list):
for i, n_layers in enumerate(n_layers_list):
    n_pol_list = np.ones(n_layers - 1) * n_pol
    wlmax = np.zeros(len(n_out_list))
    for k, n_out in enumerate(n_out_list):
        R = np.zeros(T.n_wls)
        for j in range(T.n_wls):
            T.trf_routine(n_in, n_pol_list, n_out, d, j)
            T.get_R()
            R[j] = T.R*100
        for n in range(1, T.n_wls-1):
            if R[n-1] <= R[n] and R[n+1] <= R[n]:
                wlmax[k] = wl_list[n]
        print(wlmax)
        axs[i//2, i%2].plot(wl_list, R, label=r"$n_f$={:.2f}".format(n_out))
#        axs[i].legend()
#        axs[i].set_title(r"$n_i$ = {}".format(n_in))
#plt.xticks(ticks=np.arange(2, n_layers + 1, 4, dtype=int))
fig.add_subplot(111, frameon=False)
#plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
#plt.xlabel(r"Wavelenght (nm)")
#plt.ylabel("Reflectance (%)", labelpad=10)
#plt.title(r"{} layers | d = 10nm".format(n_layers, n_in), pad=30)
#plt.savefig("Plots/R_vs_WL/R_vs_wl_for_nlayers_ni1.46.pdf".format(n_layers), dpi=800, bbox_inches="tight")
plt.show()
