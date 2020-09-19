#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt

#n_in = 1.8
n_in_list = np.linspace(1, 2, 20)
n_pol = 1.6
n_layers_list = [30, 35, 40]

d = 10
wl_list1 = np.linspace(550, 750, 800)
wl_list2 = np.linspace(650, 900, 800)
wl_list3 = np.linspace(800, 1200, 800)

wl_lists = np.stack((wl_list1, wl_list2, wl_list3))

n_out1 = 1
n_out2 = 3



#fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, gridspec_kw={'hspace':0.14, 'wspace':0.03})
symbs = [".", "+", "^"]
for i, n_layers in enumerate(n_layers_list):
    wl_list = wl_lists[i]
    T = trf.TransferMatrix(300, -0.2, wl_list)
    T.get_graphene_xi()
    n_pol_list = np.ones(n_layers - 1) * n_pol
    focci_wl = np.zeros(len(n_in_list))
    focci_R = np.zeros(len(n_in_list))
    wl_0 = (wl_list[-1] + wl_list[0]) / 2
    dl = 10
    gaussian = np.exp(-(wl_list[(wl_list >= wl_0 - 6*dl)&(wl_list <= wl_0 + 6*dl)] - wl_0)**2 / 2 / dl) / np.sqrt(2 * np.pi * dl**2)
    for k, n_in in enumerate(n_in_list):
        R1 = np.zeros(T.n_wls)
        R2 = np.zeros(T.n_wls)
        for j in range(T.n_wls):
            T.trf_routine(n_in, n_pol_list, n_out1, d, j)
            T.get_R()
            R1[j] = T.R*100
            T.trf_routine(n_in, n_pol_list, n_out2, d, j)
            T.get_R()
            R2[j] = T.R*100
        R1 = np.convolve(R1, gaussian, mode="same")[1:-1]
        R2 = np.convolve(R2, gaussian, mode="same")[1:-1]
        wl_min = wl_list[np.argmin(np.absolute(R2-R1))]
        R_min = R1[np.argmin(np.absolute(R2-R1))]
        focci_wl[k] = wl_min
        focci_R[k] = R_min
    plt.plot(focci_wl, focci_R, symbs[i], label="{} layers".format(n_layers))
    plt.legend()
#    axs[i].plot(focci_wl, focci_R, ".")
#   axs[i].legend([r"$\Delta n_i$ = {:.2f}".format(n_in_list[2] - n_in_list[1])])
#    axs[i].set_title(r"{} layers".format(n_layers))
#plt.xticks(ticks=np.arange(2, n_layers + 1, 4, dtype=int))
#fig.add_subplot(111, frameon=False)
#plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
plt.xlabel(r"Wavelenght (nm)")
plt.ylabel("Reflectance (%)", labelpad=10)
plt.title(r"$n_i \in [1, 2]$ | $\Delta n_i$ = {:.2f} | d = 10nm".format(n_in_list[1] - n_in_list[0], n_in))
#plt.savefig("Plots/R_vs_WL/R_vs_wl_for_nlayers_ni1.46.pdf".format(n_layers), dpi=800, bbox_inches="tight")
plt.show()
