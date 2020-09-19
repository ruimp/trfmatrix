#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("default")

n_in = 1.8
n_in_list = [1.3, 1.5]
n_pol = 1.6
n_out_list = [1, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8]

n_layers_list = [25, 30, 39]
d = 10

#d_list = np.random.normal(size = n_layers - 1, scale = .8, loc = 1.0) * d

wl_list = np.linspace(400, 1000, 200)
#wl_list = np.linspace(570, 610, 200)

T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()

fig, axs = plt.subplots(2, 3, sharex=True, sharey=True, gridspec_kw={'hspace':0.0, 'wspace':0.00})
#for i, n_layers in enumerate(n_layers_list):
#for i, n_in in enumerate(n_in_list):

wl_0 = 589.3
dl = 10
gaussian = np.exp(-(wl_list[(wl_list >= wl_0 - 6 * dl)&(wl_list <= wl_0 + 6* dl)] - wl_0)**2 / 2 / dl) / np.sqrt(2 * np.pi * dl**2)

for i in range(6):
    n_layers = n_layers_list[i%3]
    n_in = n_in_list[i//3]
    n_pol_list = np.ones(n_layers - 1) * n_pol
    for n_out in n_out_list:
        R = np.zeros(T.n_wls)
        for j in range(T.n_wls):
#            T.trf_routine_disordered(n_in, n_pol_list, n_out, d_list, j)
            T.trf_routine(n_in, n_pol_list, n_out, d, j)
            T.get_R()
            R[j] = T.R*100
        R = np.convolve(R, gaussian, mode="same")
        axs[i//3, i%3].plot(wl_list, R, label=r"$n_f$={:.2f}".format(n_out))
        if i == 2:
            axs[i//3, i%3].legend(loc="lower right")
        if i < 3:
            axs[i//3, i%3].set_title(r"{} layers".format(n_layers), pad = -20)
axs[0, 0].text(5, 450, "text")
#        axs[i].plot(wl_list, R, label=r"$n_f$={:.2f}".format(n_out))
#        axs[i].legend()
#        axs[i].set_title(r"$n_i$ = {}".format(n_in))
plt.xticks(ticks=np.arange(wl_list[0], wl_list[-1], 125))
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
plt.xlabel(r"Wavelenght (nm)")
plt.ylabel("Reflectance (%)", labelpad=10)
plt.text(.4, .4, "Text.")
#plt.title(r"d = 10nm", pad=30)
#plt.savefig("Plots/R_vs_WL/R_vs_wl_for_nlayers_ni1.46.pdf".format(n_layers), dpi=800, bbox_inches="tight")
plt.show()
