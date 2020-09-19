#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

plt.style.use("default")
plt.rc("text", usetex=True)
plt.rc('font', family='serif', size=18)
plt.rc('xtick', labelsize=17)
plt.rc('ytick', labelsize=17)
plt.rc('legend', fontsize=16)

n_in = 1.46
n_pol = 1.6
n_layers_list = [27, 32, 40]
d = 10
Ds = [20, 40, 60, 80, 100, 200]
n_inter_list = [1, 1.2, 1.4, 1.6]
n_layers = 25

wl_list = np.linspace(400, 1000, 400)
#wl_list = np.linspace(570, 610, 200)

T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()

fig, axs = plt.subplots(2, 3, sharex=True, sharey=True, gridspec_kw={'hspace':0.0, 'wspace':0.00}, figsize=(32, 10))
#for i, n_inter in enumerate(n_inter_list):
for i, D in enumerate(Ds):
    n_pol_list = np.ones(n_layers - 1) * n_pol
#    for D in Ds:
    for n_inter in n_inter_list:
        R = np.zeros(T.n_wls)
        for j in range(T.n_wls):
            T.trf_routine_sym(n_in, n_pol_list, n_in, d, n_inter, D, j)
            T.get_R()
            R[j] = T.R*100
        axs[i//3, i%3].plot(wl_list, R, label=r"$n_D$={:.2f}".format(n_inter))
    if i == 0:
        axs[i//3, i%3].legend()
    axs[i//3, i%3].set_xticks(np.arange(wl_list[0], wl_list[-1], 125))
    axs[i//3, i%3].margins(0)
#        axs[i].plot(wl_list, R, label=r"$n_f$={:.2f}".format(n_out))
#        axs[i].legend()
#        axs[i].set_title(r"$n_i$ = {}".format(n_in))
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
plt.xlabel(r"Wavelength (nm)")
plt.ylabel("Reflectance (\%)", labelpad=10)
plt.yticks(visible=True)
#plt.text(0.01, 0.525, "$n_i$ = {}".format(n_in_list[0]), fontsize=18)
#plt.text(0.01, 0.45, "$n_i$ = {}".format(n_in_list[1]), fontsize=18)
plt.title(r"$Ds \in [20, 200]$ \quad d = {}nm \quad $n_i$ = {} \quad $n_p$ = {} \quad \# layers = {}".format(d, n_in, n_pol_list[0], n_layers), pad=30)
plt.savefig("sym_D.pdf".format(n_layers), dpi = 600, bbox_inches="tight")
plt.show()
