#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
plt.style.use("default")
plt.rc("text", usetex=True)
plt.rc('font', family='serif', size=18)
plt.rc('xtick', labelsize=17)
plt.rc('ytick', labelsize=17)
plt.rc('legend', fontsize=16)

n_in = 1.8
n_in_list = [1.46, 1.67]
n_pol = 1.6
n_out_list = [1, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8]
n_layers_list = [27, 30, 40]
d = 10
w = 0

wl_list = np.linspace(650, 900, 200)
#wl_list = np.linspace(570, 610, 200)

T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()

fig, axs = plt.subplots(2, 3, sharex=True, sharey=True, gridspec_kw={'hspace':0.0, 'wspace':0.00}, figsize=(32, 10))
#axins = zoomed_inset_axes(axs[1, 0], 3.8, loc="upper center")
#plt.yticks(visible=False)
#plt.xticks(ticks=np.array([589.3]), fontsize=13)
for i in range(6):
    n_layers = n_layers_list[i%3]
    d_list = np.random.normal(size = n_layers, scale = w, loc = 1.0) * d
    n_in = n_in_list[i//3]
    n_pol_list = np.ones(n_layers - 1) * n_pol
    for n_out in n_out_list:
        R = np.zeros(T.n_wls)
        for j in range(T.n_wls):
            T.trf_routine(n_in, n_pol_list, n_out, d_list, j)
            T.get_R()
            R[j] = T.R*100
        axs[i//3, i%3].plot(wl_list, R, label=r"$n_f$={:.2f}".format(n_out))
        if i == 3:
#            axs[i//3, i%3].axvline(x=589.6, linestyle="--", color="black", alpha=.1, lw= .5)
#            axins.plot(wl_list[(wl_list>560)&(wl_list<630)], R[(wl_list>560)&(wl_list<630)])
            pass
    if i == 2:
        axs[i//3, i%3].legend(loc="lower right")
    if i < 3:
        axs[i//3, i%3].set_title(r"{} layers".format(n_layers), pad = -20)
    axs[i//3, i%3].set_xticks(np.arange(wl_list[0], wl_list[-1], 125))
    axs[i//3, i%3].margins(0)
#        axs[i].plot(wl_list, R, label=r"$n_f$={:.2f}".format(n_out))
#        axs[i].legend()
#        axs[i].set_title(r"$n_i$ = {}".format(n_in))
#    axins.axvline(589.3, linestyle="--", linewidth=.5, color=".5")
#    axins.margins(0)
#    mark_inset(axs[1, 0], axins, loc1=4, loc2=2, fc="none", ec="0.5")
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
plt.xlabel(r"Wavelength (nm)")
plt.ylabel("Reflectance (\%)", labelpad=10)
#plt.yticks(visible=True)
plt.text(0.01, 0.525, "$n_i$ = {}".format(n_in_list[0]), fontsize=18)
plt.text(0.01, 0.45, "$n_i$ = {}".format(n_in_list[1]), fontsize=18)
#plt.title(r"d = 10nm", pad=30)
#plt.savefig("plot1.pdf".format(n_layers), dpi = 600, bbox_inches="tight")
plt.show()
