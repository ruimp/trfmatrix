#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("default")
plt.rc("text", usetex=True)
plt.rc('font', family='serif', size=18)
plt.rc('xtick', labelsize=17)
plt.rc('ytick', labelsize=17)
plt.rc('legend', fontsize=16)

n_in = 1.46
n_pol = 1.6
n_out1 = 1.3
n_out2 = 1.7
n_layers_list = np.arange(15, 41, 5)
ds = np.arange(6, 31)
wl_list = np.linspace(300, 1400, 200)
w1 = .05
w2 = .20

fig, axs = plt.subplots(2, 3, sharex=True, sharey=True, gridspec_kw={'hspace':0.0, 'wspace':0.00}, figsize=(32, 10))
T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()
for k, n_layers in enumerate(n_layers_list):
    n_pol_list = np.ones(n_layers-1) * n_pol
    focci = np.zeros(len(ds))
    for i, d in enumerate(ds):
        d_list1 = d * np.random.normal(loc=1, scale=w, size=n_layers-1)
        d_list2 = d * np.random.normal(loc=1, scale=w, size=n_layers-1)
        R1 = np.zeros(T.n_wls)
        R2 = np.zeros(T.n_wls)
        for j in range(T.n_wls):
            T.trf_routine_disordered(n_in, n_pol_list, n_out1, d_list, j)
            T.get_R()
            R1[j] = T.R*100
            T.trf_routine_disordered(n_in, n_pol_list, n_out2, d_list, j)
            T.get_R()
            R2[j] = T.R*100
        R = R2 - R1
#        plt.plot(wl_list, R, label="d = {}".format(d))
#        plt.legend()
#        plt.show()
        focci[i] = len(np.where(np.diff(np.sign(R)))[0])
    axs[k//3, k%3].plot(ds, focci, "^", label="{} layers".format(n_layers), color="black")
    axs[k//3, k%3].legend()
    axs[k//3, k%3].set_yticks(ticks=np.arange(0, np.max(focci)+1, 2, dtype=int))
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
plt.xlabel(r"Polymer width (nm)")
plt.ylabel("Number of focci", labelpad=16)
plt.savefig("Plots/Focci/focci_number_w0.05.pdf", dpi=600, bbox_inches="tight")
plt.show()
