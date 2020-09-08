#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt

n_in = 1.46
n_pol = 1.6
n_out_list = [1, 1.4, 1.8, 2.2, 2.6, 3]
n_layers = 40
n_pol_list = np.ones(n_layers - 1) * n_pol
d = 10
wl_list = np.arange(600, 2600, 200)

T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()

fig, axs = plt.subplots(3, 2, sharex=True, sharey=True, gridspec_kw={'hspace':0.8, 'wspace':0.03})
for i, n_out in enumerate(n_out_list):
    for j, wl in enumerate(wl_list):
        R = np.zeros(n_layers-1)
        for n in range(1, n_layers):
            T.trf_routine(n_in, n_pol_list[:n], n_out, d, j)
            T.get_R()
            R[n-1] = T.R*100
        axs[i//2, i%2].plot(np.arange(2, n_layers + 1, dtype=int), R, label=r"$\lambda$={}".format(wl))
        if i == 0:
            axs[i//2, i%2].legend()

plt.xticks(ticks=np.arange(2, n_layers + 1, 4, dtype=int))
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
plt.xlabel("Number of layers")
plt.ylabel("Reflectance (%)", labelpad=10)
plt.title(r"$\lambda \in$ [600, 2600, 200] | n $\in$ [1, 3, .4]")
plt.savefig("Week2/Plot_R_vs_40Layers.pdf", dpi=800, bbox_inches="tight")
plt.show()
