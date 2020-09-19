#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

plt.style.use("default")
#plt.rc("text", usetex=True)
plt.rc('font', family='serif', size=18)
plt.rc('xtick', labelsize=17)
plt.rc('ytick', labelsize=17)
plt.rc('legend', fontsize=16)

n_in = 1.46
n_pol = 1.6
n_out1 = 1.3
n_out2 = 1.7
#n_layers_list = np.arange(15, 41, 5)
n_layers = 30
ds = np.arange(6, 31)
wl_list = np.linspace(300, 1400, 200)
w1 = .05
w2 = .20

T = trf.TransferMatrix(300, -0.2, wl_list)
T.get_graphene_xi()
n_pol_list = np.ones(n_layers-1) * n_pol
focci = np.zeros(len(ds))

def f(x, a, b):
    return a*x + b

for w in [0, .05, .2, .5]:
    for i, d in enumerate(ds):
        d_list = d * np.random.normal(loc=1, scale=w, size=n_layers-1)
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
        focci[i] = len(np.where(np.diff(np.sign(R)))[0])
    popt, pcov = opt.curve_fit(f, ds, focci)
    plt.plot(ds, focci, "^")
    plt.plot(ds, f(ds, popt[0], popt[1]), label="w = {:.2f}".format(w))
plt.legend()
plt.yticks(ticks=np.arange(0, np.max(focci)+1, 2, dtype=int))
plt.xlabel(r"Polymer width (nm)")
plt.ylabel("Number of focci", labelpad=16)
#plt.savefig("Plots/Focci/focci_numb.pdf", dpi=600, bbox_inches="tight")
plt.show()
