#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt

n_in = 1.46
n_pol = 1.6
n_out = 1.0
xi = 0.149
n_3layer = np.ones(2) * n_pol
n_10layer = np.ones(3) * n_pol
n_20layer = np.ones(19) * n_pol
d = 10
d_list = [5, 10, 15, 20]
wl_list = [600, 800, 1000, 1200]

T = trf.TransferMatrix(300, -0.2, wl_list)

T.get_constant_xi(xi)
T.plot_field_vs_d(n_in, n_3layer, n_out, d_list, wl_ind = 2)
plt.savefig("Plots/Plot_Field/field_vs_d_3layer_Xic_wl1000.pdf", dpi=800, bbox_inches="tight")
T.plot_field_vs_d(n_in, n_10layer, n_out, d_list, wl_ind = 2)
plt.savefig("Plots/Plot_Field/field_vs_d_10layer_Xic_wl1000.pdf", dpi=800, bbox_inches="tight")
T.plot_field_vs_d(n_in, n_20layer, n_out, d_list, wl_ind = 2)
plt.savefig("Plots/Plot_Field/field_vs_d_20layer_Xic_wl1000.pdf", dpi=800, bbox_inches="tight")

T.get_graphene_xi()
T.plot_field_vs_d(n_in, n_3layer, n_out, d_list, wl_ind = 2)
plt.savefig("Plots/Plot_Field/field_vs_d_3layer_Xig_wl1000.pdf", dpi=800, bbox_inches="tight")
T.plot_field_vs_d(n_in, n_10layer, n_out, d_list, wl_ind = 2)
plt.savefig("Plots/Plot_Field/field_vs_d_10layer_Xig_wl1000.pdf", dpi=800, bbox_inches="tight")
T.plot_field_vs_d(n_in, n_20layer, n_out, d_list, wl_ind = 2)
plt.savefig("Plots/Plot_Field/field_vs_d_20layer_Xig_wl1000.pdf", dpi=800, bbox_inches="tight")

T.get_constant_xi(xi)
T.plot_field_vs_wl(n_in, n_3layer, n_out, d)
plt.savefig("Plots/Plot_Field/field_vs_wl_3layer_Xic_d10.pdf", dpi=800, bbox_inches="tight")
T.plot_field_vs_wl(n_in, n_10layer, n_out, d)
plt.savefig("Plots/Plot_Field/field_vs
_wl_10layer_Xic_d10.pdf", dpi=800, bbox_inches="tight")
T.plot_field_vs_wl(n_in, n_20layer, n_out, d)
plt.savefig("Plots/Plot_Field/field_vs_wl_20layer_Xic_d10.pdf", dpi=800, bbox_inches="tight")

T.get_graphene_xi()
T.plot_field_vs_wl(n_in, n_3layer, n_out, d)
plt.savefig("Plots/Plot_Field/field_vs_wl_3layer_Xig_d10.pdf", dpi=800, bbox_inches="tight")
T.plot_field_vs_wl(n_in, n_10layer, n_out, d)
plt.savefig("Plots/Plot_Field/field_vs_wl_10layer_Xig_d10.pdf", dpi=800, bbox_inches="tight")
T.plot_field_vs_wl(n_in, n_20layer, n_out, d)
plt.savefig("Plots/Plot_Field/field_vs_wl_20layer_Xig_d10.pdf", dpi=800, bbox_inches="tight")
