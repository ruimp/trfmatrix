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
ds = np.linspace(10, 40, 200)
wls = np.linspace(400, 2000, 200)
T = trf.TransferMatrix(300, -0.2, wls)

"""Reflectance"""

#T.get_constant_xi(xi)
#T.plot_r_heatmap(n_in, n_3layer, n_out, ds)
#plt.savefig("Plots/Plot_Heatmap/R_heatmap_Xic_3layer.pdf", dpi=800, bbox_inches="tight")

#T.plot_r_heatmap(n_in, n_10layer, n_out, ds)
#plt.savefig("Plots/Plot_Heatmap/R_heatmap_Xic_10layer.pdf", dpi=800, bbox_inches="tight")

#T.plot_r_heatmap(n_in, n_20layer, n_out, ds)
#plt.savefig("Plots/Plot_Heatmap/R_heatmap_Xic_20layer.pdf", dpi=800, bbox_inches="tight")

#T.get_graphene_xi()
#T.plot_r_heatmap(n_in, n_3layer, n_out, ds)
#plt.savefig("Plots/Plot_Heatmap/R_heatmap_Xig_3layer.pdf", dpi=800, bbox_inches="tight")

#T.plot_r_heatmap(n_in, n_10layer, n_out, ds)
#plt.savefig("Plots/Plot_Heatmap/R_heatmap_Xig_10layer.pdf", dpi=800, bbox_inches="tight")

#T.plot_r_heatmap(n_in, n_20layer, n_out, ds)
#plt.savefig("Plots/Plot_Heatmap/R_heatmap_Xig_20layer.pdf", dpi=800, bbox_inches="tight")

"""Transmittance"""

T.get_constant_xi(xi)
T.plot_t_heatmap(n_in, n_3layer, n_out, ds)
plt.savefig("Plots/Plot_Heatmap/T_heatmap_Xic_3layer.pdf", dpi=800, bbox_inches="tight")

T.plot_t_heatmap(n_in, n_10layer, n_out, ds)
plt.savefig("Plots/Plot_Heatmap/T_heatmap_Xic_10layer.pdf", dpi=800, bbox_inches="tight")

T.plot_t_heatmap(n_in, n_20layer, n_out, ds)
plt.savefig("Plots/Plot_Heatmap/T_heatmap_Xic_20layer.pdf", dpi=800, bbox_inches="tight")

T.get_graphene_xi()
T.plot_t_heatmap(n_in, n_3layer, n_out, ds)
plt.savefig("Plots/Plot_Heatmap/T_heatmap_Xig_3layer.pdf", dpi=800, bbox_inches="tight")

T.plot_t_heatmap(n_in, n_10layer, n_out, ds)
plt.savefig("Plots/Plot_Heatmap/T_heatmap_Xig_10layer.pdf", dpi=800, bbox_inches="tight")

T.plot_t_heatmap(n_in, n_20layer, n_out, ds)
plt.savefig("Plots/Plot_Heatmap/T_heatmap_Xig_20layer.pdf", dpi=800, bbox_inches="tight")

"""Absorption"""

T.get_constant_xi(xi)
T.plot_a_heatmap(n_in, n_3layer, n_out, ds)
plt.savefig("Plots/Plot_Heatmap/A_heatmap_Xic_3layer.pdf", dpi=800, bbox_inches="tight")

T.plot_a_heatmap(n_in, n_10layer, n_out, ds)
plt.savefig("Plots/Plot_Heatmap/A_heatmap_Xic_10layer.pdf", dpi=800, bbox_inches="tight")

T.plot_a_heatmap(n_in, n_20layer, n_out, ds)
plt.savefig("Plots/Plot_Heatmap/A_heatmap_Xic_20layer.pdf", dpi=800, bbox_inches="tight")

T.get_graphene_xi()
T.plot_a_heatmap(n_in, n_3layer, n_out, ds)
plt.savefig("Plots/Plot_Heatmap/A_heatmap_Xig_3layer.pdf", dpi=800, bbox_inches="tight")

T.plot_a_heatmap(n_in, n_10layer, n_out, ds)
plt.savefig("Plots/Plot_Heatmap/A_heatmap_Xig_10layer.pdf", dpi=800, bbox_inches="tight")

T.plot_a_heatmap(n_in, n_20layer, n_out, ds)
plt.savefig("Plots/Plot_Heatmap/A_heatmap_Xig_20layer.pdf", dpi=800, bbox_inches="tight")
