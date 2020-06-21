#!/usr/bin/env python3

import trf_matrix as trf
import numpy as np
import matplotlib.pyplot as plt

n_in = 1.46
n_pol = 1.6
n_out = 1.0
d_list = [10, 20, 30, 40]
wls = np.linspace(40, 2000, 1000)
T = trf.TransferMatrix(300, -0.2, wls)

"""R vs WL - 2 layer"""

#T.get_constant_xi(0.149)
#T.plot_r_vs_wl_2layer(n_in, n_pol, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/R_vs_WL_Xic_2layer.pdf", dpi=1000, bbox_inches="tight")

#T.get_graphene_xi()
#T.plot_r_vs_wl_2layer(n_in, n_pol, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/R_vs_WL_Xig_2layer.pdf", dpi=1000, bbox_inches="tight")

"""T vs WL - 2 layer"""

#T.get_constant_xi(0.149)
#T.plot_t_vs_wl_2layer(n_in, n_pol, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/T_vs_WL_Xic_2layer.pdf", dpi=1000, bbox_inches="tight")

#T.get_graphene_xi()
#T.plot_t_vs_wl_2layer(n_in, n_pol, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/T_vs_WL_Xig_2layer.pdf", dpi=1000, bbox_inches="tight")

"""A vs WL - 2 layer"""

#T.get_constant_xi(0.149)
#T.plot_a_vs_wl_2layer(n_in, n_pol, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/A_vs_WL_Xic_2layer.pdf", dpi=1000, bbox_inches="tight")

#T.get_graphene_xi()
#T.plot_a_vs_wl_2layer(n_in, n_pol, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/A_vs_WL_Xig_2layer.pdf", dpi=1000, bbox_inches="tight")


"""R vs WL - 3 layer"""

n_3layer = np.ones(2) * n_pol

#T.get_constant_xi(0.149)
#T.plot_r_vs_wl(n_in, n_3layer, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/R_vs_WL_Xic_3layer.pdf", dpi=1000, bbox_inches="tight")

#T.get_graphene_xi()
#T.plot_r_vs_wl(n_in, n_3layer, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/R_vs_WL_Xig_3layer.pdf", dpi=1000, bbox_inches="tight")

"""T vs WL - 3 layer"""

#T.get_constant_xi(0.149)
#T.plot_t_vs_wl(n_in, n_3layer, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/T_vs_WL_Xic_3layer.pdf", dpi=1000, bbox_inches="tight")

#T.get_graphene_xi()
#T.plot_t_vs_wl(n_in, n_3layer, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/T_vs_WL_Xig_3layer.pdf", dpi=1000, bbox_inches="tight")

"""A vs WL - 3 layer"""

#T.get_constant_xi(0.149)
#T.plot_a_vs_wl(n_in, n_3layer, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/A_vs_WL_Xic_3layer.pdf", dpi=1000, bbox_inches="tight")

#T.get_graphene_xi()
#T.plot_a_vs_wl(n_in, n_3layer, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/A_vs_WL_Xig_3layer.pdf", dpi=1000, bbox_inches="tight")


"""R vs WL - 10 layer"""

n_10layer = np.ones(9) * n_pol

#T.get_constant_xi(0.149)
#T.plot_r_vs_wl(n_in, n_10layer, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/R_vs_WL_Xic_10layer.jpg", dpi=1000, bbox_inches="tight")

#T.get_graphene_xi()
#T.plot_r_vs_wl(n_in, n_3layer, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/R_vs_WL_Xig_10layer.pdf", dpi=1000, bbox_inches="tight")

"""T vs WL - 10 layer"""

#T.get_constant_xi(0.149)
#T.plot_t_vs_wl(n_in, n_3layer, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/T_vs_WL_Xic_10layer.pdf", dpi=1000, bbox_inches="tight")

#T.get_graphene_xi()
#T.plot_t_vs_wl(n_in, n_3layer, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/T_vs_WL_Xig_10layer.pdf", dpi=1000, bbox_inches="tight")

"""A vs WL - 10 layer"""

#T.get_constant_xi(0.149)
#T.plot_a_vs_wl(n_in, n_3layer, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/A_vs_WL_Xic_10layer.pdf", dpi=1000, bbox_inches="tight")

#T.get_graphene_xi()
#T.plot_a_vs_wl(n_in, n_3layer, n_out, d_list)
#plt.savefig("Plots/Plot_vs_Wavelength/A_vs_WL_Xig_10layer.pdf", dpi=1000, bbox_inches="tight")


plt.show()
