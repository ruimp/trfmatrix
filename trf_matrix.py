#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
_ = np.newaxis
plt.style.use("seaborn-deep")

class TransferMatrix:
    def __init__(self, temperature, chem_potential, wavelength_list):
        self.wl_list = np.array(wavelength_list)
        self.n_wls = len(wavelength_list)
        self.chem_potential = chem_potential
        self.temperature = temperature
        self.R = 0
        self.T = 0
        self.A = 0
        self.trf_matrix = np.diag(np.ones(2, dtype=complex))
        self.xi = np.zeros(len(wavelength_list), dtype=complex)

    def get_graphene_xi(self):
        hbar = 1.0545e-34
        kB = 1.3806e-23
        t = 2.7 * 1.602e-19
        e = 1.602e-19
        mu0 = np.pi * 4e-7
        c = 3.0e8
        omega = 2.0 * c * np.pi / (self.wl_list * 1e-9)
        mu = self.chem_potential * e
        sigma_gr_r = np.zeros(self.n_wls, dtype=complex)
        sigma_gr_i = np.zeros(self.n_wls, dtype=complex)
        sigma_gr_r += np.pi * e * e / (2.0*hbar)
        sigma_gr_r *= (0.5 + (hbar*omega)**2/(72*t*t))
        sigma_gr_r *=(np.tanh((hbar*omega + 2.0*mu)/(4.0*kB*self.temperature)) + np.tanh((hbar*omega - 2.0*mu) /
                                                                (4.0*kB*self.temperature)))
        sigma_gr_i += (2.0*mu*e*e/(hbar*hbar*omega))*(1 - (2.0*mu*mu/(9.0*t*t)))
        sigma_gr_i -= (0.5*e*e/hbar)*np.log(np.abs((hbar*omega + 2.0*mu)/(hbar*omega - 2.0*mu)))
        sigma_gr_i -= ((1/72.0)*e*e/hbar)*((hbar*omega/t)**2)*np.log( np.abs((hbar*omega + 2.0*mu) /
                                                                (hbar*omega - 2.0*mu)))
        self.sigma_gr = sigma_gr_r + 1j*sigma_gr_i
        self.xi = mu0 * c * self.sigma_gr

    def get_constant_xi(self, xi):
        self.xi = np.ones(self.n_wls, dtype=complex) * xi

    def init_trf_matrix(self):
        self.trf_matrix = np.diag(np.ones(2, dtype=complex))

    def multiply_by_layer(self, n1, n2, wl_ind):
        A = (n1 / n2) * (1.0 - self.xi[wl_ind] / n1)
        B = (n1 / n2) * (1.0 + self.xi[wl_ind] / n1)
        t_layer = 0.5 * np.array([[1.0 + A, 1.0 - B], [1.0 - A, 1.0 + B]], dtype=complex)
        self.trf_matrix = np.dot(t_layer, self.trf_matrix)

    def multiply_by_chunk(self, n1, n2, d, wl_ind):
        A = (n1 / n2) * (1.0 - self.xi[wl_ind] / n1)
        B = (n1 / n2) * (1.0 + self.xi[wl_ind] / n1)
        EXP = np.exp(2j * np.pi * d * n1 / self.wl_list[wl_ind])
        t_chunk = 0.5*np.array([[(1.0 + A)*EXP, (1.0 - B)*(EXP.conjugate())], [(1.0 - A)*EXP,
                                                                (1.0 + B)*(EXP.conjugate())]], dtype=complex)
        self.trf_matrix = np.dot(t_chunk, self.trf_matrix)

    def get_R(self):
        r = -self.trf_matrix[1, 0] / self.trf_matrix[1, 1]
        self.R = np.absolute(r)**2

    def get_T(self, n_in, n_out):
        t = self.trf_matrix[0, 0] * self.trf_matrix[1, 1] - (self.trf_matrix[1 ,0] * self.trf_matrix[0, 1])
        a = np.absolute(t)
        b = np.absolute(self.trf_matrix[1, 1])
        self.T = ((n_out / n_in)*(a/b))**2

    def get_A(self):
        self.A = 1.0 - self.R - self.T

    def get_coeffs(self, n_in, n_out):
        self.get_R()
        self.get_T(n_in, n_out)
        self.get_A()


    """R vs WL - 2 layer"""

    def plot_r_vs_wl_2layer(self, n_in, n_pol, n_out, d_list):
        n_plots = len(d_list)
        x = np.tile(self.wl_list, (n_plots, 1))
        y = np.zeros((n_plots, self.n_wls))
        labels = ["D = {} nm".format(d) for d in d_list]
        for i, d in enumerate(d_list):
            for j in range(self.n_wls):
                self.init_trf_matrix()
                self.multiply_by_layer(n_in, n_pol, j)
                self.multiply_by_chunk(n_pol, n_out, d, j)
                self.get_R()
                y[i, j] = 100*self.R
        self.plot_routine(x, y, labels, xl="Wavelength (nm)", yl="Reflectance")

    def plot_t_vs_wl_2layer(self, n_in, n_pol, n_out, d_list):
        n_plots = len(d_list)
        x = np.tile(self.wl_list, (n_plots, 1))
        y = np.zeros((n_plots, self.n_wls))
        labels = ["D = {} nm".format(d) for d in d_list]
        for i, d in enumerate(d_list):
            for j in range(self.n_wls):
                self.init_trf_matrix()
                self.multiply_by_layer(n_in, n_pol, j)
                self.multiply_by_chunk(n_pol, n_out, d, j)
                self.get_T(n_in, n_out)
                y[i, j] = 100*self.T
        self.plot_routine(x, y, labels, xl="Wavelength (nm)", yl="Transmittance")

    def plot_a_vs_wl_2layer(self, n_in, n_pol, n_out, d_list):
        n_plots = len(d_list)
        x = np.tile(self.wl_list, (n_plots, 1))
        y = np.zeros((n_plots, self.n_wls))
        labels = ["D = {} nm".format(d) for d in d_list]
        for i, d in enumerate(d_list):
            for j in range(self.n_wls):
                self.init_trf_matrix()
                self.multiply_by_layer(n_in, n_pol, j)
                self.multiply_by_chunk(n_pol, n_out, d, j)
                self.get_coeffs(n_in, n_out)
                y[i, j] = 100*self.A
        self.plot_routine(x, y, labels, xl="Wavelength (nm)", yl="Absorption")


    """R vs WL - multi-layer"""

    def plot_r_vs_wl(self, n_in, n_pol_list, n_out, d_list):
        n_plots = len(d_list)
        n_layer = len(n_pol_list)
        x = np.tile(self.wl_list, (n_plots, 1))
        y = np.zeros((n_plots, self.n_wls))
        labels = ["D = {} nm".format(d) for d in d_list]
        for i, d in enumerate(d_list):
            for j in range(self.n_wls):
                self.init_trf_matrix()
                self.multiply_by_layer(n_in, n_pol_list[0], j)
                for k in range(1, n_layer):
                    self.multiply_by_chunk(n_pol_list[k-1], n_pol_list[k], d, j)
                self.multiply_by_chunk(n_pol_list[-1], n_out, d, j)
                self.get_R()
                y[i, j] = 100*self.R
        self.plot_routine(x, y, labels, xl="Wavelength (nm)", yl="Reflectance")

    def plot_t_vs_wl(self, n_in, n_pol_list, n_out, d_list):
        n_plots = len(d_list)
        n_layer = len(n_pol_list)
        x = np.tile(self.wl_list, (n_plots, 1))
        y = np.zeros((n_plots, self.n_wls))
        labels = ["D = {} nm".format(d) for d in d_list]
        for i, d in enumerate(d_list):
            for j in range(self.n_wls):
                self.init_trf_matrix()
                self.multiply_by_layer(n_in, n_pol_list[0], j)
                for k in range(1, n_layer):
                    self.multiply_by_chunk(n_pol_list[k-1], n_pol_list[k], d, j)
                self.multiply_by_chunk(n_pol_list[-1], n_out, d, j)
                self.get_T(n_in, n_out)
                y[i, j] = 100*self.T
        self.plot_routine(x, y, labels, xl="Wavelength (nm)", yl="Transmittance")

    def plot_a_vs_wl(self, n_in, n_pol_list, n_out, d_list):
        n_plots = len(d_list)
        n_layer = len(n_pol_list)
        x = np.tile(self.wl_list, (n_plots, 1))
        y = np.zeros((n_plots, self.n_wls))
        labels = ["D = {} nm".format(d) for d in d_list]
        for i, d in enumerate(d_list):
            for j in range(self.n_wls):
                self.init_trf_matrix()
                self.multiply_by_layer(n_in, n_pol_list[0], j)
                for k in range(1, n_layer):
                    self.multiply_by_chunk(n_pol_list[k-1], n_pol_list[k], d, j)
                self.multiply_by_chunk(n_pol_list[-1], n_out, d, j)
                self.get_coeffs(n_in, n_out)
                y[i, j] = 100*self.A
        self.plot_routine(x, y, labels, xl="Wavelength (nm)", yl="Absorption")


    """Plotting routine"""

    def plot_routine(self, xs, ys, labels, xl=None, yl=None):
        h = len(ys) // 2
        fig, axs = plt.subplots(h, 2, sharex=True, sharey=True, gridspec_kw={'hspace':0, 'wspace':0.03})
        for i, (x, y, label) in enumerate(zip(xs, ys, labels)):
            axs[i//2, i%2].plot(x, y, label=label)
            axs[i//2, i%2].legend()
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
        plt.xlabel(xl)
        plt.ylabel(yl, labelpad=10)
