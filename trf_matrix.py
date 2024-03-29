#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
plt.style.use("seaborn-white")

# Enable for data export
#os.chdir("/home/rmp/Projects/Transfer_Matrix/Data")

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
        a = (n1 / n2) * (1 - self.xi[wl_ind] / n1)
        b = (n1 / n2) * (1 + self.xi[wl_ind] / n1)
        t_layer = .5 * np.array(
            [[1 + a, 1 - b],
            [1 - a, 1 + b]])
        self.trf_matrix = np.dot(t_layer, self.trf_matrix)

    def multiply_by_chunk(self, n1, n2, d, wl_ind):
        a = (n1 / n2) * (1 - self.xi[wl_ind] / n1)
        b = (n1 / n2) * (1 + self.xi[wl_ind] / n1)
        c = np.exp(2j * np.pi * d * n1 / self.wl_list[wl_ind])
        t_chunk = .5 * np.array([
            [(1 + a) * c, (1 - b) / c],
            [(1 - a) * c, (1 + b) / c]])
        self.trf_matrix = np.dot(t_chunk, self.trf_matrix)

    def get_r(self):
        return -self.trf_matrix[1, 0] / self.trf_matrix[1, 1]

    def get_R(self):
        r = -self.trf_matrix[1, 0] / self.trf_matrix[1, 1]
        self.R = np.absolute(r)**2

    def get_T(self, n_in, n_out):
        t = self.trf_matrix[0, 0] * self.trf_matrix[1, 1] - (self.trf_matrix[1 ,0] * self.trf_matrix[0, 1])
        a = np.absolute(t)
        b = np.absolute(self.trf_matrix[1, 1])
        self.T = ((n_out / n_in)*(a/b)**2)

    def get_A(self):
        self.A = 1.0 - self.R - self.T

    def get_coeffs(self, n_in, n_out):
        self.get_R()
        self.get_T(n_in, n_out)
        self.get_A()

    def trf_routine(self, n_in, n_pol_list, n_out, d_list, wl_ind):
        """
        n_pol_list      size: n_layers - 1
        d_list          size: n_layers
        """
        self.init_trf_matrix()
        self.multiply_by_layer(n_in, n_pol_list[0], wl_ind)
        for i in range(1, n_pol_list.size):
            self.multiply_by_chunk(n_pol_list[i-1], n_pol_list[i], d_list[i-1], wl_ind)
        self.multiply_by_chunk(n_pol_list[-1], n_out, d_list[-1], wl_ind)

    def trf_routine_sym(self, n_in, n_pol_list, n_out, d, n_interface, D, wl_ind):
        self.init_trf_matrix()
        self.multiply_by_layer(n_in, n_pol_list[0], wl_ind)
        for i in range(1, n_pol_list.size):
            self.multiply_by_chunk(n_pol_list[i-1], n_pol_list[i], d, wl_ind)
        self.multiply_by_chunk(n_pol_list[-1], n_interface, d, wl_ind)
        self.multiply_by_chunk(n_interface, n_pol_list[-1], D, wl_ind)
        for i in range(1, n_pol_list.size):
            self.multiply_by_chunk(n_pol_list[-i], n_pol_list[-i-1], d, wl_ind)
        self.multiply_by_chunk(n_pol_list[0], n_out, d, wl_ind)

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
        x = np.tile(self.wl_list, (n_plots, 1))
        y = np.zeros((n_plots, self.n_wls))
        labels = ["D = {} nm".format(d) for d in d_list]
        for i, d in enumerate(d_list):
            for j in range(self.n_wls):
                self.trf_routine(n_in, n_pol_list, n_out, d, j)
                self.get_R()
                y[i, j] = 100*self.R
    #        np.savetxt("out_python.dat", np.array([x[1], y[1]]).transpose())
        self.plot_routine(x, y, labels, xl="Wavelength (nm)", yl="Reflectance")

    def plot_t_vs_wl(self, n_in, n_pol_list, n_out, d_list):
        n_plots = len(d_list)
        x = np.tile(self.wl_list, (n_plots, 1))
        y = np.zeros((n_plots, self.n_wls))
        labels = ["D = {} nm".format(d) for d in d_list]
        for i, d in enumerate(d_list):
            for j in range(self.n_wls):
                self.trf_routine(n_in, n_pol_list, n_out, d, j)
                self.get_T(n_in, n_out)
                y[i, j] = 100*self.T
        self.plot_routine(x, y, labels, xl="Wavelength (nm)", yl="Transmittance")

    def plot_a_vs_wl(self, n_in, n_pol_list, n_out, d_list):
        n_plots = len(d_list)
        x = np.tile(self.wl_list, (n_plots, 1))
        y = np.zeros((n_plots, self.n_wls))
        labels = ["D = {} nm".format(d) for d in d_list]
        for i, d in enumerate(d_list):
            for j in range(self.n_wls):
                self.trf_routine(n_in, n_pol_list, n_out, d, j)
                self.get_coeffs(n_in, n_out)
                y[i, j] = 100*self.A
        self.plot_routine(x, y, labels, xl="Wavelength (nm)", yl="Absorption")


    """Heatmap plots"""

    def plot_r_heatmap(self, n_in, n_pol_list, n_out, d_list):
        R = np.zeros((len(d_list), self.n_wls))
        for i, d in enumerate(d_list):
            for j in range(self.n_wls):
                self.trf_routine(n_in, n_pol_list, n_out, d, j)
                self.get_R()
                R[i, j] = 100*self.R
        fig, ax = plt.subplots()
        c = ax.pcolormesh(self.wl_list, d_list, R)
        fig.colorbar(c, ax=ax)
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Polymer Width (nm)")

    def plot_t_heatmap(self, n_in, n_pol_list, n_out, d_list):
        R = np.zeros((len(d_list), self.n_wls))
        for i, d in enumerate(d_list):
            for j in range(self.n_wls):
                self.trf_routine(n_in, n_pol_list, n_out, d, j)
                self.get_T(n_in, n_out)
                R[i, j] = 100*self.T
        fig, ax = plt.subplots()
        c = ax.pcolormesh(self.wl_list, d_list, R)
        fig.colorbar(c, ax=ax)
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Polymer Width (nm)")

    def plot_a_heatmap(self, n_in, n_pol_list, n_out, d_list):
        R = np.zeros((len(d_list), self.n_wls))
        for i, d in enumerate(d_list):
            for j in range(self.n_wls):
                self.trf_routine(n_in, n_pol_list, n_out, d, j)
                self.get_coeffs(n_in, n_out)
                R[i, j] = 100*self.A
        fig, ax = plt.subplots()
        c = ax.pcolormesh(self.wl_list, d_list, R)
        fig.colorbar(c, ax=ax)
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Polymer Width (nm)")


    """Field Reconstruction"""

    def get_field(self, n_in, n_pol_list, n_out, d_list, wl_ind):
        self.trf_routine(n_in, n_pol_list, n_out, d_list, wl_ind)
        n_layers = n_pol_list.size
        v = np.zeros((n_layers + 1, 2), dtype=complex)
        v[0] = np.array([1, self.get_r()])
        self.init_trf_matrix()
        self.multiply_by_layer(n_in, n_pol_list[0], wl_ind)
        v[0] = self.trf_matrix.dot(v[0])
        for i in range(1, n_layers):
            self.init_trf_matrix()
            self.multiply_by_chunk(n_pol_list[i-1], n_pol_list[i], d_list[i-1], wl_ind)
            v[i] = self.trf_matrix.dot(v[i-1])
        self.init_trf_matrix()
        self.multiply_by_chunk(n_pol_list[-1], n_out, d_list[-1], wl_ind)
        v[-1] = self.trf_matrix.dot(v[-2])
        return v.sum(axis=-1)

    def plot_field_vs_d(self, n_in, n_pol_list, n_out, d_list, wl_ind = 0):
        n_plots = len(d_list)
        n_layers = n_pol_list.size
        xs = np.tile(np.arange(n_layers + 1), (n_plots, 1))
        ys = np.zeros((n_plots, n_layers + 1))
        labels = ["D = {} nm".format(d) for d in d_list]
        for i, d in enumerate(d_list):
            E = self.get_field(n_in, n_pol_list, n_out, d, wl_ind)
            ys[i] = np.absolute(E)
        self.plot_routine(xs, ys, labels, xl="Layers", yl=r"$|E|$ - Absolute electric Field",
                            xticks=np.arange(0, n_layers + 2, 5))

    def plot_field_vs_wl(self, n_in, n_pol_list, n_out, d):
        n_layers = n_pol_list.size
        xs = np.tile(np.arange(n_layers + 1), (self.n_wls, 1))
        ys = np.zeros((self.n_wls, n_layers + 1))
        labels = [r"$\lambda$ = {} nm".format(wl) for wl in self.wl_list]
        for i in range(self.n_wls):
            E = self.get_field(n_in, n_pol_list, n_out, d, i)
            ys[i] = np.absolute(E)
        self.plot_routine(xs, ys, labels, xl="Layers", yl=r"$|E|$ - Absolute Electric Field",
                            xticks=np.arange(0,n_layers + 2, 5),)

    def get_field_sum(self, n_in, n_pol_list, n_out, d, wl_ind):
        E = np.absolute(self.get_field(n_in, n_pol_list, n_out, d, wl_ind))
        return np.sum(E**2)

    def plot_field_sum_vs_layers(self, n_in, n_pol_list, n_out, d):
        n_layers = n_pol_list.size
        xs = np.tile(np.arange(n_layers), (self.n_wls, 1))
        ys = np.zeros((self.n_wls, n_layers))
        labels = [r"$\lambda$ = {} nm".format(wl) for wl in self.wl_list]
        for i in range(self.n_wls):
            for j in range(1, n_layers):
                ys[i, j] = self.get_field_sum(n_in, n_pol_list[:j], n_out, d, i)
        self.plot_routine(xs, ys, labels, xl="Layers", yl=r"$\sum |E|^2$ - Sum of eletric field squared",
                            xticks=np.arange(0,n_layers + 2, 20))
        return xs, ys

    def get_a_list(self, n_in, n_pol_list, n_out, d):
        n_layers = n_pol_list.size
        ys = np.zeros((self.n_wls, n_layers))
        for i in range(self.n_wls):
            for j in range(1, n_layers):
                self.trf_routine(n_in, n_pol_list[:j], n_out, d, i)
                self.get_coeffs(n_in, n_out)
                ys[i, j] = self.A
        return ys

    def plot_r_vs_layers(self, n_in, n_pol_list, n_out, d):
        n_layers = n_pol_list.size
        xs = np.tile(np.arange(n_layers), (self.n_wls, 1))
        ys = np.zeros((self.n_wls, n_layers))
        labels = [r"$\lambda$ = {} nm".format(wl) for wl in self.wl_list]
        for i in range(self.n_wls):
            for j in range(1, n_layers):
                self.trf_routine(n_in, n_pol_list[:j], n_out, d, i)
                self.get_coeffs(n_in, n_out)
                ys[i, j] = self.R*100
        self.plot_routine(xs, ys, labels, xl="Layers", yl=r"Reflectance",
                            xticks=np.arange(0,n_layers + 2, 1))


    def plot_a_vs_layers(self, n_in, n_pol_list, n_out, d):
        n_layers = n_pol_list.size
        xs = np.tile(np.arange(n_layers), (self.n_wls, 1))
        ys = np.zeros((self.n_wls, n_layers))
        labels = [r"$\lambda$ = {} nm".format(wl) for wl in self.wl_list]
        for i in range(self.n_wls):
            for j in range(1, n_layers):
                self.trf_routine(n_in, n_pol_list[:j], n_out, d, i)
                self.get_coeffs(n_in, n_out)
                ys[i, j] = self.A
        self.plot_routine(xs, ys, labels, xl="Layers", yl=r"Absorption",
                            xticks=np.arange(0,n_layers + 2, 20))


    def plot_a_vs_n_out(self, n_in, n_pol_list, n_out_list, d):
        xs = np.tile(np.arange(len(n_out_list)), (self.n_wls, 1))
        ys = np.zeros((self.n_wls, len(n_out_list)))
        labels = [r"$\lambda$ = {} nm".format(wl) for wl in self.wl_list]
        for i in range(self.n_wls):
            for j, n_out in enumerate(n_out_list):
                self.trf_routine(n_in, n_pol_list, n_out, d, i)
                self.get_coeffs(n_in, n_out)
                ys[i, j] = self.A
        self.plot_routine(xs, ys, labels, xl="n_out", yl=r"Absorption")

    """Plotting routine"""

    def plot_routine(self, xs, ys, labels, xl=None, yl=None, xticks = None):
        h = len(ys) // 2
        fig, axs = plt.subplots(h, 2, sharex=True, sharey=True, gridspec_kw={'hspace':0, 'wspace':0.03})
        for i, (x, y, label) in enumerate(zip(xs, ys, labels)):
            axs[i//2, i%2].plot(x, y, label=label)
            axs[i//2, i%2].legend()
        plt.xticks(ticks=xticks)
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
        plt.xlabel(xl)
        plt.ylabel(yl, labelpad=10)
