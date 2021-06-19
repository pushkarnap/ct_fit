"""
AUTHOR: Paarangat Pushkarna
DATE: 26/05/2021
MODIFIED: 30/05/2021
"""
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import numpy as np
from lmfit import fit_report

GRID_Y = 4
GRID_X = 1

"""
=======================================
Various plotting functions for reporting on the fitting process.
=======================================
"""
def plot_reconstruction(energy, intensity, uncertainty, title, name):
    fig = plt.figure(constrained_layout = True)
    grid = fig.add_gridspec(GRID_Y, GRID_X)

    fig_ax1 = fig.add_subplot(grid[0:-1, :])
    fig_ax1.set_title(title)
    fig_ax1.set_xlabel("Energy [eV]")
    fig_ax1.set_ylabel("Intensity [arb. units]")
    fig_ax1.plot(energy, intensity, color = "black", linewidth = 0.75)

    fig_ax2 = fig.add_subplot(grid[-1, :])
    fig_ax2.set_title("Uncertainty")
    fig_ax2.set_xlabel("Energy [eV]")
    fig_ax2.plot(energy, uncertainty, color = "blue", linewidth = 0.75)

    fig.savefig(name)
    return

def plot_sticks(spectra, chosen, colours, xmin, xmax, title, name):
    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.set_ylabel(r"Line strength [arb. units]")
    ax.set_xlabel("Energy [eV]")
    ax.set_xlim([xmin, xmax])

    for i in range(len(chosen)):
        S = spectra.spectra[chosen[i]].line_strength
        E = spectra.spectra[chosen[i]].energies

        markerline, stemline, baseline = ax.stem(E, S, use_line_collection = True,
                                                    linefmt = colours[i], label = f"{chosen[i]}")
        markerline.set_markersize(0)
        stemline.set_linewidths(np.array([0.5]))
        baseline.set_linewidth(0)
        ax.set_ylim(0)
        ax.legend()

    fig.savefig(name)

    return

def plot_optimal(energy, intensity_expt, intensity_thy, residual, chi_red, title, name):
    fig, ax = plt.subplots()
    ax.plot(energy, intensity_expt,
        color = "black", label = "experiment", linewidth = 0.75)
    ax.plot(energy, intensity_thy,
        color = "blue", label = "theory", linewidth = 0.75)
    ax.set_title(title)
    ax.set_ylabel("Relative Intensity [arb. units]")
    ax.set_xlabel("Energy [eV]")
    ax.text(4450, 0.95, f"$\chi_r^2 = {chi_red:.2f}$")
    ax.legend()
    fig.savefig(name)

    #Plotting code in case residuals are also required.
    """
    fig = plt.figure(constrained_layout = True)
    grid = fig.add_gridspec(GRID_Y, GRID_X)

    fig_ax1 = fig.add_subplot(grid[:, :])
    fig_ax1.set_title(title)
    fig_ax1.set_xlabel("Energy [eV]")
    fig_ax1.set_ylabel("Intensity [arb. units]")
    fig_ax1.plot(energy, intensity_thy,
                    color = "blue", label = "theory", linewidth = 0.75)
    fig_ax1.plot(energy, intensity_expt,
                    color = "black", label = "experiment", linewidth = 0.75)
    fig_ax1.text(4450, 0.9, f"$\chi_r^2 = {chi_red:.2f}$")
    fig_ax1.legend()


    fig_ax2 = fig.add_subplot(grid[-1, :])
    fig_ax2.set_title("Residual")
    fig_ax2.set_xlabel("Energy [eV]")
    fig_ax2.plot(energy, residual, color = "green", linewidth = 0.75)
    """
    return
