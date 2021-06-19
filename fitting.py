"""
AUTHOR: Paarangat Pushkarna
DATE: 26/05/2021
MODIFIED: 30/05/2021
"""
from ct_process import spectra, spectrum, characterisation
from lmfit import Minimizer
from time import perf_counter
from plotting import *
import numpy as np

"""
=======================================
Functions which assist fitting routines
=======================================
"""

def use_characterisation(char_path, energies):
    """
    Run the '.reconstruct' and '.uncertainties' methods to compute deconvolved
    experimental spectrum and associated uncertainties.
    Inputs
    -----------------------------------
    char_path: path to csv file containing characterisation information. str.
    energies: energy values as defined by the chosen energy range. np array.
    Returns
    -----------------------------------
    intensities_expt: intensity values reconstructed from charactersiation.
    np.array
    uncertainties_expt: uncertainty values reconstructed from characterisation.
    np.array
    """
    print("LOADING CHARACTERISATION...")
    expt = characterisation(char_path)
    intensities_expt = expt.reconstruct(energies)
    uncertainties_expt = expt.uncertainties(energies)
    print("CHARACTERISATION LOADED")
    return (intensities_expt, uncertainties_expt)

def run_fit(spectra, data, uncert, energies):
    """
    Run the fitting algorithm for a chosen set of ct files.
    Inputs:
    -----------------------------------
    spectra: spectra object containing information on all ct files. spectra obj.
    data: intensity values reconstructed from experiment. np.array.
    uncert: uncertainties in intensity computed from experiment. np.array.
    energies: energy values defined by chosen energy range. np.array.
    Returns
    -----------------------------------
    result: lmfit MinimizerResult object which reports the optimal fit
    parameters.
    """
    thy = spectra
    params = thy.create_parameters()
    minner = Minimizer(thy.objective, params, fcn_args = (energies, data, uncert))
    start_fit = perf_counter()
    result = minner.minimize() #fitting algorithm executed here.
    end_fit = perf_counter()
    print(f"TOOK: {end_fit - start_fit} seconds")
    return result

def auto_fit(ct_path, recon_path, res, folder_name, window):
    """
    Run automated fitting routine that fits every combination of ct files.
    Inputs
    -----------------------------------
    ct_path: path to ct file storage directory. str.
    recon_path: path to csv file storing experimental charactersiation. str.
    res: number of points over which to run fitting algorithm
    folder_name: folder in which to store output pngs
    window: 'padding' around the minimum and maximum transition energies
    """
    print("BEGIN FITTING ROUTINE")
    print("LOADING CTs")
    thy = spectra(ct_path)
    min, max = thy.extrema
    ENERGY_MIN = min - window
    ENERGY_MAX = max + window
    ENERGIES = np.linspace(ENERGY_MIN, ENERGY_MAX, res)
    fit_combos = thy.routine
    print("DONE LOADING CTs")

    inten_expt, uncert_expt = use_characterisation(recon_path, ENERGIES)

    fout = open("chis.txt", 'w')
    fout.write("satellites, chi sq. red.\n")

    for fit_combo in fit_combos:
        print(f"NOW FITTING: {fit_combo}")
        thy.chosen = fit_combo
        fit_out = run_fit(thy, inten_expt, uncert_expt, ENERGIES)
        opt_params = fit_out.params
        chi_red = fit_out.redchi
        residual = fit_out.residual
        intensity_thy_opt = thy.model(opt_params, ENERGIES)
        
        fit_plot_name = ""
        fit_plot_title = ""
        for name in fit_combo:
            name_split = name.split(sep = "_", maxsplit = -1)
            fit_plot_name += name_split[2]
            fit_plot_title += name_split[2]
        fit_plot_name += ".png"

        fout.write(f"{fit_plot_title}, {chi_red:.4f}\n")

        #plot_optimal(ENERGIES, inten_expt, intensity_thy_opt, residual,
                                #chi_red, fit_plot_title, fit_plot_name)
        print(f"DONE FITTING: {fit_combo}")
    fout.close()
    print("FITTING COMPLETE!")
    return
