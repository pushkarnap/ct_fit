"""
AUTHOR: Paarangat Pushkarna
DATE: 26/05/2021
MODIFIED: 30/05/2021
"""
import os
from fitting import *
from plotting import *
from ct_process import spectra, spectrum, characterisation
import numpy as np
import matplotlib.pyplot as plt
from lmfit import fit_report

RES = int(1e4)
WINDOW = 2

###USER_CHOICE###
char_plot_title = r""
char_plot_name = "" #Reconstruction of experiment
stick_plot_title = r""#
stick_plot_name = ""#Visualisation of theory
fit_plot_title = r""#
fit_plot_name = ""#Visualisation of fit
######
def main():
    work_path = os.getcwd()
    ###PUT CHARACTERISATION .HERE
    to_reconstruct = os.path.join(work_path, "")
    ###USER_CHOICE###
    ct_folder_name = "" ###Choose folder where ct files are stored
    ######
    cts = os.path.join(work_path, ct_folder_name)
    stick_plot_path = os.path.join(cts, stick_plot_name)
    fit_plot_path = os.path.join(cts, fit_plot_name)

    #create spectra object to store data from ct files
    thy = spectra(cts)
    ###USER_CHOICE###
    to_fit = []
    ###.ct files to be included in fit
    ######
    thy.chosen = to_fit

    min, max = thy.extrema
    #prepare the energy range to be considered on the basis of smallest
    #and largest transition energies present in the theory
    ENERGY_MIN = min-WINDOW
    ENERGY_MAX = max+WINDOW
    ENERGIES = np.linspace(ENERGY_MIN, ENERGY_MAX, RES)

    #retrieve the experimental spectrum from the charactersation
    intensity_expt, uncert_expt = use_characterisation(to_reconstruct, ENERGIES)

    #run the fitting routine and collect results
    fit_out = run_fit(thy, intensity_expt, uncert_expt, ENERGIES)
    opt_params = fit_out.params
    chi_red = fit_out.redchi
    residual = fit_out.residual
    intensity_thy_opt = thy.model(opt_params, ENERGIES)


    #visualise the fitting results, reconstruction and theoretical spectra
    ###USER_CHOICE###
    stick_colours = [] ###Colours to represent satellites in stick plot
    #MUST BE THE SAME LENGTH AS 'to_fit' (3 .ct files = 3 colours)
    ######

    plot_sticks(thy, to_fit, stick_colours, ENERGY_MIN, ENERGY_MAX,
                        stick_plot_title, stick_plot_path)

    plot_reconstruction(ENERGIES, intensity_expt, uncert_expt,
                          char_plot_title, char_plot_name)

    plot_optimal(ENERGIES, intensity_expt, intensity_thy_opt, residual, chi_red,
                        fit_plot_title, fit_plot_path)

    print(fit_report(fit_out, show_correl = False))


if __name__ == "__main__":
    main()
