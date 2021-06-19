"""
AUTHOR: Paarangat Pushkarna
DATE: 26/05/2021
MODIFIED: 11/06/2021
"""
### FILE PROCESSING MODULES ###
import re
import csv
import os

### DATA ANALYSIS MODULES ###
import pandas as pd
import numpy as np
from scipy import constants
from random import uniform

##FITTING MODULES##
from lmfit import Parameters
from itertools import combinations

class spectra:
    """
    ===================================
    Umbrella class to store information about multiple spectra.
    ===================================
    """
    def __init__(self, working_dir):
        """
        spectra class initialiser.
        Inputs
        ---------------------------------
        working_dir (str):
        path/directory where the ct files are placed.
        ---------------------------------
        NOTE: ct files must be named so that filenames do not contain dots
        or dashes (. or -), excepting the filename.
        """
        self.wd = working_dir #storage directory for ct files
        self.ct_names = self.retrieve_names()
        self.spectra = self.create_spectra() #dictionary of spectrum objects
        self.routine = self.define_routine()
        self.extrema = self.find_minmax()
        self.chosen = None

    def retrieve_names(self):
        """
        Retrieves ct filenames from working directory.
        Inputs
        -------------------------------
        self.wd: working directory where ct files are. str.
        Returns
        -------------------------------
        ct_names: list of ct filenames without extension. list of str.
        """
        ct_files = os.listdir(self.wd)
        ct_names = []
        for ct_file in ct_files:
            if ".ct" in ct_file: #only choose ct files, skip other file types
                name, ext = os.path.splitext(ct_file)
                ct_names.append(name) #only append name, and not extension
        return ct_names

    def create_spectra(self):
        """
        Creates a dictionary of spectrum objects from ct filenames.
        Inputs
        -------------------------------
        self.wd: workding directory where ct files are. str.
        self.ct_names: list of ct filenames in self.wd. list of str.
        Returns
        -------------------------------
        spectra: dictionary of spectrum objects.
            key: ct filename without extension
            value: 'spectrum' object
        """
        spectra = {}
        for ct_name in self.ct_names:
            spectra[ct_name] = spectrum(working_dir = self.wd, ct_file_name = ct_name)
        return spectra

    def define_routine(self):
        """
        Produces fitting routine for auto_fit using combinations.
        Inputs
        -------------------------------
        self.ct_names: list of ct filenames in self.wd. list of str.
        Returns
        -------------------------------
        all_combs: list of ct file combinations with diagram line. list of str.
        """
        all_combs = []
        names = self.ct_names
        for i in range(len(names)+1):
            combs = combinations(names, i)
            for comb in combs:
                comb = list(comb)
                dia_flag = 0
                for name in comb:
                    if "diagram" in name: #combination must include diagram line
                        dia_flag = 1
                if dia_flag:
                    all_combs.append(comb)
        return all_combs

    def create_parameters(self):
        """
        Initialise free parameters for fit.
        Inputs
        -------------------------------
        self.chosen: the ct files chosen to be included in fit. list of str.
        WARNING: follow naming convention for ct files outlined in
        /documentation/usage_READ_ME.txt
        Returns
        -------------------------------
        params: initialised free parameters for fit. lmfit 'Parameters' object.
        """
        params = Parameters()
        #'global' parameters, which are true for every fit
        params.add("fwhm", value = uniform(0,1), min = 0) #Lorentzian FWHM
        params.add("vert_scale", value = uniform(0,1), min = 0) #Overall vertical scaling
        params.add("E_offset", value = 0, vary = True) #Overall x offset
        params.add("alpha", value = 0, vary = False) #horizontal scaling parameter
        #fit specific parameters depending on the satellites chosen
        #Satellite weights/scaling factors
        for choice in self.chosen:
            params.add(choice, value = uniform(0,1), min = 0)
        return params

    def model(self, params, E):
        """
        Creates model from chosen ct files/chosen multiplets.
        Inputs
        -------------------------------
        params: Initialised free parameters from 'create_parameters'.
        lmfit 'Parameters' object
        E: array containing energies over specified energy range. np array.
        Returns:
        -------------------------------
        model: sum of spectra due to each satellite included. np array.
        """
        p = params.valuesdict()
        model = 0
        for choice in self.chosen: #loop over chosen ct files/multiplets
            _spectrum = self.spectra[choice]
            model += _spectrum.multiplet_thy(E, p["fwhm"], p["E_offset"],
                                            _spectrum.energies,
                                            _spectrum.line_strength,
                                            p["alpha"],
                                            p[choice] ) #compute multiplet
                                            #and add to model
        model = model*p["vert_scale"]
        return model

    def objective(self, params, E, data, uncert):
        """
        Objective function to be minimised. ((data - model)/uncertainty)
        Parameters
        ----------------------------------
        params: Initialised free parameters from 'create_parameters'.
        lmfit 'Parameters' object.
        E: array containing energies over specified energy range. np array.
        Returns:
        data: experimental data, reconstructed from characterisation. np array.
        uncert: total uncertainty associated with experimental data. np array.
        """
        model = self.model(params, E)
        return (data - model)/uncert

    def find_minmax(self):
        """
        Computes the global min/max energy value present in spectra.
        Inputs
        -------------------------------
        spectra: dictionary of spectrum objects earlier defined.
        Returns
        -------------------------------
        min: minimum transition energy for all spectrum obj in spectra. float.
        max: maximum transition energy for all spectrum obj in spectra. float.
        """
        min = 1e28
        max = 0
        for name in self.spectra.keys():
            energies = self.spectra[name].energies
            name_min = np.amin(energies)
            name_max = np.amax(energies)
            if name_max > max:
                max = name_max
            if name_min < min:
                min = name_min
        return min, max

class spectrum:
    """
    ===================================
    Stores information regarding a particular spectrum instance.
    ===================================
    """
    def __init__(self, ct_file_name, working_dir):
        """
        spectrum class initialiser.
        Inputs:
        -------------------------------
        ct_file_name: name of the ct file WITHOUT extension. str.
        working_dir: path where ct file is stored. str (path).
        """
        ### FILE I/O VARS ###
        self.wd = working_dir
        self.ct_fname = ct_file_name #w/o extension!
        self.ct_fpath = os.path.join(self.wd, self.ct_fname + '.ct')
        self.hdr_rm_name = '{}_hdr_rm.txt'.format(self.ct_fname)
        self.bab_name_txt = '{}_b.txt'.format(self.ct_fname)
        self.bab_name_csv = '{}_b.csv'.format(self.ct_fname)
        self.coul_name_txt = '{}_c.txt'.format(self.ct_fname)
        self.coul_name_csv = '{}_c.csv'.format(self.ct_fname)
        ### READING IN DATA TO A PD DATAFRAME ###
        self.clean_ct()
        self.csvb_path = os.path.join(self.wd, self.bab_name_csv)
        self.df_b = self.read_csv_ct(self.csvb_path)
        self.csvc_path = os.path.join(self.wd, self.coul_name_csv)
        self.df_c = self.read_csv_ct(self.csvc_path) ##CLASS VAR

        ### ATTRIBUTES OF EACH SPECTRUM GLEANED FROM DATAFRAME ###
        self.energies = self.df_c["E(eV)"].to_numpy() #transition energies
        self.tr_rate_b = self.df_b["A(s-1)"].to_numpy() #transition rates
        self.tr_rate_c = self.df_c["A(s-1)"].to_numpy()
        self.oscillator_strength = self.df_c["gf"].to_numpy()
        self.line_strength = self.df_c["S"].to_numpy()

    ## FILE PROCESSING METHODS ##
    def clean_ct(self):
        """
        Master function which controls the processing of the ct file
        Inputs
        -------------------------------
        misc. filenames: all filenames for data processing in __init__. str.
        self.wd: working directory where ct files are. str.
        """
        self.remove_header()
        self.split_gauge()
        self.ct_csv(self.bab_name_txt, self.bab_name_csv)
        self.ct_csv(self.coul_name_txt, self.coul_name_csv)

        #Delete intermediate files. Can turn off this code if intermediate files
        #are to be retained.
        os.remove(os.path.join(self.wd, self.hdr_rm_name))
        os.remove(os.path.join(self.wd, self.bab_name_txt))
        os.remove(os.path.join(self.wd, self.coul_name_txt))
        return None

    def ct_csv(self, inp_fname, out_fname):
        """
        Write the csv file.
        Input
        -------------------------------
        inp_fname: name of text file to be read. str.
        out_fname: name of csv file to be produced. str.
        """
        inp_fpath = os.path.join(self.wd, inp_fname)
        out_fpath = os.path.join(self.wd, out_fname)
        with open(out_fpath, mode = 'w') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter = ',')
            with open(inp_fpath, 'r') as ct_text:
                for line in ct_text.readlines():
                    line = line.split()
                    csv_writer.writerow(line)
        return None

    def remove_header(self):
        """
        Removes header content of ct file.
        Inputs
        -------------------------------
        self.wd: working directory where ct files are. str.
        self.hdr_rm_name: name of ct with header removed.
        self.ct_fpath: full path of ct file.
        """
        out_fpath = os.path.join(self.wd, self.hdr_rm_name)
        inp_fpath = self.ct_fpath

        with open(out_fpath, 'w') as inter_file:
            with open(inp_fpath, 'r') as ct_file:
                line_generator = (line for line in ct_file.readlines())
                for line in line_generator:
                    if 'File' in line:
                        inter_file.write(self.edit_header(line))
                        #also fix the column naming
                        break #once you find the column names for the data table
                        #record that to the file and break out of the generator
                for line in line_generator:
                    #Accessing generator again starts the line read sequence from
                    #where we left off
                    inter_file.write(line)
        return None

    def split_gauge(self):
        """
        Splits file into two separate text files by gauge difference.
        Inputs
        -------------------------------
        self.wd: working directory where ct files are. str.
        self.bab_name_txt: name of the output file in Babushkin gauge. str.
        self.coul_name_txt: name of the output file in Coulomb gauge. str.
        self.hdr_rm_name: name of file with header removed. str.
        """
        out_fpath_b = os.path.join(self.wd, self.bab_name_txt)
        out_fpath_c = os.path.join(self.wd, self.coul_name_txt)
        inp_fpath = os.path.join(self.wd, self.hdr_rm_name)

        babushkin_file = open(out_fpath_b, 'w')
        coulomb_file = open(out_fpath_c, 'w')

        pivot = re.compile(r'([0-9]|\s)\s+[BC]') #pinpoints pattern that marks the
        #split between coulomb and babushkin gauges in the ct file.

        with open(inp_fpath, 'r') as ct_file_noh:
            line_reader = ct_file_noh.readlines()
            #A list is created of every line in the ct file because we want to
            #'peek' ahead to see what the values of the amplitudes are in babushkin
            #gauge, so reference each line by index 'i'
            for i in range(len(line_reader)):
                line = line_reader[i]
                if not i: #If we are dealing with column names i = 0
                    babushkin_file.write(line)
                    coulomb_file.write(line)
                matches = pivot.finditer(line) #setup the regex pattern matcher
                for match in matches:
                    start = match.span()[0]
                    end = match.span()[1]
                    common = line[:end-1] #this is the common energy between both gauges
                    match_str = match.group()
                    if 'C' in match_str:
                        #if we have coulomb gauge, peek ahead to see what the
                        #babushkin gauge values are
                        col = line_reader[i][end:]
                        bab = line_reader[i+1][end:]
                        col = self.change_sci_not(col) #change scientific notation
                        bab = self.change_sci_not(bab)
                        coulomb_file.write(common + col)
                        babushkin_file.write(common + bab)

        babushkin_file.close()
        coulomb_file.close()

        return None

    def edit_header(self, line):
        """
        Edits header to comply with pd Dataframe standards.
        Inputs
        -------------------------------
        line: str being considered in loop over lines of text file. str.
        Returns
        -------------------------------
        line: modified header line.
        """
        #identify replacements to be made
        diff_file = re.compile(r'(File)\s+(Lev)\s+(J)\s+(P)\s+(File)\s+(Lev)\s+(J)\s+(P)')
        kays = re.compile(r'E\s\(Kays\)')
        einstein = re.compile(r'A\s\(')
        #make replacements
        line = diff_file.sub(r'File1 Lev1 J1 P1 File2 Lev2 J2 P2', line)
        line = kays.sub(r'E(Kays)', line)
        line = einstein.sub(r'A(', line)
        #handle eV case
        if 'eV' in line:
            ev = re.compile(r'E\s\(\seV\s\)')
            line = ev.sub(r'E(eV)', line)
        return line

    def change_sci_not(self, line):
        """
        Change old FORTRAN scientific notation, to modern 'e' notation.
        Inputs
        -------------------------------
        line: line to be modified in text file. str.
        Returns
        -------------------------------
        new_line: modified line.
        """
        pattern = re.compile(r'D')
        new_line = pattern.sub(r'e', line)
        return new_line

    ##DATA READING METHODS##
    def read_csv_ct(self, inp_csv_path):
        """
        Reads csv file into a pandas Dataframe.
        Inputs
        -------------------------------
        inp_csv_path: path where candidate csv file is stored.
        Returns
        -------------------------------
        df: pandas Dataframe with spectrum information
        """
        df = pd.read_csv(inp_csv_path)
        header_list = df.columns.to_numpy()
        einstein_coeff = df['A(s-1)'].to_numpy()
        if 'E(Kays)' in header_list: #make unit conversion
            wavenums = df['E(Kays)'].to_numpy()
            energy = self.kays_to_eV(wavenums)
            df.rename(columns = {'E(Kays)': 'E(eV)'}, inplace = True)
            df['E(eV)'] = energy
        return df

    def kays_to_eV(self, wavenum_values):
        """
        Converts units of wavelength (Kays or cm^-1) into units of energy (eV)
        """
        try:
            h_ev = constants.physical_constants["Planck constant in eV s"][0]
        except:
            h_ev = constants.physical_constants["Planck constant in eV/Hz"][0]
        c = constants.c
        wavelength_values_cm = (1.)/wavenum_values
        wavelength_values_m = wavelength_values_cm/100
        energy_values = (h_ev*c)/wavelength_values_m
        return energy_values

    ##BUILD MULTIPLET USING SPECTRUM VALUES##
    def multiplet_thy(self, E, fwhm, E_offset, E_stick, A_stick, alpha, A_mult):
        """
        Broadens sharp transition energies with lorentzian broadening.
        Sums all broadened lorentzians to create the 'multiplet'.
        Inputs
        -------------------------------
        E: energy values for defined energy range. np array.
        fwhm: Free lorentzian FWHM parameter. float.
        E_offset: Free lorentzian horizontal offset parameter. float.
        E_stick: Transition energies for the spectrum. np array.
        A_stick: Line strengths for each transition energy. np array.
        alpha: Horizontal scaling free parameter. float.
        A_mult: Vertical scaling free parameter. float.
        Returns
        -------------------------------
        sum*A_mult: Sum of broadened lorentzians, scaled by A_mult. np array.
        """
        #uses np matrix manipulation to omit for loop over transition energies.
        #matrix has lorentzian for each transition energy as row.
        E_matrix = E - E_stick.reshape(E_stick.size, 1)
        energy_term = E_matrix*(1+alpha) - E_offset
        lambda_term = 0.5*fwhm

        numerator = np.multiply(lambda_term, A_stick)
        denominator = np.square(energy_term) + lambda_term**2
        unscaled = np.divide(1.0, denominator)
        lorentzians = np.multiply(unscaled, numerator[:, np.newaxis])
        sum = np.sum(lorentzians, axis = 0)

        return sum*A_mult

class characterisation:
    """
    ===================================
    Retrieve the characterisation of experimental data.
    ===================================
    """
    def __init__(self, csv_path):
        """
        charactersation class initialiser.
        Inputs
        -------------------------------
        csv_path: path to csv file with characterisation.
        """
        self.df = pd.read_csv(csv_path) #df defined by reading csv
        self.centroids = self.df["C_i"].to_numpy()
        self.centroid_us = self.df["u C_i"].to_numpy() #uncertainty in centroid
        self.fwhms = self.df["W_i"].to_numpy()
        self.fwhm_us = self.df["u W_i"].to_numpy() #uncertainty in FWHM
        self.areas = self.df["A_i"].to_numpy()
        self.area_us = self.df["u A_i"].to_numpy() #uncertainty in int. Area

    def reconstruct(self, E):
        """
        Uses lorentzian broadening to 'reconstruct' the deconvolved
        experimental spectrum.
        Inputs
        -------------------------------
        E: energy values over defined energy range. np array.
        Returns
        -------------------------------
        sum: sum of lorentzians broadened by specified FWHM in characterisation
        """
        OFFSET = 0
        ALPHA = 0
        A_MULT = 1
        sum = self.multiplet_expt(E, self.fwhms, OFFSET, self.centroids,
                                    self.areas, ALPHA, A_MULT)
        return sum

    def uncertainties(self, E):
        """
        Compute total uncertainty associated with intensity values in
        reconstruction.
        Inputs
        -------------------------------
        E: energy values over defined energy range. np array.
        Returns
        -------------------------------
        np.sqrt(sum_sq): Error in each intensity value from characterisation.
        """
        sum_sq = 0
        for i in range(len(self.centroids)):
            area_err = (self.delRdelA(E, self.centroids[i], self.fwhms[i])**2)*(self.area_us[i]**2)
            fwhm_err = (self.delRdelW(E, self.centroids[i], self.fwhms[i], self.areas[i])**2)*(self.fwhm_us[i]**2)
            centroid_err = (self.delRdelW(E, self.centroids[i], self.fwhms[i], self.areas[i])**2)*(self.centroid_us[i]**2)
            sum_sq += area_err + fwhm_err + centroid_err
        return np.sqrt(sum_sq)

    def delRdelA(self, E, C, W):
        """
        Compute uncertainty in integrated area.
        """
        numerator = (W/2)**2
        denominator = (E - C)**2 + numerator
        scale = 1/np.pi
        return np.divide(numerator, denominator)*scale

    def delRdelW(self, E, C, W, A):
        """
        Compute uncertainty in FWHM.
        """
        scale = A/np.pi
        numerator = W*((E - C)**2)
        denominator = ((E - C)**2 + (W/2)**2)**2
        return np.divide(numerator, denominator)*scale

    def delRdelC(self, E, C, W, A):
        """
        Compute uncertainty in centroid.
        """
        scale = A/(2*np.pi)
        numerator = (W**2)*(E - C)
        denominator = ((E - C)**2 + (W/2)**2)**2
        return np.divide(numerator, denominator)*scale

    def multiplet_expt(self, E, fwhm, E_offset, E_stick, A_stick, alpha, A_mult):
        """
        Compute sum over broadened centroid values. Identical to 'multiplet_thy'
        but treats 'fwhm' as an np array instead of float.
        """
        E_matrix = E - E_stick.reshape(E_stick.size, 1)
        energy_term = E_matrix*(1+alpha) - E_offset
        lambda_term = 0.5*fwhm

        numerator = np.multiply(lambda_term, A_stick)
        denominator = np.square(energy_term) + np.square(lambda_term[:, np.newaxis])
        unscaled = np.divide(1.0, denominator)
        lorentzians = np.multiply(unscaled, numerator[:, np.newaxis])
        sum = np.sum(lorentzians, axis = 0)

        return sum*A_mult
