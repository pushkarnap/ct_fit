AUTHOR: Paarangat Pushkarna
DATE: 29/05/2021

.ct FILE NAMING CONVENTIONS
--------------------------------------------------------------------------------
FORMAT: '[element]_[diagram line]_[satellite]_[calculation level].ct'
To execute the automated fitting routine, .ct files must be named with the
above format. This is because the fitting package used does not accept
dots or dashes for parameter names. It is recommended that this
naming convention be used throughout.

USAGE INSTRUCTIONS
--------------------------------------------------------------------------------
AUTO-FITTING with 'sc_kb_fit_auto.py'
---------------------------------------
To initiate full fitting process, open 'sc_kb_fit_auto.py' and on line 14,
change 'ct_folder_name' to one of:
--- 'sc_kb_ct_4s' (.ct files for calculations up to the 4s level)
--- 'sc_kb_ct_4f' (.ct files for calculations up to the 4f level)
--- 'sc_kb_ct_5s' (.ct files for calculations up to the 5s level)
and run 'sc_kb_fit_auto.py' in the terminal.
This will fit all possible combinations of the files within the chosen
directory that involve the diagram line, to the experimental spectra.
WARNING: This process may take ~24Hrs, i.e. should be left to run overnight.
It will provide coarse results, and is intended to give a cursory indication
of the satellites which contribute.

'sc_kb_fit_auto.py' will produce 'chis.txt' in the 'ct_pythoncode' directory.
This is a record of reduced chi squared values for each combination of ct
files. It is recommended that 'chis.txt' is moved to another directory.

MANUAL FITTING with 'sc_kb_fit_manual.py'
---------------------------------------
To initiate the fitting of one particular combination of ct files,
use 'sc_kb_fit_manual.py'. This gives the user finer control over plotting and
reporting of fitting results.
Sections marked '###USER_CHOICE###' indicate where in 'sc_kb_fit_manual.py'
user input is required. Once all these sections are complete, run
'sc_kb_fit_manual.py' in the terminal.
If fitting report is to be saved (chi-squared, function evals, etc.),
redirect terminal output to a text file.
NOTE: Manual fitting is intended for use after the coarse fitting is complete,
although can be run prior without issues. Once interesting fits in the
automated run are pinpointed, particular fits can be re-run manually and
viewed in more detail.

'sc_kb_fit_manual.py' will produce a reconstructed plot in the 'ct_pythoncode'
directory, a stick plot in the ct file directory, and a fitted plot in the
ct file directory. The names of the latter plots are user defined.

CHANGING FITTING PARAMETERS in 'ct_process.py'
---------------------------------------
To change values of initial guesses, which parameters to fix and vary,
and the bounds on each parameter, navigate to lines 103-113 in 'ct_process.py',
under the function 'create_parameters'. Here, parameters for the fit
will be added using the 'params.add()' function.

To change the initial guess, change the 'value' input
within the 'params.add()' function call.
To change the bounds, change the 'min'/'max' input within the 'params.add()'
function call.
To fix/vary a parameter, set 'vary' to False/True within the 'params.add()'
function call.

WARNING: Change the minimum bounds on scaling parameters at your own risk.
Removing bounds may cause fits to have negative scaling parameters, which is
unphysical. 

For more information on the lmfit 'Parameters' object visit:
https://lmfit.github.io/lmfit-py/parameters.html
For more information on the lmfit package, visit:
https://lmfit.github.io/lmfit-py/intro.html

MISCELLANEOUS
--------------------------------------------------------------------------------
Please ignore files in the 'superseded' directory. These include
obsolete/deprecated code.
Different notation for Planck's constant in eV is used in different versions
of scipy.
If scipy.constants.physical_constants raises a 'KeyError' exception
for the value of 'Planck constant in eV s', try replacing this with
'Planck constant in eV/Hz'.
