#!/usr/bin/python3
# Standard Packages
from copy import deepcopy

# 3rd Party Packages
import numpy as np

# Local Packages
import settings
import modules.options as options
import modules.controls as controls
import modules.utils as utils
import modules.adjustments as adjustments
import modules.parse_scans as parse_scans
import modules.mmm as mmm
import plotting.modules.profiles as profiles
from modules.enums import ShotType, ScanType, ProfileType


def _execute_basic_run(mmm_vars, input_controls):
    '''
    Executes a single MMM run, without varying any input parameters

    Creates an input file for the MMM driver using mmm_vars.  The MMM driver
    is then ran, which produces an output file.  This output file is parsed
    and a CSV of both the input and output data are stored, and an output
    profile PDF is created.

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * input_controls (InputControls): Specifies input control values in the MMM input file
    '''

    output_vars = mmm.run_wrapper(mmm_vars, input_controls)
    output_vars.save_all_vars(options.instance)
    profiles.plot_profiles(ProfileType.OUTPUT, output_vars)


def _execute_variable_scan(mmm_vars, input_controls):
    '''
    Executes an input variable scan, where the values of an input variable are
    varied over a specified range and are then sent to the MMM driver for
    each value of the range

    Create a copy of mmm_vars as modified_vars to keep variables that are
    modified over the course of the scan separate from base MMM input
    variables.  For each factor of the scan_range, we modify the value of the
    specified var_to_scan, and then adjust any dependent variables. The MMM
    driver is ran each time var_to_scan is adjusted, and all input and output
    variable data is saved to a subfolder named after var_to_scan.
    Afterwards, the saved CSV data is reshaped into data dependent on the
    scanned parameter, and is saved to another set of CSV within a new
    subfolder labeled rho.

    Parameter scan PDFs are not produced here, and the output data is intended
    to be plotted by a separate process after the scan is complete.

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write the MMM input file
    * input_controls (InputControls): Specifies input control values in the MMM input file
    '''

    var_to_scan = options.instance.var_to_scan
    scan_range = options.instance.scan_range

    for i, scan_factor in enumerate(scan_range):
        print(f'Executing variable scan {i + 1} of {len(scan_range)} for variable {var_to_scan}')
        adjusted_vars = adjustments.adjust_scanned_variable(mmm_vars, var_to_scan, scan_factor)
        adjusted_vars.save_all_vars(options.instance, scan_factor)
        output_vars = mmm.run_wrapper(adjusted_vars, input_controls)
        output_vars.save_all_vars(options.instance, scan_factor)

    parse_scans.create_rho_files()  # Creates CSVs in the rho folder

    print('\nVariable scan complete!')


def _execute_control_scan(mmm_vars, input_controls):
    '''
    Executes an input control scan, where the values of an input control are
    varied over a specified range and are then sent to the MMM driver for
    each value of the range

    Parameter scan PDFs are not produced here, and the output data is intended
    to be plotted by a separate process after the scan is complete.

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * input_controls (InputControls): Specifies input control values in the MMM input file
    '''

    var_to_scan = options.instance.var_to_scan
    scan_range = options.instance.scan_range

    # Create a reference to control being scanned in InputControls. Modifying
    # scanned_control values will modify its corresponding values in
    # adjusted_controls

    adjusted_controls = deepcopy(input_controls)
    scanned_control = getattr(adjusted_controls, var_to_scan)
    base_control = getattr(input_controls, var_to_scan)

    for i, scan_factor in enumerate(scan_range):
        print(f'Executing control scan {i + 1} of {len(scan_range)} for control {var_to_scan}')
        scanned_control.values = scan_factor * base_control.values
        mmm_vars.save_all_vars(options.instance, scan_factor)
        adjusted_controls.save_to_csv(options.instance, scan_factor)
        output_vars = mmm.run_wrapper(mmm_vars, adjusted_controls)
        output_vars.save_all_vars(options.instance, scan_factor)

    parse_scans.create_rho_files()  # Creates CSVs in the rho folder

    print('\nInput control scan complete!')


def main(scanned_vars, input_controls):
    '''
    Runs the MMM controller

    Needed output folders are created and a unique scan number is chosen for
    storing output data. All input variable objects are initialized and
    corresponding plot PDFs are created.  The MMM driver is then ran once,
    and then an optional variable scan can be ran afterwards.

    Parameters:
    * scanned_vars (Dict): Dictionary of variables being scanned
        - keys (str or None): The variable being scanned
        - values (np.ndarray or None): The range of factors to scan over
    * input_controls (InputControls): Specifies input control values in the MMM input file
    '''

    # TODO: Add validation for all items in scanned_vars

    for var_to_scan, scan_range in scanned_vars.items():
        options.instance.set(var_to_scan=var_to_scan, scan_range=scan_range)

        print(f'Running MMM Controller...\n')

        utils.clear_temp_folder()
        utils.init_output_dirs(options.instance)

        mmm_vars, cdf_vars, __ = utils.initialize_variables()

        options.instance.save()  # Need to be saved after variable initialization
        input_controls.update_from_options(options.instance)
        input_controls.save_to_csv(options.instance)
        mmm_vars.save_all_vars(options.instance)

        profiles.plot_profiles(ProfileType.INPUT, mmm_vars)
        profiles.plot_profiles(ProfileType.ADDITIONAL, mmm_vars)
        profiles.plot_profiles(ProfileType.COMPARED, mmm_vars, cdf_vars)

        # Basic runs create output profile plots and save output profile CSVs
        _execute_basic_run(mmm_vars, input_controls)

        if options.instance.scan_type == ScanType.VARIABLE:
            _execute_variable_scan(mmm_vars, input_controls)
        elif options.instance.scan_type == ScanType.CONTROL:
            _execute_control_scan(mmm_vars, input_controls)


# Run this file directly to plot variable profiles and run the MMM driver
if __name__ == '__main__':
    scanned_vars = {}

    '''
    TRANSP Data:
    * Uncomment the line you wish to use
    '''
    # runid, shot_type, input_time = '120968A02', ShotType.NSTX, 0.5
    # runid, shot_type, input_time = '120982A09', ShotType.NSTX, 0.5
    # runid, shot_type, input_time = '129041A10', ShotType.NSTX, 0.5
    # runid, shot_type, input_time = '138536A01', ShotType.NSTX, 0.630
    # runid, shot_type, input_time = '132017T01', ShotType.DIII_D, 2.1
    # runid, shot_type, input_time = '141552A01', ShotType.DIII_D, 2.1
    runid, shot_type, input_time = 'TEST', ShotType.NSTX, 0.5

    '''
    Scanned Variables:
    * Uncomment the lines you wish to include in scanned_vars
    * Using None as the scanned variable will skip the variable scan
    '''
    scanned_vars[None] = None
    # scanned_vars['gni'] = np.arange(start=0.5, stop=5 + 1e-6, step=0.5)

    # scanned_vars['betae'] = np.arange(start=0.05, stop=6 + 1e-6, step=0.05)
    # scanned_vars['btor'] = np.arange(start=0.025, stop=3 + 1e-6, step=0.025)
    # scanned_vars['etae'] = np.arange(start=0.025, stop=3 + 1e-6, step=0.025)
    # scanned_vars['etgm_kyrhoe'] = np.arange(start=0.05, stop=6 + 1e-6, step=0.05)
    # scanned_vars['etgm_kyrhos'] = np.arange(start=0.05, stop=6 + 1e-6, step=0.05)
    # scanned_vars['gnh'] = np.arange(start=0.025, stop=3 + 1e-6, step=0.025)
    # scanned_vars['gnz'] = np.arange(start=0.05, stop=9 + 1e-6, step=0.05)
    # scanned_vars['gte'] = np.arange(start=0.025, stop=6 + 1e-6, step=0.05)
    # scanned_vars['nuei'] = np.arange(start=0.025, stop=3 + 1e-6, step=0.025)
    # scanned_vars['q'] = np.arange(start=0.6, stop=2.4 + 1e-6, step=0.015)
    # scanned_vars['shear'] = np.arange(start=-6.0, stop=6 + 1e-6, step=0.1)
    # scanned_vars['tau'] = np.arange(start=0.025, stop=3 + 1e-6, step=0.025)
    # scanned_vars['zeff'] = np.arange(start=0.02, stop=4 + 1e-6, step=0.02)**2

    '''
    Options:
    * input_points is the number of points to use when making the MMM input file
    * Set input_points = None to match the number of points used in the CDF
    * Set uniform_rho = True to interpolate to a grid of evenly spaced rho values (takes longer)
    * apply_smoothing enables smoothing of all variables that have a smooth value set in the Variables class
    '''
    options.instance.set(
        runid=runid,
        shot_type=shot_type,
        input_time=input_time,
        input_points=101,
        uniform_rho=True,
        apply_smoothing=True,
    )

    '''
    Input Controls:
    * cmodel enables (disables) the corresponding model if set to 1 (0)
    '''
    input_controls = controls.InputControls()
    input_controls.set(
        cmodel_weiland=0,
        cmodel_dribm=0,
        cmodel_etg=0,
        cmodel_etgm=1,
        cmodel_mtm=0,
        etgm_kyrhoe=0.25,
        etgm_kyrhos=0.33,
        etgm_cl=1,  # etgm_cl=0 is collisionless, etgm_cl=1 is collisional
        etgm_exbs=1,
    )

    settings.AUTO_OPEN_PDFS = 0

    main(scanned_vars, input_controls)
