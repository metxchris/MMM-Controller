# Standard Packages
from copy import deepcopy

# 3rd Party Packages
import numpy as np

# Local Packages
from main import *
from main.enums import ShotType, ScanType
from main.options import Options
from main.controls import InputControls
from plots import plot_profiles


def execute_basic_run(mmm_vars, controls):
    '''
    Executes a single MMM run, without varying any input parameters

    Creates an input file for the MMM driver using mmm_vars.  The MMM driver is then
    ran, which produces an output file.  This output file is parsed and a CSV of both
    the input and output data are stored, and an output profile PDF is created.

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * controls (InputControls): Specifies input control values in the MMM input file
    '''

    controls.save_controls(Options.instance)
    write_inputs.write_input_file(mmm_vars, controls)
    run_driver.run_mmm_driver()
    output_vars = read_output.read_output_file()
    plot_profiles.plot_output_profiles(output_vars)


def execute_variable_scan(mmm_vars, controls):
    '''
    Executes an input variable scan, where the values of an input variable are varied
    over a specified range and are then sent to the MMM driver for each value of the range

    Create a copy of mmm_vars as modified_vars to keep variables that are modified over the
    course of the scan separate from base MMM input variables.  For each factor of the scan_range,
    we modify the value of the specified var_to_scan, and then adjust any dependent variables.
    The MMM driver is ran each time var_to_scan is adjusted, and all input and output variable data
    is saved to a subfolder named after var_to_scan.  Afterwards, the saved CSV data is reshaped
    into data dependent on the scanned parameter, and is saved to another set of CSV within a
    new subfolder labeled rho.

    Parameter scan PDFs are not produced here, and the output data is intended to be plotted by
    a separate process after the scan is complete.

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * controls (InputControls): Specifies input control values in the MMM input file
    '''

    modified_vars = deepcopy(mmm_vars)
    var_to_scan = Options.instance.var_to_scan
    scan_range = Options.instance.scan_range

    # Create references to variable being scanned in mmm_vars and modified_vars
    # Modifying scanned_var values will modify its corresponding values in modified_vars
    base_var = getattr(mmm_vars, var_to_scan)
    scanned_var = getattr(modified_vars, var_to_scan)
    controls.save_controls(Options.instance)

    for i, scan_factor in enumerate(scan_range):
        print(f'Executing variable scan {i + 1} of {len(scan_range)} for variable {var_to_scan}')

        # Modifiy values of variable being scanned
        # Note: Dependent variables will be handled on a case-by-case basis
        scanned_var.values = scan_factor * base_var.values
        write_inputs.write_input_file(modified_vars, controls)
        run_driver.run_mmm_driver()
        read_output.read_output_file(scan_factor)

    # Reshaped scanned CSV into new CSV dependent on the scanned parameter
    parse_scans.parse_scan_csv()

    print('\nVariable scan complete!')


def execute_control_scan(mmm_vars, controls):
    '''
    Executes an input control scan, where the values of an input control are varied
    over a specified range and are then sent to the MMM driver for each value of the range

    Parameter scan PDFs are not produced here, and the output data is intended to be plotted by
    a separate process after the scan is complete.

    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    * controls (InputControls): Specifies input control values in the MMM input file
    '''

    var_to_scan = Options.instance.var_to_scan
    scan_range = Options.instance.scan_range

    # Create references to control being scanned in InputControls
    # Modifying scanned_control values will modify its corresponding values in controls
    scanned_control = getattr(controls, var_to_scan)
    base_control = deepcopy(scanned_control)

    for i, scan_factor in enumerate(scan_range):
        print(f'Executing control scan {i + 1} of {len(scan_range)} for control {var_to_scan}')

        # Modifiy values of variable being scanned
        # Note: Dependent variables will be handled on a case-by-case basis
        scanned_control.values = scan_factor * base_control.values
        controls.save_controls(Options.instance, scan_factor)
        write_inputs.write_input_file(mmm_vars, controls)
        run_driver.run_mmm_driver()
        read_output.read_output_file(scan_factor)

    # Reshaped scanned CSV into new CSV dependent on the scanned parameter
    parse_scans.parse_scan_csv()

    print('\nVariable scan complete!')


def initialize_variables():
    '''
    Initializes all input variables needed to run the MMM Driver and plot variable profiles

    Returns:
    * mmm_vars (InputVariables): All calculated variables, interpolated onto a grid of size input_points
    * input_vars (InputVariables): All calculated variables, interpolated onto XB+1 obtained from the CDF
    * cdf_vars (InputVariables): All CDF variables, interpolated onto XB+1 obtained from the CDF
    * raw_cdf_vars (InputVariables): All unedited CDF variables (saved for troubleshooting)
    '''

    raw_cdf_vars = read_cdf.read_cdf()
    cdf_vars = conversions.initial_conversion(raw_cdf_vars)
    input_vars = calculations.calculate_inputs(cdf_vars)
    mmm_vars = conversions.final_interpolation(input_vars)

    return mmm_vars, input_vars, cdf_vars, raw_cdf_vars


def main(controls):
    '''
    Main function which controls the MMM driver

    Needed output folders are created and a unique scan number is chosen for storing output data.
    All input variable objects are initialized and corresponding plot PDFs are created.  The MMM driver
    is then ran once, and then an optional variable scan can be ran afterwards.  Note that raw_cdf_vars
    does not exist on the same grid as other variable objects created here, and is only saved for
    debugging purposes.  Both input_vars and cdf_vars are guaranteed to be on the same grid, so
    these are used for profile comparisons.  mmm_vars is only on the same grid as input_vars if
    input_points is set to the same value as the size of XB+1 from the CDF, and if uniform_rho = False.

    Parameters:
    * controls (InputControls): Specifies input control values in the MMM input file
    '''

    print('Running MMM Controller...\n')

    utils.clear_temp_folder()
    utils.init_output_dirs(Options.instance)

    mmm_vars, input_vars, cdf_vars, raw_cdf_vars = initialize_variables()

    Options.instance.save_options()  # TODO: Create an event to save Options
    plot_profiles.plot_profile_comparison(cdf_vars, input_vars)
    plot_profiles.plot_input_profiles(mmm_vars)
    plot_profiles.plot_additional_profiles(mmm_vars)

    execute_basic_run(mmm_vars, controls)

    if Options.instance.scan_type == ScanType.VARIABLE:
        execute_variable_scan(mmm_vars, controls)
    elif Options.instance.scan_type == ScanType.CONTROL:
        execute_control_scan(mmm_vars, controls)


# Run this file directly to plot variable profiles and run the MMM driver
if __name__ == '__main__':
    '''
    CDF Options:
    * Uncomment the line you wish to use
    '''
    # cdf_name, shot_type, input_time = '129041A10', ShotType.NSTX, 0.5
    # cdf_name, shot_type, input_time = '120968A02', ShotType.NSTX, 0.5
    cdf_name, shot_type, input_time = '120982A09', ShotType.NSTX, 0.5
    # cdf_name, shot_type, input_time = '132017T01', ShotType.DIII_D, 2.1
    # cdf_name, shot_type, input_time = '141552A01', ShotType.DIII_D, 2.1

    '''
    Scan Options:
    * Uncomment the line you wish to use
    '''
    # var_to_scan, scan_range = None, None
    var_to_scan, scan_range = 'gti', np.arange(start=0.5, stop=3 + 1e-6, step=0.5)
    # var_to_scan, scan_range = 'gte', np.arange(start=0.025, stop=5 + 1e-6, step=0.025)
    # var_to_scan, scan_range = 'kyrhoe_etgm', np.arange(start=0, stop=5 + 1e-6, step=0.025)

    '''
    Input Options:
    * input_points is the number of points to use when making the MMM input file
    * Set input_points = None to match the number of points used in the CDF
    * Set uniform_rho = True to interpolate to a grid of evenly spaced rho values (takes longer)
    * apply_smoothing enables smoothing of all variables that have a smooth value set in the Variables class
    '''
    Options.instance.set(
        runid=cdf_name,
        shot_type=shot_type,
        input_time=input_time,
        input_points=201,
        uniform_rho=False,
        apply_smoothing=True,
        var_to_scan=var_to_scan,
        scan_range=scan_range,
    )

    '''
    Input Control Options:
    * cmodel controls enable the corresponding model if set to 1, and disable it if set to 0
    '''
    controls = InputControls(Options.instance)
    controls.set(
        cmodel_weiland=0,
        cmodel_dribm=0,
        cmodel_etg=0,
        cmodel_etgm=1,
        cmodel_mtm=0,
        etgm_kyrhoe=0.25,
        etgm_kyrhos=0.33,
        etgm_cl=1,
    )

    main(controls)
