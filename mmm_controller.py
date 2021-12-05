# Standard Packages
from copy import deepcopy

# 3rd Party Packages
import numpy as np

# Local Packages
from main import *
from main.enums import ShotType
from main.options import Options
from plots import plot_profiles
from tests import test


def execute_basic_run(mmm_vars):
    '''
    Executes a single MMM run, without varying any input parameters
  
    Creates an input file for the MMM driver using mmm_vars.  The MMM driver is then
    ran, which produces an output file.  This output file is parsed and a CSV of both
    the input and output data are stored, then an output profile PDF is created.
  
    Parameters:
    * mmm_vars (InputVariables): Contains all variables needed to write MMM input file
    '''

    write_inputs.write_input_file(mmm_vars)
    run_driver.run_mmm_driver()
    output_vars = read_output.read_output_file()
    plot_profiles.plot_output_profiles(output_vars)

def execute_variable_scan(mmm_vars):
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
    '''

    modified_vars = deepcopy(mmm_vars)
    input_options = Options.instance
    var_to_scan = input_options.var_to_scan
    scan_range = input_options.scan_range

    # Create references to variable being scanned in mmm_vars and modified_vars
    # Modifying scanned_var values will modify its corresponding values in modified_vars
    base_var = getattr(mmm_vars, var_to_scan)
    scanned_var = getattr(modified_vars, var_to_scan)

    for i, scan_factor in enumerate(scan_range):
        print(f'Executing variable scan {i + 1} of {len(scan_range)} for variable {var_to_scan}')

        # Modifiy values of variable being scanned, and store the scan_factor
        # Note: Dependent variables will be handled on a case-by-case basis
        scanned_var.set_variable(scan_factor * base_var.values)
        input_options.scan_factor_str = scan_factor

        write_inputs.write_input_file(modified_vars)
        run_driver.run_mmm_driver()
        read_output.read_output_file()

    # Reshaped scanned CSV into new CSV dependent on the scanned parameter
    parse_scans.parse_scan_csv()

    print('\nVariable scan complete!')

def initialize_controller():
    '''
    Initializes all input variables needed to run the MMM Driver and plot variable profiles
  
    The temp folder is used to store individual PDF sheets before they are merged, and is 
    cleaned out at the start of each initialization.  The CDF is then read for variables needed
    to create the MMM input file, and then various conversions and calculations are made to
    create the mmm_vars object needed to write the MMM Driver input file.

    Returns:
    * mmm_vars (InputVariables): All calculated variables, interpolated onto a grid of size input_points 
    * input_vars (InputVariables): All calculated variables, interpolated onto XB+1 obtained from the CDF
    * cdf_vars (InputVariables): All CDF variables, interpolated onto XB+1 obtained from the CDF
    * raw_cdf_vars (InputVariables): All unedited CDF variables (saved for troubleshooting)
    '''

    utils.clear_temp_folder()
    raw_cdf_vars = read_cdf.read_cdf()
    cdf_vars = convert_inputs.initial_conversion(raw_cdf_vars)
    input_vars = calculate_inputs.calculate_inputs(cdf_vars)
    mmm_vars = convert_inputs.final_interpolation(input_vars)

    return mmm_vars, input_vars, cdf_vars, raw_cdf_vars

def run_controller():
    '''
    Main function which controls the MMM driver

    All input variable objects are initialized and corresponding plot PDFs are created.  The MMM driver
    is then ran once, and then an optional variable scan can be ran afterwards.  Note that raw_cdf_vars
    does not exist on the same grid as other variable objects created here, and is only saved for
    debugging purposes.  Both input_vars and cdf_vars are guaranteed to be on the same grid, so
    these are used for profile comparisons.  mmm_vars is only on the same grid as input_vars if 
    input_points is set to the same value as the size of XB+1 from the CDF, and if uniform_rho = False.
    '''

    mmm_vars, input_vars, cdf_vars, raw_cdf_vars = initialize_controller()
    plot_profiles.plot_profile_comparison(cdf_vars, input_vars)
    plot_profiles.plot_input_profiles(mmm_vars)
    plot_profiles.plot_additional_profiles(mmm_vars)
    
    execute_basic_run(mmm_vars)

    if Options.instance.var_to_scan is not None:
        execute_variable_scan(mmm_vars)


# Run this file directly to plot variable profiles and run the MMM driver
if __name__ == '__main__':
    '''
    CDF Options: 
    * Uncomment the line you wish to use
    * Edit enums.py to view or add additional ShotTypes
    '''
    cdf_name, shot_type, input_time = '129041A10', ShotType.NSTX, 0.5
    # cdf_name, shot_type, input_time = '120982A09', ShotType.NSTX, 0.5
    # cdf_name, shot_type, input_time = '132017T01', ShotType.D3D, 2.1
    # cdf_name, shot_type, input_time = '141552A01', ShotType.D3D, 2.1

    '''
    Input Options:
    * input_points is the number of points to use when making the MMM input file
    * Set uniform_rho = True to interpolate to a grid of evenly spaced rho values (takes longer)
    * Set input_points = None to match the number of points used in the CDF
    * Set var_to_scan = var_name (str) to run a scan of the specified variable
    * Set var_to_scan = None to skip the variable scan
    * E.g.: var_to_scan = 'te'
    '''
    Options.instance.set_options(
        runid = cdf_name,
        shot_type = shot_type,
        input_time = input_time,
        input_points = 51,
        uniform_rho = False,
        apply_smoothing = False,
        var_to_scan = 'gte',
        scan_range = np.arange(start=0.1, stop=3.1, step=0.1))

    run_controller()
