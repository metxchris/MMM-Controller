# Standard Packages
from copy import deepcopy
# 3rd Party Packages
import numpy as np
# Local Packages
from main import *
from plots import plot_profiles, plot2d

# Run the MMM Driver once, and show output profile plots
def execute_basic_run(input_vars, input_options):
    # Write variables to input file for MMM Driver
    write_inputs.write_input_file(input_vars, input_options)

    # Run MMM driver to produce output file
    run_driver.run_mmm_driver(input_options)

    # Read output variables from output file and save values to a CSV
    output_vars = read_output.read_output_file(input_options)

    # Plot output profiles
    plot_profiles.plot_output_profiles(output_vars, input_options)

# Run the MMM Driver multiple times while varying values of an input variable by a multiplicative factor
def execute_variable_scan(input_vars, input_options):
    # Duplicate input_vars to modified_vars for the variable scan
    modified_vars = deepcopy(input_vars)

    # Create references to variable being scanned in input_vars and modified_vars
    # Modifying scanned_var values will modify its corresponding values in modified_vars
    base_var = getattr(input_vars, input_options.var_to_scan)
    scanned_var = getattr(modified_vars, input_options.var_to_scan)

    for i, scan_factor in enumerate(input_options.scan_range):
        print('Executing variable scan {0} of {1} for variable {2}'
            .format(i + 1, len(input_options.scan_range), input_options.var_to_scan))

        # Modifiy values of variable being scanned, and store the scan_factor.
        # Note: We are intentionally not recalculating dependent variables in this step
        scanned_var.set_variable(scan_factor * base_var.values)
        input_options.scan_factor_str = scan_factor

        # Write modified variables to input file for MMM Driver
        write_inputs.write_input_file(modified_vars, input_options)

        # Run MMM driver to produce output file
        run_driver.run_mmm_driver(input_options)

        # Read output variables from output file and save values to a CSV
        read_output.read_output_file(input_options)
    
    # TODO: Plot results of variable scan

    print('Variable scan complete!')

def execute_profile_comparison(input_options):
    # Clear temp folder
    utils.clear_temp_folder()

    # Read variables from specified CDF
    cdf_vars = read_cdf.read_cdf(input_options)

    # Initial conversion of variables from CDF format to MMM format
    cdf_vars = convert_inputs.initial_conversion(cdf_vars, input_options)
    input_vars = deepcopy(cdf_vars)

    # Calculate new variables from CDF variables
    calculate_inputs.calculate_inputs(input_vars)

    plot_profiles.plot_profile_comparison(cdf_vars, input_vars, input_options)

    # cdf_vars.nh.values *= 10**5
    # plot2d.plot(input_options, cdf_vars.xb, cdf_vars.nh, r' (NH $\times 10^5$)', input_vars.xb, input_vars.nh)

# Initializes all input variables needed to run the MMM Driver
def initialize_controller(input_options):
    # Clear temp folder
    utils.clear_temp_folder()

    # Read variables from specified CDF
    cdf_vars = read_cdf.read_cdf(input_options)

    # Initial conversion of variables from CDF format to MMM format
    input_vars = convert_inputs.initial_conversion(cdf_vars, input_options)

    # Calculate new variables from CDF variables
    calculate_inputs.calculate_inputs(input_vars)

    # Final conversion: Interpolate onto larger grid of points
    convert_inputs.final_conversion(input_vars, input_options)

    # Plot input profiles being sent to the MMM driver and save as PDF
    plot_profiles.plot_input_profiles(input_vars, input_options)

    # TODO: add option to plot all variables calculatd from input variables

    # Execute a basic run to plot output profiles and base output values csv
    execute_basic_run(input_vars, input_options)

    # Condition for executing a variable scan
    if input_options.var_to_scan is not None:
        execute_variable_scan(input_vars, input_options)
    
if __name__ == '__main__':
    # CDF Options
    cdf_name = '132017T01'
    shot_type = 'DIII-D'
    input_time = 2.1

    # Interpolation points used for creating mmm input file
    input_points = 51

    """
    Variable Scan Options
    Set var_to_scan = var_name (str) to run a scan of the specified variable
    Set var_to_scan = None to skip the variable scan
    E.g.: var_to_scan = 'te'
    """
    var_to_scan = 'te'
    scan_range = np.arange(start=0.5, stop=2.1, step=0.1)

    # Save Input Options
    input_options = variables.InputOptions(cdf_name, shot_type, input_time, input_points)
    input_options.set_scan_values(var_to_scan, scan_range)

    # Run MMM Controller
    initialize_controller(input_options)

    # (Optional) Compare calculated profiles with those found in the CDF
    execute_profile_comparison(input_options)
