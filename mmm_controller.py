# Standard Packages
from copy import deepcopy
# 3rd Party Packages
import numpy as np
# Local Packages
from main import *
from main.enums import ShotType
from plots import plot_profiles

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

# Initializes all input variables needed to run the MMM Driver
def initialize_controller(input_options):
    # Clear temp folder
    utils.clear_temp_folder()

    # Read variables from specified CDF, exactly as they are in the CDF
    raw_cdf_vars = read_cdf.read_cdf(input_options)

    # Initial conversion of variables (converts units, interpolates to XBO, and applies smoothing)
    cdf_vars = convert_inputs.initial_conversion(raw_cdf_vars, input_options)

    # Calculate new variables from CDF variables
    input_vars = calculate_inputs.calculate_inputs(cdf_vars)

    # Final conversion: Interpolate onto larger grid of points
    mmm_vars = convert_inputs.final_conversion(input_vars, input_options)

    return mmm_vars, input_vars, cdf_vars, raw_cdf_vars

def run_controller(input_options):
    # Initialize variable objects
    mmm_vars, input_vars, cdf_vars, raw_cdf_vars = initialize_controller(input_options)

    # (Optional) Create PDFs of various variable profile plots
    plot_profiles.plot_profile_comparison(cdf_vars, input_vars, input_options)
    plot_profiles.plot_input_profiles(mmm_vars, input_options)
    plot_profiles.plot_additional_profiles(mmm_vars, input_options)

    # Execute a basic MMM run to plot output profiles and base output values csv
    execute_basic_run(mmm_vars, input_options)

    # Execute a variable scan in MMM
    if input_options.var_to_scan is not None:
        execute_variable_scan(mmm_vars, input_options)

if __name__ == '__main__':
    '''
    CDF Options: Uncomment the line you wish to use
    '''
    cdf_name, shot_type, input_time = '129041A10', 'NSTX', 0.5
    # cdf_name, shot_type, input_time = '120982A09', 'NSTX', 0.5
    # cdf_name, shot_type, input_time = '132017T01', 'NSTX', 2.1
    # cdf_name, shot_type, input_time = '141552A01', 'NSTX', 2.1

    """
    Input Options:
    * input_points is the number of points to use when making the MMM input file
    * Set input_points = None to match the number of points used in the CDF
    * Set var_to_scan = var_name (str) to run a scan of the specified variable
    * Set var_to_scan = None to skip the variable scan
    * E.g.: var_to_scan = 'te'
    """
    input_options = variables.InputOptions(
        cdf_name=cdf_name,
        shot_type=shot_type,
        input_time=input_time,
        input_points=None,
        var_to_scan=None,
        scan_range=np.arange(start=0.5, stop=2.1, step=0.1))

    run_controller(input_options)
