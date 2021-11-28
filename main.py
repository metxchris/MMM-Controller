# 3rd Party Packages
import numpy as np
# Local Packages
from main import read_cdf, convert_inputs, calculate_inputs, variables, constants, utils, write_inputs, run_driver
from plots import plot_input_profiles

def main(input_options):
    # Clear temp folder
    utils.clear_temp_folder()

    # Read variables from specified CDF
    cdf_vars = read_cdf.read_cdf(input_options)

    # Initial conversion variables from CDF format to MMM format
    input_vars = convert_inputs.initial_conversion(cdf_vars, input_options)

    # TODO: add step to vary input_vars over a specified range

    # Calculate new variables from CDF variables
    calculate_inputs.calculate_inputs(input_vars)

    # Final conversion: Interpolate onto larger grid of points
    convert_inputs.final_conversion(input_vars, input_options)

    # Plot input profiles being sent to the MMM driver and save as PDF
    plot_input_profiles.make_plots(input_vars, input_options)

    # Write variables into input for MMM Driver
    write_inputs.write_input_file(input_vars, input_options)

    # Run MMM driver to produce output file
    run_driver.run_mmm_driver(input_options)

    # TODO: read results from output file
    
    # input_vars.print_nonzero_variables()

if __name__ == '__main__':
    cdf_name = '132017T01'
    shot_type = 'DIII-D'
    input_time = 2.1
    input_points = 201

    main(variables.InputOptions(cdf_name, shot_type, input_time, input_points))
