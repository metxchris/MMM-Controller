# 3rd Party Packages
import numpy as np
# Local Packages
from main import *
from plots import plot_profiles, plot2d

def main(input_options):
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

    # TODO: add step to vary input_vars over a specified range

    # Plot input profiles being sent to the MMM driver and save as PDF
    plot_profiles.plot_input_profiles(input_vars, input_options)

    # Write variables to input file for MMM Driver
    write_inputs.write_input_file(input_vars, input_options)

    # Run MMM driver to produce output file
    run_driver.run_mmm_driver(input_options)

    # Read output variables from output file
    output_vars = read_output.read_output_file(input_options)

    # Plot output profiles
    plot_profiles.plot_output_profiles(output_vars, input_options)
    
    # plot2d.plot(output_vars.rho.values, output_vars.xtiW20.values)
    # input_vars.print_nonzero_variables()

if __name__ == '__main__':
    cdf_name = '132017T01'
    shot_type = 'DIII-D'
    input_time = 2.1
    input_points = 51

    main(variables.InputOptions(cdf_name, shot_type, input_time, input_points))
