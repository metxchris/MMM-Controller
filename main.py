# Standard Packages
import time
# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt
# Local Packages
from mmm_package import read_cdf, convert_inputs, calculate_inputs, variables, constants, utils
from tests import test
from plots import plot2d, plot_input_profiles

def main(input_options):
    # Clear temp folder
    utils.clear_temp_folder()

    # Read variables from specified CDF
    cdf_vars = read_cdf.read_cdf(input_options)

    # Convert variables from CDF format to MMM format
    input_vars = convert_inputs.convert_inputs(cdf_vars, input_options)

    # Calculate new variables from CDF variables
    calculate_inputs.calculate_inputs(input_vars)

    # TODO: add step to vary input_vars over a specified range
    # input_vars.print_nonzero_variables()

    # Plot input profiles being sent to the MMM driver and save as PDF
    plot_input_profiles.make_plots(input_vars, input_options)

if __name__ == '__main__':
    cdf_name = '132017T01'
    shot_type = 'DIII-D'
    input_time = 2.1
    input_points = 200

    main(variables.InputOptions(cdf_name, shot_type, input_time, input_points))
