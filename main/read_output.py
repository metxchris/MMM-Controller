# Standard Packages
import time
import sys
sys.path.insert(0, '../')
# 3rd Party Packages
import numpy as np
# Local Packages
from main import utils, variables

# The number of comment lines in the output file before data starts
# These will need to be updated if the mmm output file format changes
num_input_comment_lines = 4
num_output_comment_lines = 3

# Read output file from MMM driver and store values to OutputVariables object
def read_output_file(input_options):
    output_vars = variables.OutputVariables()

    output_file = utils.get_temp_path('output')
    with open(output_file) as file:
        lines = file.readlines()

    # Start and end line numbers for output data
    data_start = num_input_comment_lines + num_output_comment_lines + input_options.interp_points
    data_end = data_start + input_options.interp_points - 1

    # Store lines of units, vars, and data
    vars_line = lines[data_start - 1]
    units_line = lines[data_start - 2]
    data_lines = lines[data_start: data_end + 1]

    # Strip extra characters from units and vars line, and drop the first item (comment character)
    vars_list = vars_line.replace('    ', ' ').replace('   ', ' ').replace('  ', ' ').replace('\n', '').split(' ')[1:]

    # The units list may have an extra item (due to an output formatting error), so we are matching the length of the variable list
    units_list = units_line.replace('    ', ' ').replace('   ', ' ').replace('  ', ' ').replace('\n', '').split(' ')[1: len(vars_list) + 1]

    # Split data lines into numpy array
    data_array = np.array([line.replace('   ', '').replace(' -', '  -').replace('\n', '').split('  ') 
        for line in data_lines], dtype=float)
    
    # Save data to OutputVariables object
    for n, var in enumerate(vars_list):
        getattr(output_vars, var).set_variable(data_array[:, n], units_list[n], ['RMIN'])

    # Calculate and save rho
    output_vars.rho.set_variable(output_vars.rmin.values / output_vars.rmin.values[-1])

    # Save output data to csv
    utils.create_directory(utils.get_output_path(input_options.runid))
    np.savetxt(utils.get_output_path('{0}\\{1} Output Profiles.csv'.format(input_options.runid, input_options.runid)), 
        data_array, header=','.join(vars_list) + '\n' + ','.join(units_list), fmt='%.4e', delimiter=',')

    return output_vars

if __name__ == '__main__':
    # For testing purposes
    input_options = variables.InputOptions('132017T01', input_points=51)
    input_options.interp_points = input_options.input_points
    input_options.runid = input_options.cdf_name
    read_output_file(input_options)
