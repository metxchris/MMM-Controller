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

# Save output values to csv format
def save_output_csv(data_array, vars_list, units_list, input_options):
    # Set output directory for a basic run
    if input_options.scan_factor_str is None:
        output_dir = utils.get_output_path(input_options.runid)
        file_name = '{0}\\{1} Output Profiles.csv'.format(output_dir, input_options.runid)
    # Set output directory for a variable scan (creates an additional sub folder)
    else:
        output_dir = utils.get_output_path('{0}\\{1}'.format(input_options.runid, input_options.var_to_scan))
        file_name = '{0}\\{1} = {2}.csv'.format(output_dir, input_options.var_to_scan, input_options.scan_factor_str)

    utils.create_directory(output_dir)

    # When doing a variable scan, clear output directory at the start of each scan
    if input_options.scan_factor_str is not None and float(input_options.scan_factor_str) == input_options.scan_range.min():
        utils.clear_folder(output_dir, '*.csv')

    # Save output data to csv
    csv_header = ','.join(vars_list) + '\n' + ','.join(units_list)
    np.savetxt(file_name, data_array, header=csv_header, fmt='%.4e', delimiter=',')

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

    save_output_csv(data_array, vars_list, units_list, input_options)

    #TODO: Save input vars as well to CSV

    return output_vars

if __name__ == '__main__':
    # For testing purposes, make sure input_points is correct for existing output file in temp folder
    input_options = variables.InputOptions('132017T01', input_points=51)
    input_options.interp_points = input_options.input_points
    input_options.runid = input_options.cdf_name
    read_output_file(input_options)
