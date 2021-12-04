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

# Save input and output values to csv format
def save_data_csvs(data_input, data_output, vars_input, vars_output, units_input, units_output, input_options):
    # Set save_dir directory for a basic run
    if input_options.scan_factor_str is None:
        save_dir = utils.get_output_path(input_options.runid)
        file_name_output = f'{save_dir}\\{input_options.runid} Output Profiles.csv'
        file_name_input = f'{save_dir}\\{input_options.runid} Input Profiles.csv'
    # Set save_dir directory for a variable scan (creates an additional sub folder)
    else:
        save_dir = utils.get_output_path(f'{input_options.runid}\\{input_options.var_to_scan}')
        file_name_output = f'{save_dir}\\Output {input_options.var_to_scan} = {input_options.scan_factor_str}.csv'
        file_name_input = f'{save_dir}\\Input {input_options.var_to_scan} = {input_options.scan_factor_str}.csv'

    utils.create_directory(save_dir)

    # When doing a variable scan, clear save_dir directory at the start of each scan
    if input_options.scan_factor_str is not None and float(input_options.scan_factor_str) == input_options.scan_range.min():
        utils.clear_folder(save_dir, '*.csv')

    # Save each input and output data to a csv
    csv_header_input = ','.join(vars_input) + '\n' + ','.join(units_input)
    csv_header_output = ','.join(vars_output) + '\n' + ','.join(units_output)
    np.savetxt(file_name_input, data_input, header=csv_header_input, fmt='%.4e', delimiter=',')
    np.savetxt(file_name_output, data_output, header=csv_header_output, fmt='%.4e', delimiter=',')

# Returns a list of names for either vars or units. The first item is a comment character, so it is skipped
def split_names_str(names_str):
    return names_str.replace('    ', ' ').replace('   ', ' ').replace('  ', ' ').replace('\n', '').split(' ')[1:]

# Returns a list of values.  The first value is an empty string, so it is skipped
def split_values_str(values_str):
    return values_str.replace('   ', '  ').replace(' -', '  -').replace('\n', '').split('  ')[1:]

# Read output file from MMM driver and store values to OutputVariables object
# This code is incredibly hacky, and would be much cleaner if the MMM driver produced out CSV files
def read_output_file(input_options):
    output_vars = variables.OutputVariables()

    output_file = utils.get_temp_path('output')
    with open(output_file) as file:
        lines = file.readlines()

    # Start and end line numbers for input and output data
    data_start_input = num_input_comment_lines
    data_end_input = data_start_input + input_options.input_points - 1
    data_start_output = num_input_comment_lines + num_output_comment_lines + input_options.input_points
    data_end_output = data_start_output + input_options.input_points - 1

    # Store var names, unit values, and data
    vars_input = split_names_str(lines[data_start_input - 1])
    vars_output = split_names_str(lines[data_start_output - 1])
    units_input = split_names_str(lines[data_start_input - 2])
    units_output = split_names_str(lines[data_start_output - 2])
    data_lines_input = lines[data_start_input: data_end_input + 1]
    data_lines_output = lines[data_start_output: data_end_output + 1]

    data_input = np.array([split_values_str(line) for line in data_lines_input], dtype=float)
    data_output = np.array([split_values_str(line) for line in data_lines_output], dtype=float)

    # Save data to OutputVariables object (No need to store InputVariables again)
    for n, var in enumerate(vars_output):
        getattr(output_vars, var).set_variable(data_output[:, n], units_output[n], ['RMIN'])

    # Calculate and save rho for OutputVariables only
    output_vars.rho.set_variable(output_vars.rmin.values / output_vars.rmin.values[-1])

    # Save both input and output variables to CSV
    save_data_csvs(data_input, data_output, vars_input, vars_output, units_input, units_output, input_options)

    return output_vars

if __name__ == '__main__':
    # For testing purposes, make sure input_points is correct for existing output file in temp folder
    input_options = variables.InputOptions('132017T01', input_points=41)
    input_options.input_points = input_options.input_points
    input_options.runid = input_options.cdf_name
    read_output_file(input_options)
