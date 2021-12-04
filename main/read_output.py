# Standard Packages
import sys
sys.path.insert(0, '../')
# 3rd Party Packages
import numpy as np
# Local Packages
from main import utils, variables

# The number of comment lines in the output file before data starts
# These will need to be updated if the mmm output file format changes
NUM_INPUT_COMMENT_LINES = 4
NUM_OUTPUT_COMMENT_LINES = 3

def save_data_csvs(data_input, data_output, vars_input, vars_output, units_input, units_output, input_options):
    '''
    Save input and output values from the MMM driver output file to separate CSVs.

    The output directories (save_dir) are given as follows, where both are relative to the top level directory.
    * For a basic MMM run, the output directory is './output/cdf_name'.
    * When doing parameter scans, the output directory is './output/cdf_name/var_to_scan'.

    The chosen output directory is cleared at the start of each scan, so that any files there will only
    be created as a result of the most recent scan.

    Parameters:
    * data_input (np.ndarray): Array of all input values
    * data_output (np.ndarray): Array of all output values
    * vars_input (list): List of input variable names
    * vars_output (list): List of output variable names
    * units_input (list): List of input units
    * units_output (list): List of output units
    * input_options (InputOptions): Stores options for the scan
    '''

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

    # Creates two header rows in each CSV
    csv_header_input = ','.join(vars_input) + '\n' + ','.join(units_input)
    csv_header_output = ','.join(vars_output) + '\n' + ','.join(units_output)
    np.savetxt(file_name_input, data_input, header=csv_header_input, fmt='%.4e', delimiter=',')
    np.savetxt(file_name_output, data_output, header=csv_header_output, fmt='%.4e', delimiter=',')
 
def split_names(names_str):
    '''
    Returns a list of names for either vars or units. The first item is a comment character, so it is skipped

    Parameters:
    * names_str (str): String of either variable or unit names from the output file

    Returns:
    * (list): List of split names from names_str
    '''

    return names_str.replace('    ', ' ').replace('   ', ' ').replace('  ', ' ').replace('\n', '').split(' ')[1:]

def split_values(values_str):
    '''
    Returns a list of values.  The first value is an empty string, so it is skipped

    Parameters:
    * names_str (str): String of variable values from the output file

    Returns:
    * (list): List of split variable values from names_str
    '''

    return values_str.replace('   ', '  ').replace(' -', '  -').replace('\n', '').split('  ')[1:]

def read_output_file(input_options):
    '''
    Read output file from MMM driver and store values

    Output data from the output file is stored to an OutputVariables object. Both input and output 
    data are also written to CSV, so that they can later be quickly parsed when plotting parameter scans.
    Note: This code is incredibly hacky, and would be much cleaner if the MMM driver produced out CSV files

    Parameters:
    * input_options (InputOptions): Stores options for the scan

    Returns:
    * output_vars (OutputVariables): All output variable data read from the output file of the MMM driver
    '''

    output_vars = variables.OutputVariables()

    output_file = utils.get_temp_path('output')
    with open(output_file) as file:
        lines = file.readlines()

    # Start and end line numbers for input and output data
    data_start_input = NUM_INPUT_COMMENT_LINES
    data_end_input = data_start_input + input_options.input_points - 1
    data_start_output = NUM_INPUT_COMMENT_LINES + NUM_OUTPUT_COMMENT_LINES + input_options.input_points
    data_end_output = data_start_output + input_options.input_points - 1

    # Store var names, unit values, and data
    # TODO: Input var units are not being parsed correctly due to formatting issues in the input file
    vars_input = split_names(lines[data_start_input - 1])
    vars_output = split_names(lines[data_start_output - 1])
    units_input = split_names(lines[data_start_input - 2])
    units_output = split_names(lines[data_start_output - 2])
    data_lines_input = lines[data_start_input: data_end_input + 1]
    data_lines_output = lines[data_start_output: data_end_output + 1]

    data_input = np.array([split_values(line) for line in data_lines_input], dtype=float)
    data_output = np.array([split_values(line) for line in data_lines_output], dtype=float)

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
    input_options.runid = input_options.cdf_name
    read_output_file(input_options)
