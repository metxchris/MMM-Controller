# Standard Packages
import sys
sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np

# Local Packages
from main import utils, variables
from main.options import Options
from main.enums import DataType


# The number of comment lines in the output file before data starts
# These will need to be updated if the mmm output file format changes
NUM_INPUT_COMMENT_LINES = 4
NUM_OUTPUT_COMMENT_LINES = 3


def save_data_csvs(scan_factor, data_input, data_output, vars_input, vars_output, units_input, units_output):
    '''
    Save input and output values from the MMM driver output file to separate CSVs.

    The output directories (save_dir) are given as follows, where both are relative to the top level directory.
    * For a basic MMM run, the output directory is './output/cdf_name'.
    * When doing parameter scans, the output directory is './output/cdf_name/var_to_scan'.

    The chosen output directory is cleared at the start of each scan, so that any files there will only
    be created as a result of the most recent scan.

    Parameters:
    * scan_factor (float or None): The current factor when doing a parameter scan
    * data_input (np.ndarray): Array of all input values
    * data_output (np.ndarray): Array of all output values
    * vars_input (list): List of input variable names
    * vars_output (list): List of output variable names
    * units_input (list): List of input units
    * units_output (list): List of output units
    '''

    opts = Options.instance
    output_str = DataType.OUTPUT.name.capitalize()
    input_str = DataType.INPUT.name.capitalize()

    # Set save_dir directory for a basic run
    if scan_factor is None:
        save_dir = utils.get_scan_num_path(opts.runid, opts.scan_num)
        file_name_output = f'{save_dir}\\{opts.runid} {output_str} Profiles.csv'
        file_name_input = f'{save_dir}\\{opts.runid} {input_str} Profiles.csv'
    # Set save_dir directory for a variable scan (creates an additional sub folder)
    else:
        scan_factor_str = '{:.3f}'.format(scan_factor)
        save_dir = utils.get_var_to_scan_path(opts.runid, opts.scan_num, opts.var_to_scan)
        file_name_output = f'{save_dir}\\{output_str} {opts.var_to_scan} = {scan_factor_str}.csv'
        file_name_input = f'{save_dir}\\{input_str} {opts.var_to_scan} = {scan_factor_str}.csv'

    # When doing a variable scan, clear save_dir directory at the start of each scan
    if scan_factor is not None and scan_factor == opts.scan_range.min():
        utils.clear_folder(save_dir, '*.csv')

    # Creates two header rows in each CSV
    csv_header_input = ','.join(vars_input) + '\n' + ','.join(units_input)
    csv_header_output = ','.join(vars_output) + '\n' + ','.join(units_output)
    np.savetxt(file_name_input, data_input, header=csv_header_input, fmt='%.4e', delimiter=',')
    np.savetxt(file_name_output, data_output, header=csv_header_output, fmt='%.4e', delimiter=',')

    print(f'Output data saved to:\n    {file_name_input}\n    {file_name_output}\n')
 
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

def read_output_file(scan_factor=None):
    '''
    Read output file from MMM driver and store values

    Output data from the output file is stored to an OutputVariables object. Both input and output 
    data are also written to CSV, so that they can later be quickly parsed when plotting parameter scans.
    Note: This code is incredibly hacky, and would be much cleaner if the MMM driver produced out CSV files

    Returns:
    * output_vars (OutputVariables): All output variable data read from the output file of the MMM driver
    '''

    opts = Options.instance
    output_vars = variables.OutputVariables()

    output_file = utils.get_temp_path(DataType.OUTPUT.name.lower())
    with open(output_file) as file:
        lines = file.readlines()

    # Start and end line numbers for input and output data
    data_start_input = NUM_INPUT_COMMENT_LINES
    data_end_input = data_start_input + opts.input_points - 1
    data_start_output = NUM_INPUT_COMMENT_LINES + NUM_OUTPUT_COMMENT_LINES + opts.input_points
    data_end_output = data_start_output + opts.input_points - 1

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
    save_data_csvs(scan_factor, data_input, data_output, vars_input, vars_output, units_input, units_output)

    return output_vars


'''
For testing purposes:
* There needs to be an existing MMM output file in the temp folder
* Options.instance.runid needs to match an existing ./output/runid folder
* Options.instance.scan_num needs to match an existing scan number in the ./output/runid folder
* Options.instance.input_points needs to match points used in existing output file in the temp folder
'''
if __name__ == '__main__':
    Options.instance.runid = '129041A10'
    Options.instance.input_points = 101
    Options.instance.scan_num = 1
    read_output_file()
