# Standard Packages
import sys; sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np

# Local Packages
from main import utils, variables
from main.options import Options
from main.enums import SaveType


# The number of comment lines in the output file before data starts
# These will need to be updated if the MMM output file format changes
NUM_INPUT_COMMENT_LINES = 4
NUM_OUTPUT_COMMENT_LINES = 3


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
    Read output file from MMM driver and store values in an OutputVariables object

    Returns:
    * output_vars (OutputVariables): All output variable data read from the output file of the MMM driver
    '''

    opts = Options.instance
    output_vars = variables.OutputVariables()

    output_file = utils.get_temp_path(SaveType.OUTPUT.name.lower())
    with open(output_file) as file:
        lines = file.readlines()

    # Start and end line numbers for input and output data
    data_start_output = NUM_INPUT_COMMENT_LINES + NUM_OUTPUT_COMMENT_LINES + opts.input_points
    data_end_output = data_start_output + opts.input_points - 1

    # Store var names, unit values, and data
    vars_output = split_names(lines[data_start_output - 1])
    units_output = split_names(lines[data_start_output - 2])
    data_lines_output = lines[data_start_output: data_end_output + 1]
    data_output = np.array([split_values(line) for line in data_lines_output], dtype=float)

    # Save data to OutputVariables object
    for n, var in enumerate(vars_output):
        getattr(output_vars, var).set_variable(data_output[:, n], units_output[n], ['RMIN'])

    output_vars.set_rho_values()

    return output_vars


'''
For testing purposes:
* There needs to be an existing MMM output file in the temp folder
* Options.instance.runid needs to match an existing ./output/runid folder
* Options.instance.scan_num needs to match an existing scan number in the ./output/runid folder
* Options.instance.input_points needs to match points used in existing output file in the temp folder
'''
if __name__ == '__main__':
    Options.instance.runid = 'TEST'
    Options.instance.input_points = 200
    Options.instance.scan_num = 1
    output_vars = read_output_file()
    output_vars.print_nonzero_variables()
