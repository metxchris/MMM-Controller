"""Runs MMM using the provided wrapper

This module controls all operation of MMM.  Input data is first written to an
input file in MMM format.  Next, the MMM driver is ran using a terminal
command, which produces an output CSV upon completion.  Afterwards, the
output data is read into an OutputVariables object.
"""

# Standard Packages
import os
import subprocess

# Local Packages
import settings
import modules.utils as utils
import modules.options as options
import modules.variables as variables
import modules.constants as constants
from modules.enums import SaveType


# Cache paths needed to run MMM
_tmp_path = utils.get_temp_path()
_input_file = utils.get_temp_path('input')  # input has no file type
_output_file = utils.get_temp_path('output.csv')


def run_wrapper(input_vars, controls):
    '''
    Controls operation of the MMM wrapper

    Steps:
    * Write input file to the temp folder
    * Run MMM wrapper, which produces output file in the temp folder
    * Check that the output file exists and is not empty
    * Output file is read into OutputVariables object
    * Delete output file, to ensure that error checking is accurate on the next run

    Parameters:
    * input_vars (InputVariables): contains all data needed to write MMM input file
    * controls (InputControls): contains all data needed to write control values in the input file

    Returns:
    * output_vars (OutputVariables): contains all data read in from the MMM output file

    Raises:
    * RuntimeError: If MMM has a runtime error
    * FileNotFoundError: If MMM does not produce an output file
    * ValueError: If MMM produces an empty output file
    '''

    # Create input file in temp directory
    f = open(_input_file, 'w')
    f.write(controls.get_mmm_header())

    # Loop through MMM variables and write input variable labels and values
    var_names = input_vars.get_vars_of_type(SaveType.INPUT)
    for var_name in var_names:
        var = getattr(input_vars, var_name)
        units_str = f' [{var.units}]' if var.units else ''
        f.write(f'! {var.name}{units_str}\n')
        f.write(f'{var_name} = \n')

        values = var.values[:, options.instance.time_idx]
        for value in values:
            f.write(f'   {constants.INPUT_VARIABLE_VALUE_FMT_STR}\n'.format(value))
        f.write('\n')

    f.write('/\n')  # Needed for the MMM wrapper to know that the input file has ended
    f.close()

    # Issue terminal command to run MMM
    result = subprocess.run(settings.MMM_DRIVER_PATH, cwd=_tmp_path,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True, shell=True)

    print(result.stdout)  # Only prints after MMM finishes running

    # Error checks
    if result.stderr:
        raise RuntimeError(result.stderr)
    if not os.path.exists(_output_file):
        raise FileNotFoundError('MMM did not produce an output file')
    if not os.stat(_output_file).st_size:
        raise ValueError('MMM produced an empty output file')

    output_vars = variables.OutputVariables()
    output_vars.load_from_file_path(_output_file)
    os.remove(_output_file)  # ensure accurate error checks on next run

    return output_vars
