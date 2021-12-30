"""Runs MMM using the provided wrapper

This module controls all operation of MMM.  Input data is first written to an
input file in MMM format.  Next, the MMM wrapper is ran using a terminal
command, which produces an output CSV upon completion.  Afterwards, the
output data is read into an OutputVariables object.

TODO:
* This module can potentially be replaced by F2PY - Calling Fortran routines
  from Python, which would eliminate the overhead involved with reading and
  writing input and output files (the MMM wrapper also reads and writes
  files).  Although, it appears that MMM would need to be compiled with a
  fixed number of input points in order to be compatible with F2PY.
"""

# Standard Packages
import os
import subprocess

# Local Packages
import settings
import modules.utils as utils
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

    time_idx = input_vars.options.time_idx

    # Create input file in temp directory
    with open(_input_file, 'w') as f:
        f.write(controls.get_mmm_header())

        # Loop through MMM variables and write input variable labels and values
        var_names = input_vars.get_vars_of_type(SaveType.INPUT)
        for var_name in var_names:
            var = getattr(input_vars, var_name)
            units_str = f' [{var.units}]' if var.units else ''
            f.write(f'! {var.name}{units_str}\n')
            f.write(f'{var_name} = \n')

            values = var.values[:, time_idx]
            for value in values:
                f.write(f'   {value:{constants.INPUT_VARIABLE_VALUE_FMT}}\n')
            f.write('\n')

        f.write('/\n')  # Needed for the MMM wrapper to know that the input file has ended

    # Issue terminal command to run MMM
    result = subprocess.run(settings.MMM_DRIVER_PATH, cwd=_tmp_path,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)

    if settings.PRINT_MMM_RESPONSE:
        print(result.stdout)  # Only prints after MMM finishes running

    # Error checks
    if result.stderr:
        raise RuntimeError(result.stderr)
    if not os.path.exists(_output_file):
        raise FileNotFoundError('MMM did not produce an output file')
    if not os.stat(_output_file).st_size:
        raise ValueError('MMM produced an empty output file')

    output_vars = variables.OutputVariables(input_vars.options)
    output_vars.load_from_file_path(_output_file)
    os.remove(_output_file)  # ensure accurate error checks on next run

    return output_vars
