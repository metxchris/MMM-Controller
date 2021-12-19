# Standard Packages
import os
import subprocess

# Local Packages
import settings
import main.variables as variables
from main.options import Options
from main.enums import SaveType


def run_driver(input_vars, controls):
    '''
    Controls operation of the MMM driver

    Steps:
    * Write input file to Cygwin home path
    * Run MMM driver, which produces output file in Cygwin home path
    * Check that the output file exists and is not empty
    * Output file is read into OutputVariables object
    * Delete output file, to ensure that error checking is accurate on the next run

    Parameters:
    * input_vars (InputVariables): contains all data needed to write MMM input file
    * controls (InputControls): contains all data needed to write control values in the input file

    Returns:
    * output_vars (OutputVariables): contains all data read in from the MMM output file
    '''

    # Files are created relative to the Cygwin home path
    input_file = f'{settings.CYGWIN_HOME_PATH}input'
    output_file = f'{settings.CYGWIN_HOME_PATH}output.csv'

    # Create input file in Cygwin directory
    f = open(input_file, 'w')
    f.write(controls.get_mmm_header())

    # Loop through MMM variables and write input file labels and values
    var_names = input_vars.get_vars_of_type(SaveType.INPUT)
    for var_name in var_names:
        var = getattr(input_vars, var_name)
        units_str = f' [{var.units}]' if var.units != '' else ''
        f.write(f'! {var.name}{units_str}\n')
        f.write(f'{var_name} = \n')

        values = var.values[:, Options.instance.time_idx]
        for value in values:
            f.write('   {:.8e}\n'.format(value))
        f.write('\n')

    f.write('/\n')  # Needed for the MMM driver to know that the input file has ended
    f.close()

    # Shell/Terminal command to run MMM driver through Cygwin
    cygwin_cmd = f'\"{settings.CYGWIN_BASH_PATH}\" -c -l {settings.MMM_DRIVER_PATH}'

    # Run MMM driver
    print('Running MMM Driver...')
    result = subprocess.run(cygwin_cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True, shell=True)

    # result.stdout is only meaningful if the MMM driver ran successfully
    print(result.stdout)

    # Error check: an output file may not be created when the MMM driver fails to run
    if not os.path.exists(output_file):
        raise FileNotFoundError('No output file produced: run testmmm manually to determine the issue.')
    # Error check: bad input values can cause the MMM driver to produce a blank output file
    if os.stat(output_file).st_size == 0:
        raise ValueError('Output file was empty: run testmmm manually to determine the issue.')

    # Read then delete output.csv
    output_vars = variables.OutputVariables()
    output_vars.load_from_file_path(output_file)
    os.remove(output_file)

    return output_vars
