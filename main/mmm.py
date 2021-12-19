# Standard Packages
import os
import subprocess

# Local Packages
import settings
import main.variables as variables
from main.options import Options

MMM_LABELS = {
    'rmin': '! Half-width of the magnetic surface, r [m]',
    'rmaj': '! Major radius to geometric center of the magnetic surface, R [m]',
    'elong': '! Elongation profile',
    'ne': '! Electron density profile [m^-3]',
    'nh': '! Hydrogenic thermal particle density profile [m^-3]',
    'nz': '! Mean impurity density profile [m^-3]',
    'nf': '! Density from fast (non-thermal) ions profile [m^-3]',
    'zeff': '! Z_effective profile',
    'te': '! Electron temperature profile [keV]',
    'ti': '! Temperature of thermal ions profile [keV]',
    'q': '! Magnetic q profile',
    'btor': '! Toroidal magnetic field profile [Tesla]',
    'zimp': '! The profile of mean charge of impurities',
    'aimp': '! The profile of mean atomic mass of impurities',
    'ahyd': '! The profile of mean atomic mass of hydrogenic ions',
    'aimass': '! The profile of mean atomic mass of thermal ions',
    'wexbs': '! ExB shearing rate profile [rad/s]',
    'gne': '! Normalized electron density gradient profile',
    'gni': '! Normalized ion density gradient profile',
    'gnh': '! Normalized hydrogenic particle density gradient profile',
    'gnz': '! Normalized impurity density gradient profile',
    'gte': '! Normalized electron temperature gradient profile',
    'gti': '! Normalized ion temperature gradient profile',
    'gq': '! Normalized q gradient profile',
    'gvtor': '! Normalized toroidal velocity gradient profile',
    'vtor': '! Toroidal velocity profile [m/s]',
    'gvpol': '! Normalized poloidal velocity gradient profile',
    'vpol': '! Poloidal velocity profile [m/s]',
    'gvpar': '! Normalized parallel velocity gradient profile',
    'vpar': '! Parallel velocity profile [m/s]'
}


def run_driver(input_vars, controls):
    '''
    Controls operation of the MMM driver

    Steps:
    * Write input file to Cygwin home path
    * Runs MMM driver, which produces output file in Cygwin home path
    * Output file is read into OutputVariables object, then discarded

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
    for var_name in MMM_LABELS.keys():
        f.write(f'{MMM_LABELS[var_name]}\n')
        f.write(f'{var_name} = \n')

        values = getattr(input_vars, var_name).values[:, Options.instance.time_idx]
        for value in values:
            f.write('   {:.8e}\n'.format(value))
        f.write('\n')

    f.write('/\n')  # Needed for the MMM driver to know that the input file has ended
    f.close()

    # Shell command to run MMM driver through Cygwin
    cygwin_cmd = f'\"{settings.CYGWIN_BASH_PATH}\" -c -l {settings.MMM_DRIVER_PATH}'

    # Run MMM driver
    print('Running MMM Driver...')
    result = subprocess.run(cygwin_cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True, shell=True)

    # result.stdout is only meaningful if the MMM driver ran successfully
    print(result.stdout)

    # Error check: an output file may not be created when the MMM driver fails to run
    if not os.path.exists(output_file):
        raise FileNotFoundError('No output file produced: run the MMM Driver directly to determine the issue.')

    # Error check: bad input values can cause the MMM driver to produce a blank output file
    if os.stat(output_file).st_size == 0:
        raise ValueError('Output file was empty: run the MMM Driver directly to determine the issue.')

    # Read output.csv into output_vars
    output_vars = variables.OutputVariables()
    output_vars.load_from_file_path(output_file)

    # Delete output file from Cygwin folder (so we don't import old values if the MMM driver fails)
    os.remove(output_file)

    return output_vars
