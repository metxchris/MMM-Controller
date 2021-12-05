# Standard Packages
import os
import subprocess
import shutil
import sys
sys.path.insert(0, '../')

# Local Packages
from main import utils
import settings


def run_mmm_driver():
    # Shell command to run mmm driver through Cygwin
    cygwin_cmd = '\"{0}\" -c -l {1}{2}'.format(settings.CYGWIN_BASH_PATH, settings.CYGWIN_HOME_PATH, settings.MMM_DRIVER_PATH)

    # Copy input file from mmm directory to Cygwin directory
    shutil.copy(utils.get_temp_path('input'), settings.CYGWIN_HOME_PATH)

    # Run mmm driver
    print('Running MMM Driver...')
    result = subprocess.run(cygwin_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)

    # result.stdout is only meaningful if the mmm driver ran successfully 
    print(result.stdout)

    # Path of output file, if produced
    output_file = settings.CYGWIN_HOME_PATH + 'output'

    # Error check: no output file is produced when the mmm driver fails
    if not os.path.exists(output_file):
        raise FileNotFoundError('MMM Driver did not produce an output file.  Run the MMM Driver through bash to determine the issue.')
    
    # Copy output file from Cygwin directory to mmm directory
    shutil.copy(output_file, utils.get_temp_path())

    # Delete output file from Cygwin folder
    os.remove(output_file)


'''
For testing purposes:
* There needs to be an existing MMM input file in the temp folder
'''
if __name__ == '__main__':
    run_mmm_driver()
