"""Contains various utility functions used by the MMM controller package

Paths to all directories are stored in different functions here, in a
centralized location for other modules to use.  Adjusting the format of a
path here should not disrupt the creation of new scans, but will likely break
old scan data from being read. Several other directory and file related
operations are stored here as well.

Utils was written to be independent of classes that hold MMM related data,
other than the Options class.  See the datahelper module for utility type
functions that interface between the different data classes.
"""

# Standard Packages
import os
import glob
import logging
from math import floor, log10

# 3rd Party Packages
import numpy as np

# Local Packages
import pdftk
import output
import temp
import cdfs
import settings
import modules.constants as constants
from modules.enums import MergeType


_log = logging.getLogger(__name__)


def init_logging():
    '''Initializes logging based on settings'''
    if settings.PRINT_SAVE_MESSAGES:
        logging.basicConfig(level="INFO")
    else:
        logging.basicConfig(level="WARNING")


def get_cdf_path(file_name):
    '''Returns (str): the path to specified CDF within the CDF folder'''
    return f'{os.path.dirname(cdfs.__file__)}\\{file_name}.CDF'


def get_temp_path(file_name=''):
    '''Returns (str): the path to the temp folder'''
    return f'{os.path.dirname(temp.__file__)}\\{file_name}'


def get_pdftk_path():
    '''Returns (str): the path to the pdftk executable'''
    return f'{os.path.dirname(pdftk.__file__)}\\pdftk.exe'


def get_output_path():
    '''Returns (str): the path to the output folder'''
    return f'{os.path.dirname(output.__file__)}'


def get_runid_path(runid):
    '''Returns (str): the path to the runid folder'''
    return f'{get_output_path()}\\{runid}'


def get_scan_num_path(runid, scan_num):
    '''Returns (str): the path to the scan number folder'''
    return f'{get_runid_path(runid)}\\scan {scan_num}'


def get_options_path(runid, scan_num):
    '''Returns (str): the path to the options pickle file'''
    return f'{get_scan_num_path(runid, scan_num)}\\Options.pickle'


def get_merged_rho_path(runid, scan_num, var_to_scan):
    '''Returns (str): the path to merged rho PDF for parameter scans'''
    return f'{get_scan_num_path(runid, scan_num)}\\merged {var_to_scan} rho'


def get_merged_profile_factors_path(runid, scan_num):
    '''Returns (str): the path to merged factors PDF for parameter scans'''
    return f'{get_scan_num_path(runid, scan_num)}\\merged profile factors'


def get_var_to_scan_path(runid, scan_num, var_to_scan):
    '''Returns: (str) the path of the scanned variable'''
    return f'{get_scan_num_path(runid, scan_num)}\\{var_to_scan} factors'


def get_rho_path(runid, scan_num, var_to_scan):
    '''Returns (str): the path of the rho folder'''
    return f'{get_scan_num_path(runid, scan_num)}\\{var_to_scan} rho'


def get_rho_files(runid, scan_num, var_to_scan, save_type):
    '''Returns (list): all rho files of save_type in the rho folder'''
    return get_files_in_dir(get_rho_path(runid, scan_num, var_to_scan), f'{save_type.name.capitalize()}*')


def get_rho_strings(runid, scan_num, var_to_scan, save_type):
    '''Returns (list[str]): the rho values of all rho files in the rho folder as strings'''
    rho_files = get_rho_files(runid, scan_num, var_to_scan, save_type)
    return [file.split(f'rho{constants.RHO_VALUE_SEPARATOR}')[1].split('.csv')[0] for file in rho_files]


def get_closest_rho(options, save_type, rho_value):
    '''Returns (str): The actual saved rho value closest to the specified rho value'''
    runid = options.runid
    scan_num = options.scan_num
    var_to_scan = options.var_to_scan

    rho_values = np.array(get_rho_strings(runid, scan_num, var_to_scan, save_type), dtype=float)
    return f'{rho_values[np.argmin(np.abs(rho_values - float(rho_value)))]:{constants.RHO_VALUE_FMT}}'


def init_output_dirs(options):
    '''
    Initializes all output directories needed for storing output data

    Created Directories (relative to top-level directory):
    * ./output/runid/
    * ./output/runid/scan_num/
    * ./output/runid/scan_num/var_to_scan factors/
    * ./output/runid/scan_num/var_to_scan rho/

    Parameters:
    * options (Options): Object containing user options
    '''

    runid = options.runid
    scan_num = options.scan_num
    var_to_scan = options.var_to_scan

    clear_temp_folder()
    create_directory(get_runid_path(runid))
    create_directory(get_scan_num_path(runid, scan_num))

    if var_to_scan:
        create_directory(get_var_to_scan_path(runid, scan_num, var_to_scan))
        create_directory(get_rho_path(runid, scan_num, var_to_scan))


def get_scan_num(runid):
    '''
    Initializes the directory for the current scan by always creating a new folder

    Parameters:
    * runid (str): The name of the CDF

    Returns:
    * scan_num (int): The chosen scan number

    Raises:
    * ValueError: If there are too many scan number directories
    '''

    num_range = range(1, 10000)

    for scan_num in num_range:
        scan_num_path = get_scan_num_path(runid, scan_num)
        if not os.path.exists(scan_num_path):
            break

    if scan_num == max(num_range):
        raise ValueError(f'Maximum scan number reached {max(num_range)}! Clear some directories to continue')

    return scan_num


def create_directory(dir_name):
    '''
    Checks if output dir exists and creates it if needed

    Parameters:
    * dir_name (str): Path of directory
    '''

    if not os.path.exists(dir_name):
        os.mkdir(dir_name)


def check_exists(file_path):
    '''Returns (bool): True if the file exists'''
    return os.path.exists(file_path)


def check_filename(file_path, file_extension):
    '''
    Checks if file exists and returns a file path

    A number in the form of (#) is appended to the end of the file path in the
    event that the file already exists.  This is done so that files are not
    overwritten in an existing directory. An exception is raised if too many
    duplicate files already exist.

    Parameters:
    * file_path (str): Path to where the file might exist
    * file_extension (str): The extension of the file (e.g.: '.pdf')

    Returns:
    * file_path (str): Path to file that does not exist

    Raises:
    * ValueError: If too many duplicate files exist for the checked file
    '''

    if check_exists(file_path):
        num_range = range(2, 1000)
        for i in num_range:
            path_split = file_path.split(file_extension)
            new_file_path = f'{path_split[0]} ({i}){file_extension}'
            if not os.path.exists(new_file_path):
                file_path = new_file_path
                break

        if i == max(num_range):
            raise ValueError(f'Too many duplicate files exist to save {file_path}')

    return file_path


def check_dirname(dir_path):
    '''
    Checks if directory exists and returns a unique directory name

    A number in the form of (#) is appended to the end of the directory name
    in the event that the directory already exists.  This is done so that
    files are not overwritten in an existing directory. An exception is
    raised if too many duplicate directories already exist.

    Parameters:
    * dir_path (str): Path to directory

    Returns:
    * dir_path (str): Path to directory that does not exist yet

    Raises:
    * ValueError: If too many duplicate directories exist for the checked directory
    '''

    if os.path.isdir(dir_path):
        num_range = range(2, 1000)
        for i in num_range:
            new_dir_path = f'{dir_path} ({i})'
            if not os.path.isdir(new_dir_path):
                dir_path = new_dir_path
                break

        if i == max(num_range):
            raise ValueError(f'Too many duplicate directories exist to save directory {dir_path}')

    return dir_path


def open_file(file_path):
    '''
    Opens the specified file (likely only works on Windows)

    Parameters:
    * file_path (str): Path of file to open
    '''

    os.startfile(file_path)


def clear_folder(dir_path, file_type):
    '''
    Clears all files in a specified folder of the specified file_type.

    Parameters:
    * dir_path (str): Path of directory
    * file_type (str): Type of file to clear (include * in the string)
    '''

    folder = f'{dir_path}\\{file_type}'
    for file in glob.glob(folder):
        os.remove(file)
    _log.info(f'\n\tCleared all files of type {file_type} from {dir_path}\n')


def clear_temp_folder():
    '''Clears temporary files from the temp folder.'''
    temp_files = get_temp_path('*.pdf')
    for file in glob.glob(temp_files):
        os.remove(file)


def get_files_in_dir(dir_path, file_type='', show_warning=True):
    '''
    Lists all files in dir_path of file_type.

    Parameters:
    * dir_path (str): Path of directory
    * file_type (str): Type of file to search for (include * in the string)
    * show_warning (bool): Prints a warning if files aren't found (optional)

    Returns:
    * file_names (list): List of file names
    '''

    files = f'{dir_path}\\{file_type}'
    file_names = [file for file in glob.glob(files)]

    if len(file_names) == 0 and show_warning:
        _log.warning(f'No files found for {files}')

    return file_names


def merge_profile_sheets(options, profile_name, merge_type, scan_factor=None):
    '''
    Merges individual PDF sheets into a single PDF

    Pdftk (server edition) is a 3rd party executable that is used to merge
    individual PDF sheets into one PDF, and is called using a shell command.
    Pdftk must be correctly installed in the top-level "pdftk" folder in
    order for this function to work.  See the readme file in the pdftk folder
    for installation instructions for pdftk.exe.

    Parameters:
    * options (Options): Object containing user options
    * profile_name (str): The name of the profiles being merged
    * merge_type (MergeType): The type of PDF merge
    * scan_factor (float): The value of the scan factor

    Returns:
    * output_file (str): Path to merged PDF

    Raises:
    * NotImplementedError: If a merge_type does not have an output path
    '''

    runid = options.runid
    scan_num = options.scan_num
    var_to_scan = options.var_to_scan

    if merge_type == MergeType.PROFILES:
        output_path = get_scan_num_path(runid, scan_num)
        output_file = f'{output_path}\\{runid} {profile_name} Profiles.pdf'
    elif merge_type == MergeType.PROFILEFACTORS:
        output_path = get_merged_profile_factors_path(runid, scan_num)
        output_file = (f'{output_path}\\{runid} {profile_name} {var_to_scan}'
                       f'{constants.SCAN_FACTOR_VALUE_SEPARATOR}'
                       f'{scan_factor:{constants.SCAN_FACTOR_PDF_FMT}}.pdf')
    elif merge_type == MergeType.RHOVALUES:
        output_path = get_merged_rho_path(runid, scan_num, var_to_scan)
        output_file = f'{output_path}\\{runid} {profile_name}.pdf'
    else:
        raise NotImplementedError(f'No output path defined for {merge_type}')

    # Output directory creation only needed if sheets are being created
    # outside of mmm_controller.py execution

    create_directory(output_path)
    output_file = check_filename(output_file, '.pdf')
    temp_path = get_temp_path()
    pdftk_path = get_pdftk_path()

    # Shell command to use pdftk.exe
    # TODO: Replace os.system with subprocess.run()
    os.system(f'cd {temp_path} & {pdftk_path} *{profile_name}*.pdf cat output \"{output_file}\"')
    _log.info(f'\n\tSaved: {output_file}\n')

    return output_file


def get_sci_notation(number, precision=1):
    '''
    Converts a number into scientific notation for use with LaTeX formatting

    Example:
    * 1.233e4 converts to $1.2\\times 10^{4}$

    Parameters:
    * number (float): The number to convert into scientific notation
    * precision (int): The number of decimal digits to show (optional)

    Returns:
    * (str): The string representing the scientific notation of the number
    '''

    exponent = int(floor(log10(abs(number))))
    coeff = number / float(10**exponent)
    return f'${coeff:.{precision}f}\\times 10^{{{exponent:d}}}$'
