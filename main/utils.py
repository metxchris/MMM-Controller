# Standard Packages
import os
import glob
from math import floor, log10

# Local Packages
import pdftk
import output
import temp
import cdfs
import main.variables as variables
import main.controls as controls
import main.constants as constants
import main.calculations as calculations
import main.conversions as conversions
import main.read_cdf as read_cdf
from main.enums import ScanType, SaveType, MergeType


def get_cdf_path(file_name):
    '''Returns (str): the path to specified CDF within the CDF folder'''
    return f'{os.path.dirname(cdfs.__file__)}\\{file_name}.CDF'


def get_temp_path(file_name=''):
    '''Returns (str): the path to the temp folder'''
    return f'{os.path.dirname(temp.__file__)}\\{file_name}'


def get_output_path(file_name=''):
    '''Returns (str): the path to the output folder'''
    return f'{os.path.dirname(output.__file__)}\\{file_name}'


def get_pdftk_path():
    '''Returns (str): the path to the pdftk executable'''
    return f'{os.path.dirname(pdftk.__file__)}\\pdftk.exe'


def get_scan_num_path(runid, scan_num):
    '''Returns (str): the path to the scan number folder'''
    return get_output_path(f'{runid}\\scan {scan_num}')


def get_merged_rho_path(runid, scan_num):
    '''Returns (str): the path to merged rho PDF for parameter scans'''
    return f'{get_scan_num_path(runid, scan_num)}\\Merged Rho'


def get_merged_profile_factors_path(runid, scan_num):
    '''Returns (str): the path to merged factors PDF for parameter scans'''
    return f'{get_scan_num_path(runid, scan_num)}\\Merged Profile Factors'


def get_var_to_scan_path(runid, scan_num, var_to_scan):
    '''Returns: (str) the path of the scanned variable'''
    return get_output_path(f'{runid}\\scan {scan_num}\\{var_to_scan}')


def get_rho_path(runid, scan_num, var_to_scan):
    '''Returns (str): the path of the rho folder'''
    return f'{get_var_to_scan_path(runid, scan_num, var_to_scan)}\\rho'


def get_rho_files(runid, scan_num, var_to_scan, save_type):
    '''Returns (list): all rho files of save_type in the rho folder'''
    return get_files_in_dir(get_rho_path(runid, scan_num, var_to_scan), f'{save_type.name.capitalize()}*')


def get_rho_values(runid, scan_num, var_to_scan, save_type):
    '''Returns (list): the rho values of all rho files in the rho folder'''
    rho_files = get_rho_files(runid, scan_num, var_to_scan, save_type)
    return [file.split(f'rho{constants.RHO_VALUE_SEPARATOR}')[1].split('.csv')[0] for file in rho_files]


def init_output_dirs(options):
    '''
    Initializes all output directories needed for storing output data

    Created Directories (relative to top-level directory):
    * ./output/runid/
    * ./output/runid/scan_num/
    * ./output/runid/scan_num/var_to_scan/
    * ./output/runid/scan_num/var_to_scan/rho/

    Parameters:
    * options (OptionsData): A reference to Options.instance
    '''

    if options.runid is None:
        raise ValueError('Cannot initialize output directories since the runid has not been set in Options')

    create_directory(get_output_path(options.runid))
    options.scan_num = set_scan_num(options.runid)

    if options.var_to_scan is not None:
        create_directory(get_var_to_scan_path(options.runid, options.scan_num, options.var_to_scan))
        create_directory(get_rho_path(options.runid, options.scan_num, options.var_to_scan))


def set_scan_num(runid):
    '''
    Initializes the directory for the current scan by always creating a new folder

    Parameters:
    * runid (str): The name of the CDF

    Returns:
    * scan_num (int): The chosen scan number
    '''

    num_range = range(1, 10000)

    for scan_num in num_range:
        scan_num_path = get_scan_num_path(runid, scan_num)
        if not os.path.exists(scan_num_path):
            create_directory(scan_num_path)
            break

    if scan_num == max(num_range):
        raise NameError(f'Maximum scan number reached {max(num_range)}! Clear some directories to continue')

    return scan_num


def create_directory(dir_name):
    '''
    Checks if output dir exists and creates it if needed

    Parameters:
    * dir_name (str): Path of directory
    '''

    if not os.path.exists(dir_name):
        os.mkdir(dir_name)


def check_filename(file_path, file_extension):
    '''
    Checks if file exists and returns a new file path if the checked file exists.

    A number in the form of (#) is appended to the end of the file path in the event that
    the file already exists.  This is done so that files are not overwritten in an existing directory.
    An exception is raised if too many duplicate files already exist.

    Parameters:
    * file_path (str): Path to where the file might exist
    * file_extension (str): The extension of the file (e.g.: '.pdf')

    Returns:
    * file_path (str): Path to file that does not exist
    '''

    if os.path.exists(file_path):
        num_range = range(2, 1000)
        for i in num_range:
            path_split = file_path.split(file_extension)
            new_file_path = f'{path_split[0]} ({i}){file_extension}'
            if not os.path.exists(new_file_path):
                file_path = new_file_path
                break

        if i == max(num_range):
            raise NameError(f'Too many duplicate files exist to save {file_path}')

    return file_path


def check_dirname(dir_path):
    '''
    Checks if directory exists and returns a new directory name if the checked directory exists.

    A number in the form of (#) is appended to the end of the directory name in the event that
    the directory already exists.  This is done so that files are not overwritten in an existing directory.
    An exception is raised if too many duplicate directories already exist.

    Parameters:
    * dir_path (str): Path to directory

    Returns:
    * dir_path (str): Path to directory that does not exist yet
    '''

    if os.path.isdir(dir_path):
        num_range = range(2, 1000)
        for i in num_range:
            new_dir_path = f'{dir_path} ({i})'
            if not os.path.isdir(new_dir_path):
                dir_path = new_dir_path
                break

        if i == max(num_range):
            raise NameError(f'Too many duplicate directories exist to save directory {dir_path}')

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
    print(f'Cleared all files of type {file_type} from {dir_path}\n')


def clear_temp_folder():
    '''Clears temporary files from the temp folder.'''

    # Clear individual pdf sheets
    temp_files = get_temp_path('*.pdf')
    for file in glob.glob(temp_files):
        os.remove(file)

    # Clear input and output files for MMM Driver
    if os.path.exists(get_temp_path('input')):
        os.remove(get_temp_path('input'))
    if os.path.exists(get_temp_path('output')):
        os.remove(get_temp_path('output'))


def get_files_in_dir(dir_path, file_type='', show_warning=True):
    '''
    Lists all files in dir_path of file_type.

    Parameters:
    * dir_path (str): Path of directory
    * file_type (str): Type of file to search for (include * in the string)

    Returns:
    * file_names (list): List of file names
    '''

    files = f'{dir_path}\\{file_type}'
    file_names = [file for file in glob.glob(files)]

    if len(file_names) == 0 and show_warning:
        print(f'*** Warning: No files found for {files}')

    return file_names


def merge_profile_sheets(runid, scan_num, profile_type, merge_type, var_to_scan=None, scan_factor=None):
    '''
    Merge PDF sheets using Pdftk in the temp folder into a single PDF, then place the merged PDF in the output folder.

    Pdftk is a 3rd party executable that is used to merge individual PDF sheets into one PDF,
    and is called using a shell command.

    Parameters:
    * runid (str): The name of the CDF
    * scan_num (int): The number of the scan
    * profile_type (str): The type of profile to merge (Input, Output, etc.)
    * merge_type (MergeType): The type of PDF merge
    * var_to_scan (str): The variable being scanned
    * scan_factor (float): The value of the scan factor

    Returns:
    * output_file (str): Path to merged PDF
    '''

    if merge_type == MergeType.PROFILES:
        output_path = get_scan_num_path(runid, scan_num)
        output_file = f'{output_path}\\{runid} {profile_type} Profiles.pdf'
    elif merge_type == MergeType.PROFILEFACTORS:
        output_path = get_merged_profile_factors_path(runid, scan_num)
        output_file = (f'{output_path}\\{runid} {profile_type} {var_to_scan}'
                       f'{constants.SCAN_FACTOR_VALUE_SEPARATOR}'
                       f'{constants.SCAN_FACTOR_PDF_FMT_STR.format(scan_factor)}.pdf')
    elif merge_type == MergeType.RHOVALUES:
        output_path = get_merged_rho_path(runid, scan_num)
        output_file = f'{output_path}\\{runid} {profile_type}.pdf'
    else:
        raise TypeError(f'No output path defined for {merge_type}')

    # Output directory creation only needed if sheets are being created outside of main mmm_controller.py execution
    create_directory(output_path)
    output_file = check_filename(output_file, '.pdf')
    temp_path = get_temp_path()
    pdftk_path = get_pdftk_path()

    # Shell command to use pdftk.exe
    # TODO: Replace os.system with subprocess.run()
    os.system(f'cd {temp_path} & {pdftk_path} *{profile_type}*.pdf cat output \"{output_file}\"')
    print(f'Profiles saved to \n    {output_file}\n')

    return output_file


def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    '''
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    '''

    if exponent is None:
        exponent = int(floor(log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if precision is None:
        precision = decimal_digits

    return r'${0:.{2}f}\times 10^{{{1:d}}}$'.format(coeff, exponent, precision)


def get_all_rho_data(runid, scan_num, var_to_scan):
    '''
    Creates dictionaries that map rho values to InputVariables and OutputVariables objects

    Data is loaded from CSVs stored in the rho folder of the runid, scan_num, and var_to_scan
    currently stored in Options.instance. A list of rho values for the scan is created from
    the filenames of the CSVs.

    Returns:
    * input_vars_dict (dict): Dictionary mapping rho values (str) to InputVariables input data
    * output_vars_dict (dict): Dictionary mapping rho values (str) to OutputVariables data
    * input_controls (InputControls or None): InputControls object with np.ndarray for values
    '''

    input_vars_dict, output_vars_dict = {}, {}
    rho_values = get_rho_values(runid, scan_num, var_to_scan, SaveType.OUTPUT)

    # Stores InputVariables and OutputVariables data objects for each rho_value
    for rho in rho_values:
        input_vars = variables.InputVariables()
        output_vars = variables.OutputVariables()

        args = (runid, scan_num, var_to_scan, None, rho)
        input_vars.load_from_csv(SaveType.INPUT, *args)
        input_vars.load_from_csv(SaveType.ADDITIONAL, *args)
        output_vars.load_from_csv(SaveType.OUTPUT, *args)

        input_vars_dict[rho] = input_vars
        output_vars_dict[rho] = output_vars

    # Get control_file from rho folder (there's at most one control file, as controls are independent of rho values)
    input_controls = controls.InputControls()
    input_controls.load_from_csv(runid, scan_num, var_to_scan, use_rho=True)

    return input_vars_dict, output_vars_dict, input_controls


def get_base_data(runid, scan_num):
    '''
    Gets all data pertaining to the base value of the scanned variable

    Returns:
    * input_vars (InputVariables): Object containing base input variable data
    * output_vars (OutputVariables): Object containing base output variable data
    * input_controls (InputControls): Object containing base input control data
    '''

    input_vars = variables.InputVariables()
    output_vars = variables.OutputVariables()
    input_controls = controls.InputControls()

    input_vars.load_from_csv(SaveType.INPUT, runid, scan_num)
    input_vars.load_from_csv(SaveType.ADDITIONAL, runid, scan_num)
    output_vars.load_from_csv(SaveType.OUTPUT, runid, scan_num)
    input_controls.load_from_csv(runid, scan_num)

    return input_vars, output_vars, input_controls


def initialize_variables():
    '''
    Initializes all input variables needed to run the MMM Driver and plot variable profiles

    Returns:
    * mmm_vars (InputVariables): All calculated variables, interpolated onto a grid of size input_points
    * cdf_vars (InputVariables): All CDF variables, interpolated onto a grid of size input_points
    * raw_cdf_vars (InputVariables): All unedited CDF variables (saved for troubleshooting)
    '''

    raw_cdf_vars = read_cdf.read_cdf()
    cdf_vars = conversions.convert_variables(raw_cdf_vars)
    mmm_vars = calculations.calculate_inputs(cdf_vars)

    return mmm_vars, cdf_vars, raw_cdf_vars
