"""Reshapes data stored in factor files into data stored in rho files

At the end of a variable scan after all factor CSV files are saved, reshaper
can be called to load in all data stored in each factor file, and then
reshape and save the data in new rho files.  In each factor file, variables
are dependent on the radial dimension rmin.  Reshaper redefines the
dependence of each variable instead to the scanned variable, and saves rho
values for each value of rho (from rmin) found in the factor files.

See the docstring on _reshape_data for an example of how this work.
"""

# Standard Packages
import sys; sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np

# Local Packages
import modules.utils as utils
import modules.constants as constants
from modules.enums import SaveType


def _read_from_files(file_list, dtype):
    '''
    Reads in data from a list of files.

    Parameters:
    * file_list (list): List of files to read data from
    * dtype (str or float): The data type in each file of file_list

    Returns:
    * (np.ndarray): Array of np.ndarrays; each sub array contains data from a CSV
    '''

    data_list = []
    for file in file_list:
        data_list.append(np.genfromtxt(file, delimiter=',', dtype=dtype))

    return np.array(data_list)


def _reshape_data(data_array, var_names):
    '''
    Reshapes data as a function of rho to data as a function of the scanned
    parameter

    Variable headers remain the same in the reshaped_data list as they are in
    data_array. Each array within reshaped_data corresponds to one value of
    rho.  For example, reshaped_data[0] will correspond to all variable data
    corresponding to the first value of rho, which is usually rho = 0.

    Example:
    * Arrays loaded from two factor CSVs as:
      - var1 factor = 1.csv      - var factor = 2.csv
        rho, var1, var2            rho, var1, var2
        0.1, 2.0, 3.0              0.1, 4.0, 3.0
        0.2, 2.5, 3.1              0.2, 5.0, 3.1
      Are reshaped into new arrays to be saved as:
      - var1 rho = 0.1.csv        - var1 rho = 0.2.csv
        rho, var1, var2             rho, var1, var2
        0.1, 2.0, 3.0               0.2, 2.5, 3.1
        0.1, 4.0, 3.0               0.2, 5.0, 3.1

    Parameters:
    * data_array (np.ndarray): Array of arrays, where each sub array corresponds to a CSV
    * var_names (list): List of names from the header of the CSV

    Returns:
    * reshaped_data (list): List of arrays; each array is a different rho value
    '''

    num_radial_points = data_array.shape[1]
    num_scan_factors = data_array.shape[0]
    num_variables = len(var_names)
    reshaped_data = []

    for r in range(num_radial_points):
        data_at_r = np.empty((num_scan_factors, num_variables))
        for f in range(num_scan_factors):
            data_at_r[f, :] = data_array[f, r]

        reshaped_data.append(data_at_r)

    return reshaped_data


def _save_reshaped_csv(reshaped_data, var_names, save_dir, save_type):
    '''
    Saves data reshaped as a function of the scanned parameter to CSV files.

    One CSV file is created for each value of rho (specified in each file
    name), for both input and output data.

    Parameters:
    * reshaped_data (list): List of np.ndarray's, where each array corresponds to a CSV to save
    * var_names (list): List of variables names that will serve as the header to the CSV
    * save_dir (str): The path where the csv is to be saved
    * save_type (str): The name of the data type to be saved
    '''

    rmin_max_value = reshaped_data[-1][0, 0]
    base_file_name = f'{save_dir}\\{save_type} rho{constants.RHO_VALUE_SEPARATOR}'
    header_str = ','.join(var_names)

    for data in reshaped_data:
        rho_value = f'{data[0, 0] / rmin_max_value:{constants.RHO_VALUE_FMT}}'
        file_name = f'{base_file_name}{rho_value}.csv'
        np.savetxt(file_name, data, fmt='%.4e', delimiter=',', header=header_str)


def _save_simple_csv(data, var_names, save_dir, save_type):
    '''
    Saves data to a CSV that didn't need to be reshaped

    Parameters:
    * reshaped_data (list): List of np.ndarray's, where each array corresponds to a CSV to save
    * var_names (list): List of variables names that will serve as the header to the CSV
    * save_dir (str): The path where the csv is to be saved
    * save_type (str): The name of the data type to be saved
    '''

    file_name = f'{save_dir}\\{save_type}.csv'
    header_str = ','.join(var_names)
    np.savetxt(file_name, data, fmt='%.4e', delimiter=',', header=header_str)


def create_rho_files(options):
    '''
    Parses all CSVs from a scan of var_to_scan and creates new CSVs for each
    rho point.

    The new CSVs are stored in './output/cdf_name/var_to_scan/rho', relative
    to the top level directory. Each CSV in the rho folder will correspond to
    one rho value of the scan, where data in these CSV will be a function of
    the scanned parameter.  When doing a control scan, only one CSV of
    controls will be created in the rho folder, since input controls are
    independent of rho.

    Parameters:
    * options (Options): Object containing user options
    '''

    save_dir = utils.get_rho_path(options.runid, options.scan_num, options.var_to_scan)
    scanned_dir = utils.get_var_to_scan_path(options.runid, options.scan_num, options.var_to_scan)

    save_types = [SaveType.INPUT, SaveType.ADDITIONAL, SaveType.OUTPUT]
    for save_type in save_types:
        saved_files = utils.get_files_in_dir(scanned_dir, f'{save_type.name.capitalize()}*')
        # Obtain negative factors by checking for negative signs in the value of the factor
        negative_factors = ([file for file in saved_files
                             if '-' in file.split(constants.SCAN_FACTOR_VALUE_SEPARATOR)[1]])
        non_negative_factors = [file for file in saved_files if file not in negative_factors]
        # Sort negative factors in reverse order (e.g., -6, -5, -4, etc.), then join with non negative factors
        saved_files = negative_factors[::-1] + non_negative_factors
        saved_data = _read_from_files(saved_files, float)
        var_names = np.genfromtxt(saved_files[0], delimiter=',', names=True).dtype.names
        reshaped_data = _reshape_data(saved_data, var_names)
        _save_reshaped_csv(reshaped_data, var_names, save_dir, save_type.name.capitalize())

    '''
    Read in control data (if any):
    * Additional control data is only saved when a control value is scanned
    * Controls are unpacked differently than how input and output data were handled above
    * Control data does not need to be reshaped, since it is not a function of rho
    '''
    control_files = utils.get_files_in_dir(scanned_dir, f'{SaveType.CONTROLS.name.capitalize()}*', show_warning=False)
    if len(control_files) > 0:
        control_names_data = _read_from_files(control_files, str)
        control_names = [name[0] for name in control_names_data[0]]
        control_data = np.array([[float(data[1]) for data in data_file] for data_file in control_names_data])
        _save_simple_csv(control_data, control_names, save_dir, SaveType.CONTROLS.name.capitalize())

    print(f'Saved reshaped scan CSV to {save_dir}')


'''
For testing purposes
'''
if __name__ == '__main__':
    from modules.options import Options
    runid = 'TEST'
    scan_num = 397
    var_to_scan = 'betae'
    options = Options(runid=runid, scan_num=scan_num, var_to_scan=var_to_scan)
    utils.clear_folder(utils.get_rho_path(runid, scan_num, var_to_scan), '*.csv')
    create_rho_files(options)
