# Standard Packages
import sys; sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np

# Local Packages
import modules.options as options
import modules.utils as utils
import modules.constants as constants
from modules.enums import SaveType


def read_from_files(file_list, dtype):
    '''
    Reads in data from a list of files.

    Parameters:
    * file_list (list): List of files to read data from
    * dtype (str or float): The data type in each file of file_list

    Returns:
    * (np.ndarray): Array of np.ndarrays, where each sub array contains data from a CSV
    '''

    data_list = []
    for file in file_list:
        data_list.append(np.genfromtxt(file, delimiter=',', dtype=dtype))

    return np.array(data_list)


def reshape_data(data_array, var_names):
    '''
    Reshapes data as a function of rho to data as a function of the scanned parameter

    Variable headers remain the same in the reshaped_data list as they are in data_array.
    Each array within reshaped_data corresponds to one value of rho.  For example,
    reshaped_data[0] will correspond to all variable data corresponding to the first value
    of rho, which is usually rho = 0.

    Parameters:
    * data_array (np.ndarray): Array of arrays, where each sub array corresponds to a CSV
    * var_names (list): List of names from the header of the CSV

    Returns:
    * reshaped_data (list): List of arrays, where each array is a different rho value
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


def save_reshaped_csv(reshaped_data, var_names, save_dir, save_type):
    '''
    Saves data reshaped as a function of the scanned parameter to CSV files.

    One CSV file is created for each value of rho (specified in each file name), for both input
    and output data.

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
        rho_value = constants.RHO_VALUE_FMT_STR.format(data[0, 0] / rmin_max_value)
        file_name = f'{base_file_name}{rho_value}.csv'
        np.savetxt(file_name, data, fmt='%.4e', delimiter=',', header=header_str)


def save_simple_csv(data, var_names, save_dir, save_type):
    '''
    Saves data to a CSV that didn't need to be reshaped

    Parameters:
    * reshaped_data (list): List of np.ndarray's, where each array corresponds to a CSV to save
    * var_names (list): List of variables names that will serve as the header to the CSV
    * save_dir (str): The path where the csv is to be saved
    * save_type (str): The name of the data type to be saved
    '''

    file_name = f'{save_dir}\\{save_type}'
    header_str = ','.join(var_names)
    np.savetxt(f'{file_name}.csv', data, fmt='%.4e', delimiter=',', header=header_str)


def parse_scan_csv():
    '''
    Parses all CSVs from a scan of var_to_scan and creates new CSVs for each rho point.

    The new CSVs are stored in './output/cdf_name/var_to_scan/rho', relative to the top level directory.
    Each CSV in the rho folder will correspond to one rho value of the scan, where data in these CSV will
    be a function of the scanned parameter.  When doing a control scan, only one CSV of controls will
    be created in the rho folder, since input controls are independent of rho.
    '''

    opts = options.instance
    save_dir = utils.get_rho_path(opts.runid, opts.scan_num, opts.var_to_scan)
    scanned_dir = utils.get_var_to_scan_path(opts.runid, opts.scan_num, opts.var_to_scan)

    save_types = [SaveType.INPUT, SaveType.ADDITIONAL, SaveType.OUTPUT]
    for save_type in save_types:
        saved_files = utils.get_files_in_dir(scanned_dir, f'{save_type.name.capitalize()}*')
        # Obtain negative factors by checking for negative signs in the value of the factor
        negative_factors = [file for file in saved_files if '-' in file.split(constants.SCAN_FACTOR_VALUE_SEPARATOR)[1]]
        non_negative_factors = [file for file in saved_files if file not in negative_factors]
        # Sort negative factors in reverse order (e.g., -6, -5, -4, etc.), then join with non negative factors
        saved_files = negative_factors[::-1] + non_negative_factors
        saved_data = read_from_files(saved_files, float)
        var_names = np.genfromtxt(saved_files[0], delimiter=',', names=True).dtype.names
        reshaped_data = reshape_data(saved_data, var_names)
        save_reshaped_csv(reshaped_data, var_names, save_dir, save_type.name.capitalize())

    '''
    Read in control data (if any):
    * Additional control data is only saved when a control value is scanned
    * Controls are unpacked differently than how input and output data were handled above
    * Control data does not need to be reshaped, since it is not a function of rho
    '''
    control_files = utils.get_files_in_dir(scanned_dir, f'{SaveType.CONTROLS.name.capitalize()}*', show_warning=False)
    if len(control_files) > 0:
        control_names_data = read_from_files(control_files, str)
        control_names = [name[0] for name in control_names_data[0]]
        control_data = np.array([[float(data[1]) for data in data_file] for data_file in control_names_data])
        save_simple_csv(control_data, control_names, save_dir, SaveType.CONTROLS.name.capitalize())

    print(f'Saved reshaped scan CSV to {save_dir}')


'''
For testing purposes:
* Options.instance.runid needs to match an existing ./output/runid folder
* Options.instance.scan_num needs to match an existing scan number in the ./output/runid folder
* Options.instance.var_to_scan needs to match an existing var_to_scan folder within the scan_num folder
* Options.instance.scan_range just needs to be some np.ndarray
'''
if __name__ == '__main__':
    opts = options.instance
    opts.runid = 'TEST'
    opts.scan_num = 129
    opts.var_to_scan = 'shear'
    opts.scan_range = np.arange(1)
    utils.clear_folder(utils.get_rho_path(opts.runid, opts.scan_num, opts.var_to_scan), '*.csv')
    parse_scan_csv()
