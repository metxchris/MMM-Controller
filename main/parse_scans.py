# Standard Packages
import sys
sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np

# Local Packages
from main import utils, variables
from main.options import Options


def read_from_files(file_list):
    '''
    Reads in data from a list of files.

    Parameters:
    * file_list (list): List of files to read data from

    Returns:
    * (np.ndarray): Array of np.ndarrays, where each sub array contains data from a CSV
    '''

    data_list = [] 
    for file in file_list:
        data_list.append(np.genfromtxt(file, delimiter=',', dtype=float))

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
    * reshaped_data (list): List of arrays, where each array corresponds to a new CSV to write
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
    '''

    rmin_max_value = reshaped_data[-1][0, 0]
    base_file_name = f'{save_dir}\\{save_type} rho = '
    header_str = ','.join(var_names)

    for data in reshaped_data:
        rho_value = '{:.3f}'.format(data[0, 0] / rmin_max_value)
        file_name = f'{base_file_name}{rho_value}.csv'
        np.savetxt(file_name, data, fmt='%.4e', delimiter=',', header=header_str)

    print(f'Saved reshaped scan CSV to {save_dir}')

def parse_scan_csv():
    '''
    Parses all CSV from a scan of var_to_scan and creates new CSVs for each rho point.

    The new output CSV are stored in './output/cdf_name/var_to_scan/rho', relative to the top level directory.
    Each CSV in the rho folder will correspond to one rho value of the scan, where data in these CSV will
    be a function of the scanned parameter.  
    '''

    opts = Options.instance

    # Get list of input and output files from parameter scan
    scanned_dir = utils.get_var_to_scan_path(opts.runid, opts.scan_num, opts.var_to_scan)
    input_files = utils.get_files_in_dir(scanned_dir, 'Input*')
    output_files = utils.get_files_in_dir(scanned_dir, 'Output*')

    # Set and clear directory to save reshaped data to, if needed
    save_dir = utils.get_rho_path(opts.runid, opts.scan_num, opts.var_to_scan)
    utils.clear_folder(save_dir, '*.csv')

    # Read in data from parameter scan, then construct lists of new data to save
    input_data = read_from_files(input_files)
    output_data = read_from_files(output_files)
    input_var_names = np.genfromtxt(input_files[0], delimiter=',', names=True).dtype.names
    output_var_names = np.genfromtxt(output_files[0], delimiter=',', names=True).dtype.names
    reshaped_input_data = reshape_data(input_data, input_var_names)
    reshaped_output_data = reshape_data(output_data, output_var_names)

    # Save reshaped data to new CSVs
    save_reshaped_csv(reshaped_input_data, input_var_names, save_dir, 'Input')
    save_reshaped_csv(reshaped_output_data, output_var_names, save_dir, 'Output')


'''
For testing purposes:
* Options.instance.runid needs to match an existing ./output/runid folder
* Options.instance.scan_num needs to match an existing scan number in the ./output/runid folder
* Options.instance.var_to_scan needs to match an existing var_to_scan folder within the scan_num folder
* Options.instance.scan_range just needs to be some np.ndarray
'''
if __name__ == '__main__':
    Options.instance.runid = '129041A10'
    Options.instance.scan_num = 2
    Options.instance.var_to_scan = 'gte'
    Options.instance.scan_range = np.arange(1)
    parse_scan_csv()
