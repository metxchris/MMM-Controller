# Standard Packages
import sys
sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt

# Local Packages
from main import utils, variables
from main.enums import ShotType, DataType
from main.options import Options
from main.input_controls import InputControls
from plots.styles import rho_layout as ps
from plots.colors import mmm


def run_plotting_loop(vars_to_plot):
    '''
    Creates PDF Plots of each variable in vars_to_plot

    The x-axis variable is pulled as var_to_scan from the saved options file

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    '''

    fig = plt.figure()
    runid = Options.instance.runid
    scan_num = Options.instance.scan_num
    var_to_scan = Options.instance.var_to_scan

    input_vars_dict, output_vars_dict = get_variable_dicts()

    for var_to_plot in vars_to_plot:
        print(f'Creating scanned variable PDF for {var_to_plot}...')

        # rho_values are the same for both input and output variables
        rho_values = input_vars_dict.keys()
        profile_type = f'{var_to_plot}_{var_to_scan}'
        
        for i, rho_value in enumerate(rho_values):
            sheet_num = '{:03d}'.format(i)

            xvar = getattr(input_vars_dict[rho_value], var_to_scan)
            yvar = getattr(output_vars_dict[rho_value], var_to_plot)

            plt.plot(xvar.values, yvar.values)
            plt.xlabel(f'{xvar.label}  {xvar.units_label}')
            plt.ylabel(f'{yvar.label}  {yvar.units_label}')
            plt.title(f'{yvar.label} vs. {xvar.label} ' + r' ($\rho = {0}$)'.format(rho_value))

            fig.savefig(utils.get_temp_path(f'{profile_type} {sheet_num}.pdf'))
            plt.clf() # Clear figure

        utils.open_file(utils.merge_profile_sheets(runid, scan_num, profile_type))

    # Clear figure from memory
    plt.close('all')

def load_variables_from_file(file_path, data_type):
    '''
    Loads data from a CSV into either an InputVariables or OutputVariables object

    Parameters:
    * file_path (str): The path to the CSV to load
    * data_type (DataType): The enum pertaining to the type of data in the CSV

    Returns:
    * loaded_vars (InputVariables or OutputVariables): Object containing data from the CSV
    '''

    data_array = np.genfromtxt(file_path, delimiter=',', dtype=float, names=True)
    var_names = data_array.dtype.names

    loaded_vars = None
    if data_type == DataType.INPUT:
        loaded_vars = variables.InputVariables()
    elif data_type == DataType.OUTPUT:
        loaded_vars = variables.OutputVariables()

    for var_name in var_names:
        getattr(loaded_vars, var_name).values = data_array[var_name]

    return loaded_vars

def get_rho_values(file_list):
    '''
    Returns a list of all rho values from the specified file list

    This extracts the value from files are named as 'rho = value.csv'

    Parameters:
    * file_list (list) List of files to check

    Returns:
    * (list) List of rho values as strings from the file list
    '''
    return [file.split('rho = ')[1].split('.csv')[0] for file in file_list]

def get_variable_dicts():
    '''
    Creates dictionaries that map rho values to InputVariables and OutputVariables objects

    Data is loaded from CSVs stored in the rho folder of the runid, scan_num, and var_to_scan
    currently stored in Options.instance. A list of rho values for the scan is created from 
    the filenames of the CSVs.

    Returns:
    * input_vars_dict (dict): Dictionary mapping rho values (str) to InputVariables data
    * output_vars_dict (dict): Dictionary mapping rho values (str) to OutputVariables data
    '''

    runid = Options.instance.runid
    scan_num = Options.instance.scan_num
    var_to_scan = Options.instance.var_to_scan

    # Get path to rho folder containing scanned data
    rho_path = utils.get_rho_path(runid, scan_num, var_to_scan)

    # Get input files from rho folder
    input_file_list = utils.get_files_in_dir(rho_path, DataType.INPUT.name.capitalize() + '*')

    # Get all values of rho from the file list
    rho_values = get_rho_values(input_file_list)

    input_vars_dict, output_vars_dict = {}, {}
    input_type_name = DataType.INPUT.name.capitalize()
    output_type_name = DataType.OUTPUT.name.capitalize()

    # Stores InputVariables and OutputVariables data objects for each rho_value
    for rho in rho_values:
        input_file_path = f'{rho_path}\\{input_type_name} rho = {rho}.csv'
        output_file_path = f'{rho_path}\\{output_type_name} rho = {rho}.csv'
        input_vars_dict[rho] = load_variables_from_file(input_file_path, DataType.INPUT)
        output_vars_dict[rho] = load_variables_from_file(output_file_path, DataType.OUTPUT)

    return input_vars_dict, output_vars_dict

def verify_vars_to_plot(vars_to_plot):
    '''
    Variables are verified to be members of OutputVariables before the plotting loop begins

    Parameters:
    * vars_to_plot (list):  List of output variables to plot
    '''

    output_vars = variables.OutputVariables()
    for var_to_plot in vars_to_plot:
        if not hasattr(output_vars, var_to_plot) and not hasattr(InputControls(), var_to_plot):
            raise NameError(f'Neither OutputVariables nor InputControls contain the variable named {var_to_plot}')

def main(vars_to_plot, runid, scan_num):
    '''
    Loads options, clears the temp folder, and verifies vars_to_plot, then runs the plotting loop

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * runid (str): The runid of the CDF
    * scan_num (int): The scan number to reference
    '''

    Options.instance.load_options(runid, scan_num)
    utils.clear_temp_folder()
    verify_vars_to_plot(vars_to_plot)
    run_plotting_loop(vars_to_plot)


# Run this file directly to plot scanned variable profiles from previously created scanned data
if __name__ == '__main__':
    '''
    Input Options:
    * vars_to_plot (list): List of output variables to plot
    * runid (str): The runid of the CDF
    * scan_num (int): The scan number to reference

    Examples:
    * vars_to_plot = ['xteMTM', 'xteETGM', 'xteETG', 'gmaMTM', 'omgMTM', 'dbsqprf']
    * vars_to_plot = variables.OutputVariables().get_vars_to_plot()
    '''
    vars_to_plot = ['kyrhoe_etgm']
    runid = '129041A10'
    scan_num = 45

    main(vars_to_plot, runid, scan_num)
