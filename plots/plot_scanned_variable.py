# Standard Packages
import sys
sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt

# Local Packages
from main import utils, constants
from main.enums import DataType, ScanType, SaveType
from main.options import Options
from main.controls import InputControls
from main.variables import InputVariables, OutputVariables
from plots.styles import rho_layout as ps
from plots.colors import mmm


def plot_parameter_scan(vars_to_plot):
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
    scan_type = Options.instance.scan_type

    input_vars_dict, additional_vars_dict, output_vars_dict, input_controls = get_scanned_data()

    var_to_scan_type = getattr(InputVariables(), var_to_scan).save_type
    x_vars_dict = input_vars_dict if var_to_scan_type == SaveType.INPUT else additional_vars_dict

    for var_to_plot in vars_to_plot:
        print(f'Creating scanned variable PDF for {var_to_plot}...')

        # rho_values are the same for both input and output variables
        rho_values = x_vars_dict.keys()
        profile_type = f'{var_to_plot}_{var_to_scan}'

        for i, rho_value in enumerate(rho_values):
            sheet_num = constants.SHEET_NUM_FMT_STR.format(i)

            xvar = getattr(x_vars_dict[rho_value] if scan_type == ScanType.VARIABLE else input_controls, var_to_scan)
            yvar = getattr(output_vars_dict[rho_value], var_to_plot)

            plt.plot(xvar.values, yvar.values)
            plt.xlabel(f'{xvar.label}  {xvar.units_label}')
            plt.ylabel(f'{yvar.label}  {yvar.units_label}')
            plt.title(f'{yvar.name}' + r' ($\rho = {0}$)'.format(rho_value))

            fig.savefig(utils.get_temp_path(f'{profile_type} {sheet_num}.pdf'))
            plt.clf()  # Clear figure data

        utils.open_file(utils.merge_profile_sheets(runid, scan_num, profile_type, is_scan=True))

    # Remove figure from memory
    plt.close('all')


def load_variables_from_file(file_path, data_type):
    '''
    Loads data from a CSV into either an InputVariables or OutputVariables object

    Parameters:
    * file_path (str): The path to the CSV to load
    * data_type (DataType): The enum pertaining to the type of data in the CSV

    Returns:
    * data_object (InputVariables, OutputVariables, or InputControls): Object containing data from the CSV
    '''

    data_array = np.genfromtxt(file_path, delimiter=',', dtype=float, names=True)
    var_names = data_array.dtype.names

    data_object = None
    if data_type == DataType.INPUT or data_type == DataType.ADDITIONAL:
        data_object = InputVariables()
    elif data_type == DataType.OUTPUT:
        data_object = OutputVariables()
    elif data_type == DataType.CONTROL:
        data_object = InputControls()

    for var_name in var_names:
        getattr(data_object, var_name).values = data_array[var_name]

    return data_object


def get_scanned_data():
    '''
    Creates dictionaries that map rho values to InputVariables and OutputVariables objects

    Data is loaded from CSVs stored in the rho folder of the runid, scan_num, and var_to_scan
    currently stored in Options.instance. A list of rho values for the scan is created from
    the filenames of the CSVs.

    Returns:
    * input_vars_dict (dict): Dictionary mapping rho values (str) to InputVariables input data
    * input_vars_dict (dict): Dictionary mapping rho values (str) to InputVariables additional data
    * output_vars_dict (dict): Dictionary mapping rho values (str) to OutputVariables data
    * input_controls (InputControls or None): InputControls object with np.ndarrays for values
    '''

    input_controls = None
    input_vars_dict, additional_vars_dict, output_vars_dict = {}, {}, {}
    input_type_name = DataType.INPUT.name.capitalize()
    additional_type_name = DataType.ADDITIONAL.name.capitalize()
    output_type_name = DataType.OUTPUT.name.capitalize()

    rho_path = get_rho_path()
    rho_values = get_rho_values(rho_path)

    # Stores InputVariables and OutputVariables data objects for each rho_value
    for rho in rho_values:
        input_file_path = f'{rho_path}\\{input_type_name} rho = {rho}.csv'
        additional_file_path = f'{rho_path}\\{additional_type_name} rho = {rho}.csv'
        output_file_path = f'{rho_path}\\{output_type_name} rho = {rho}.csv'
        input_vars_dict[rho] = load_variables_from_file(input_file_path, DataType.INPUT)
        additional_vars_dict[rho] = load_variables_from_file(additional_file_path, DataType.ADDITIONAL)
        output_vars_dict[rho] = load_variables_from_file(output_file_path, DataType.OUTPUT)

    # Get control_file from rho folder (there should be at most one control file)
    control_file = utils.get_files_in_dir(rho_path, DataType.CONTROL.name.capitalize() + '*', show_warning=False)

    for file in control_file:
        input_controls = load_variables_from_file(file, DataType.CONTROL)

    return input_vars_dict, additional_vars_dict, output_vars_dict, input_controls


def get_rho_path():
    '''Returns a path to the rho folder for this runid, scan_num, and var_to_scan (str)'''

    runid = Options.instance.runid
    scan_num = Options.instance.scan_num
    var_to_scan = Options.instance.var_to_scan
    return utils.get_rho_path(runid, scan_num, var_to_scan)


def get_rho_values(rho_path):
    '''
    Returns a list of all rho values used for output files that were saved in the rho_path (list)

    Parameters:
    * rho_path (str): Path to the rho folder
    '''

    output_files = utils.get_files_in_dir(rho_path, DataType.OUTPUT.name.capitalize() + '*')
    return [file.split('rho = ')[1].split('.csv')[0] for file in output_files]


def verify_vars_to_plot(vars_to_plot):
    '''
    Variables are verified to be members of OutputVariables before the plotting loop begins

    Parameters:
    * vars_to_plot (list):  List of output variables to plot
    '''

    output_vars = OutputVariables()
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

    print('Initializing data...')
    utils.clear_temp_folder()
    Options.instance.load_options(runid, scan_num)
    verify_vars_to_plot(vars_to_plot)
    plot_parameter_scan(vars_to_plot)


# Run this file directly to plot scanned variable profiles from previously created scanned data
if __name__ == '__main__':
    '''
    Input Options:
    * vars_to_plot (list): List of output variables to plot
    * runid (str): The runid of the CDF
    * scan_num (int): The scan number to reference

    Examples:
    * vars_to_plot = ['xteMTM', 'xteETGM', 'xteETG', 'gmaMTM', 'omgMTM', 'dbsqprf']
    * vars_to_plot = OutputVariables().get_all_output_vars()
    * vars_to_plot = OutputVariables().get_etgm_vars()
    '''
    vars_to_plot = ['xteETGM', 'gmaETGM', 'omgETGM']
    # runid = '120968A02'
    # runid = '120982A09'
    # runid = '129041A10'
    runid = 'TEST'
    scan_num = 19

    main(vars_to_plot, runid, scan_num)
