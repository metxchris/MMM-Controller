# Standard Packages
import sys
sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt

# Local Packages
import settings
from main import utils, constants
from main.enums import ScanType, SaveType
from main.options import Options
from main.controls import InputControls
from main.variables import InputVariables, OutputVariables
from plots.styles import single as plotlayout
from plots.colors import mmm as plotcolors


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
    scan_type = Options.instance.scan_type

    input_vars_dict, output_vars_dict, input_controls = get_scanned_data()

    for var_to_plot in vars_to_plot:
        print(f'Creating scanned variable PDF for {var_to_plot}...')

        # rho_strs are the same for both input and output variables
        rho_strs = input_vars_dict.keys()
        profile_type = f'{var_to_plot}_{var_to_scan}'

        for i, rho_str in enumerate(rho_strs):
            sheet_num = constants.SHEET_NUM_FMT_STR.format(i)

            xvar_data = input_vars_dict[rho_str] if scan_type == ScanType.VARIABLE else input_controls
            xvar = getattr(xvar_data, var_to_scan)
            yvar = getattr(output_vars_dict[rho_str], var_to_plot)

            plt.plot(xvar.values, yvar.values)
            plt.xlabel(f'{xvar.label}  {xvar.units_label}')
            plt.ylabel(f'{yvar.label}  {yvar.units_label}')
            plt.title(f'{yvar.name}' + r' ($\rho = {0}$)'.format(rho_str))

            fig.savefig(utils.get_temp_path(f'{profile_type} {sheet_num}.pdf'))
            plt.clf()  # Clear figure data

        merged_pdf = utils.merge_profile_sheets(runid, scan_num, profile_type, is_scan=True)

        # File opening may only work on Windows
        if settings.AUTO_OPEN_PDFS:
            utils.open_file(merged_pdf)

    # Remove figure from memory
    plt.close('all')


def get_scanned_data():
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

    runid = Options.instance.runid
    scan_num = Options.instance.scan_num
    var_to_scan = Options.instance.var_to_scan

    input_vars_dict, output_vars_dict = {}, {}
    rho_values = utils.get_rho_values(runid, scan_num, var_to_scan, SaveType.OUTPUT)

    # Stores InputVariables and OutputVariables data objects for each rho_value
    for rho in rho_values:
        input_vars = InputVariables()
        output_vars = OutputVariables()

        args = (runid, scan_num, var_to_scan, None, rho)
        input_vars.load_data_from_csv(SaveType.INPUT, *args)
        input_vars.load_data_from_csv(SaveType.ADDITIONAL, *args)
        output_vars.load_data_from_csv(SaveType.OUTPUT, *args)

        input_vars_dict[rho] = input_vars
        output_vars_dict[rho] = output_vars

    # Get control_file from rho folder (there's at most one control file, as controls are independent of rho values)
    input_controls = InputControls()
    input_controls.load_from_csv(runid, scan_num, var_to_scan, use_rho=True)

    return input_vars_dict, output_vars_dict, input_controls


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


def main(vars_to_plot, scan_data):
    '''
    Verifies vars_to_plot, then runs the plotting loop for each var_to_plot, runid, and scan_num

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * scan_data (dict): Dictionary of runid to list of scan numbers
    '''

    verify_vars_to_plot(vars_to_plot)

    plotlayout.init()
    plotcolors.init()

    for runid, scan_nums in scan_data.items():
        for scan_num in scan_nums:
            print(f'Initializing data for {runid}, scan {scan_num}...')
            utils.clear_temp_folder()
            Options.instance.load_options(runid, scan_num)
            run_plotting_loop(vars_to_plot)


# Run this file directly to plot scanned variable profiles from previously created scanned data
if __name__ == '__main__':
    scan_data = {}

    '''
    Input Options:
    * vars_to_plot (list): List of output variables to plot

    Examples:
    * vars_to_plot = ['xteMTM', 'xteETGM', 'xteETG', 'gmaMTM', 'omgMTM', 'dbsqprf']
    * vars_to_plot = OutputVariables().get_all_output_vars()
    * vars_to_plot = OutputVariables().get_etgm_vars()
    '''
    vars_to_plot = ['xteETGM']

    '''
    Scan Data:
    * Uncomment the lines you wish to use
        - keys (str): The runid of the scan
        - values (list of int): The scan_numbers to plot from
    '''
    # scan_data['120968A02'] = [1]
    # scan_data['120982A09'] = [1]
    # scan_data['129041A10'] = [1]
    scan_data['TEST'] = [25, 44]

    settings.AUTO_OPEN_PDFS = True

    main(vars_to_plot, scan_data)
