# Standard Packages
import sys; sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt

# Local Packages
import settings
import modules.options as options
import modules.utils as utils
import modules.constants as constants
from modules.enums import ScanType, MergeType
from modules.controls import InputControls
from modules.variables import OutputVariables
from plotting.modules.styles import singlescan as plotlayout
from plotting.modules.colors import mmmscan as plotcolors


def run_plotting_loop(vars_to_plot):
    '''
    Creates PDF Plots of each variable in vars_to_plot

    The x-axis variable is pulled as var_to_scan from the saved options file

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    '''

    fig = plt.figure()
    runid = options.instance.runid
    scan_num = options.instance.scan_num
    var_to_scan = options.instance.var_to_scan
    scan_type = options.instance.scan_type

    input_vars_dict, output_vars_dict, input_controls = utils.get_all_rho_data(runid, scan_num, var_to_scan)
    base_input_vars, base_output_vars, base_input_controls = utils.get_base_data(runid, scan_num)

    xbase = None
    if hasattr(base_input_vars, var_to_scan):
        xbase = getattr(base_input_vars, var_to_scan)
    elif hasattr(base_input_controls, var_to_scan):
        xbase = getattr(base_input_controls, var_to_scan)

    for var_to_plot in vars_to_plot:
        print(f'Creating scanned variable PDF for {var_to_plot} vs {var_to_scan}...')

        # rho_strs are the same for both input and output variables
        rho_strs = input_vars_dict.keys()
        profile_type = f'{var_to_plot}_{var_to_scan}'
        ybase = getattr(base_output_vars, var_to_plot)

        for i, rho_str in enumerate(rho_strs):
            sheet_num = constants.SHEET_NUM_FMT_STR.format(i)

            # Plot scanned values
            xvar_data = input_vars_dict[rho_str] if scan_type == ScanType.VARIABLE else input_controls
            xvar = getattr(xvar_data, var_to_scan)
            yvar = getattr(output_vars_dict[rho_str], var_to_plot)
            plt.plot(xvar.values, yvar.values)

            # Plot base value
            xbase_values = xbase.values[i] if type(xbase.values) is np.ndarray else xbase.values
            plt.plot(xbase_values, ybase.values[i], 'o', markeredgewidth=1.25, markersize=3, alpha=0.8)

            plt.xlabel(f'{xvar.label}  {xvar.units_label}')
            plt.ylabel(f'{yvar.label}  {yvar.units_label}')
            plt.title(f'{yvar.name}'r' ($\rho = {0}$)'.format(rho_str))
            fig.savefig(utils.get_temp_path(f'{profile_type} {sheet_num}.pdf'))
            fig.clear()

        merged_pdf = utils.merge_profile_sheets(runid, scan_num, profile_type, MergeType.RHOVALUES, var_to_scan)

        # File opening may only work on Windows
        if settings.AUTO_OPEN_PDFS:
            utils.open_file(merged_pdf)

    # Remove figure from memory
    plt.close(fig)


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
            options.instance.load_options(runid, scan_num)
            if options.instance.var_to_scan:
                run_plotting_loop(vars_to_plot)
            else:
                print(f'ERROR: No variable scan detected')


# Run this file directly to plot scanned variable profiles from previously created scanned data
if __name__ == '__main__':
    scan_data = {}

    '''
    Input options:
    * vars_to_plot (list): List of output variables to plot

    Examples:
    * vars_to_plot = ['xteETGM']
    * vars_to_plot = ['xteMTM', 'xteETGM', 'xteETG', 'gmaMTM', 'omgMTM', 'dbsqprf']
    * vars_to_plot = OutputVariables().get_all_output_vars()
    * vars_to_plot = OutputVariables().get_etgm_vars()
    '''
    vars_to_plot = OutputVariables().get_etgm_vars()

    '''
    Scan Data:
    * Uncomment the lines you wish to use
        - keys (str): The runid of the scan
        - values (list of int): The scan_numbers to plot from
    '''
    # scan_data['120968A02'] = [1]
    # scan_data['120982A09'] = [1]
    # scan_data['129041A10'] = [1]
    # scan_data['TEST'] = [181]
    # scan_data['138536A01'] = [i for i in range(100, 126)]
    scan_data['138536A01'] = [124]

    settings.AUTO_OPEN_PDFS = 1

    main(vars_to_plot, scan_data)
