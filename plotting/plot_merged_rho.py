# Standard Packages
import sys; sys.path.insert(0, '../')
import logging

# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt

# Local Packages
import settings
import modules.options
import modules.utils as utils
import modules.datahelper as datahelper
import modules.constants as constants
from modules.enums import ScanType, MergeType
from modules.variables import OutputVariables
from plotting.modules.plotstyles import PlotStyles, StyleType


_log = logging.getLogger(__name__)


def run_plotting_loop(vars_to_plot, options):
    '''
    Creates PDF Plots of each variable in vars_to_plot

    The x-axis variable is pulled as var_to_scan from the saved options file

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * options (Options): Object containing user options
    '''

    fig = plt.figure()
    var_to_scan = options.var_to_scan
    scan_type = options.scan_type

    input_vars_dict, output_vars_dict, input_controls = datahelper.get_all_rho_data(options)
    base_input_vars, base_output_vars, base_input_controls = datahelper.get_data_objects(options)

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
            sheet_num = f'{i:{constants.SHEET_NUM_FMT}}'

            xbase_values = xbase.values[i] if type(xbase.values) is np.ndarray else xbase.values
            xvar_data = input_vars_dict[rho_str] if scan_type == ScanType.VARIABLE else input_controls
            xvar = getattr(xvar_data, var_to_scan)
            yvar = getattr(output_vars_dict[rho_str], var_to_plot)

            if xbase_values < 0:
                plt.plot([], [])  # Advance the cycler twice
                plt.plot([], [])

            plt.plot(xvar.values, yvar.values, dashes=[1, 0])
            plt.plot(xbase_values, ybase.values[i])

            if xbase_values < 0:
                plt.xlim(plt.xlim()[::-1])

            plt.xlabel(f'{xvar.label}  {xvar.units_label}')
            plt.ylabel(f'{yvar.label}  {yvar.units_label}')
            plt.title(f'{yvar.name}'r' ($\rho = {0}$)'.format(rho_str))
            fig.savefig(utils.get_temp_path(f'{profile_type} {sheet_num}.pdf'))
            fig.clear()

        merged_pdf = utils.merge_profile_sheets(options, profile_type, MergeType.RHOVALUES)

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

    Raises:
    * NameError: If the variable to plot is not found in OutputVariables
    '''

    output_vars = OutputVariables()
    for var_to_plot in vars_to_plot:
        if not hasattr(output_vars, var_to_plot):
            raise NameError(f'{var_to_plot} not found in OutputVariables class')


def main(vars_to_plot, scan_data):
    '''
    Verifies vars_to_plot, then runs the plotting loop for each var_to_plot, runid, and scan_num

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * scan_data (dict): Dictionary of runid to list of scan numbers
    '''

    utils.init_logging()
    verify_vars_to_plot(vars_to_plot)
    options = modules.options.Options()

    for runid, scan_nums in scan_data.items():
        for scan_num in scan_nums:
            print(f'Initializing data for {runid}, scan {scan_num}...')
            utils.clear_temp_folder()
            options.load(runid, scan_num)
            if options.var_to_scan:
                run_plotting_loop(vars_to_plot, options)
            else:
                _log.error(f'\n\tNo variable scan detected for {runid}, scan {scan_num}\n')


# Run this file directly to plot scanned variable profiles from previously created scanned data
if __name__ == '__main__':
    scan_data = {}

    PlotStyles(
        axes=StyleType.Axes.GRAY,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.SINGLE3,
    )

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
    scan_data['138536A01'] = [3]

    settings.AUTO_OPEN_PDFS = 1

    main(vars_to_plot, scan_data)
