#!/usr/bin/python3

# Standard Packages
import sys; sys.path.insert(0, '../')
import logging
import io

# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from PyQt5.QtGui import QImage
from PyQt5.QtWidgets import QApplication

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

    def on_press(event):
        if event.key == 'x':  # flip x-axis limits
            plt.xlim(plt.xlim()[::-1])
            plt.gcf().canvas.draw()

        if event.key == 'y':  # flip y-axis limits
            plt.ylim(plt.ylim()[::-1])
            plt.gcf().canvas.draw()

        if event.key == "ctrl+c":  # copy figure to clipboard
            with io.BytesIO() as buffer:
                plt.gcf().savefig(buffer)
                QApplication.clipboard().setImage(QImage.fromData(buffer.getvalue()))

    def get_base_variable_or_control(name):
        obj = None
        if hasattr(base_input_vars, name):
            obj = getattr(base_input_vars, name)
        elif hasattr(base_input_controls, name):
            obj = getattr(base_input_controls, name)
        elif hasattr(base_output_vars, name):
            obj = getattr(base_output_vars, name)
        return obj

    levels = 20

    cmap = cm.get_cmap('magma_r', 60)
    colors = np.vstack((np.array([1, 1, 1, 0]), cmap(np.arange(0, cmap.N))[:-10]))
    # colors = np.array(cmap(np.arange(0, cmap.N))[2:-10])
    cmap = LinearSegmentedColormap.from_list('magma_new', colors, N=levels)
    cmap.set_under([1, 1, 1, 0])

    cmap2 = cm.get_cmap('magma_r', 200)
    colors = np.vstack((np.array([1, 1, 1, 0]), cmap2(np.arange(0, cmap2.N))[180:]))
    cmap2 = LinearSegmentedColormap.from_list('magma2_new', colors, N=256)

    var_to_scan = options.var_to_scan
    scan_type = options.scan_type

    input_vars_dict, output_vars_dict, input_controls = datahelper.get_all_rho_data(options)
    base_input_vars, base_output_vars, base_input_controls = datahelper.get_data_objects(options)

    xbase = get_base_variable_or_control(var_to_scan)

    # xvar_data = input_vars_dict[rho_str] if scan_type == ScanType.VARIABLE else input_controls
    # xvar = getattr(xvar_data, var_to_scan)
    # yvar = getattr(output_vars_dict[rho_str], var_to_plot)

    x = np.array(list(output_vars_dict.keys()), dtype=float)
    y = options.scan_range
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)

    rho_strs = input_vars_dict.keys()

    for var_to_plot in vars_to_plot:
        fig = plt.figure()
        ax = plt.gca()

        for s in ax.spines.values():
            s.set_zorder(10)  # Put axes frame on top of everything

        fig.canvas.mpl_connect('key_press_event', on_press)

        # profile_type = f'{var_to_plot}_{var_to_scan}'
        ybase = get_base_variable_or_control(var_to_plot)

        for i, rho_str in enumerate(rho_strs):
            Z[:, i] = getattr(output_vars_dict[rho_str], var_to_plot).values

        # Z[Z>0.02*Z.max()] = 0.02*Z.max()

        # Plot filled contour areas
        plt.contourf(X, Y, Z, cmap=cmap, levels=levels, zorder=3)#, extend='min')

        # Colorbar for filled contour colormap
        cb = plt.colorbar(spacing='uniform', pad=0.02, aspect=30, fraction=0.1, drawedges=True)
        cb.ax.tick_params(size=0, labelsize=6.5)
        # cb.ax.yaxis.set_major_locator(plt.MaxNLocator(8))  # number of tick labels

        # Plot contour lines
        plt.contour(X, Y, Z, cmap=cmap2, levels=levels, linewidths=0.5, alpha=0.5, zorder=4)

        # Base value line
        plt.hlines(1, 0, 1, color="#333", lw=0.5, alpha=0.2, ls='--', zorder=2)

        plt.xlabel(r'$\rho$')
        plt.ylabel(f'{xbase.label} (factors)')

        title = f'{ybase.label} ({ybase.units_label})' if ybase.units_label else f'{ybase.label}'
        # title = fr'{title} [$g_\mathrm{{n_e}} = 0$]'
        plt.title(title)

    plt.show()
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
            options.load(runid, scan_num)
            utils.clear_temp_folder(options)
            if options.var_to_scan:
                run_plotting_loop(vars_to_plot, options)
                utils.clear_temp_folder(options)
            else:
                _log.error(f'\n\tNo variable scan detected for {runid}, scan {scan_num}\n')


# Run this file directly to plot scanned variable profiles from previously created scanned data
if __name__ == '__main__':
    scan_data = {}

    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.SINGLE1,
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
    vars_to_plot = ['gmaETGM']

    '''
    Scan Data:
    * Uncomment the lines you wish to use
        - keys (str): The runid of the scan
        - values (list of int): The scan_numbers to plot from
    '''
    # scan_data['120968A02'] = [3]
    # scan_data['120982A09'] = [1]
    # scan_data['129041A10'] = [1]
    # scan_data['TEST'] = [554]
    # scan_data['138536A01'] = [i for i in range(1, 14)]
    # scan_data['138536A01'] = [i for i in range(97, 109)]
    # scan_data['138536A01'] = [i for i in range(109, 121)]
    scan_data['138536A01'] = [98]

    settings.AUTO_OPEN_PDFS = 1

    main(vars_to_plot, scan_data)
