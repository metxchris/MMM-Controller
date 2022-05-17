#!/usr/bin/python3

"""Creates contour plots of all data from a variable scan

This module is designed to automatically save contour plots and their
associated data from a variable scan.  This code isn't as polished as the
rest of this project, so some values will need to be tweaked depending on
the variables and conditions being plotted.

The following sub-functions within run_plotting_loop may need to be updated as new
variables are added:
* get_vmax_vmin: Get the maximum and minimum values of the plotted variable
* get_ylabel: Get the ylabel of the plotted variable
* get_title: Get the title of the plotted variable
* get_smoothing_sigma: Get the smoothing used for the plotted variable
* get_savename: Get the save name of the plotted variable

Note: At some point the above functions should either be moved to the
Variables objects, or have new objects created to store the associated data.

TODO:
* This module could also be adapted to plot variables directly from a CDF
"""


# Standard Packages
import sys; sys.path.insert(0, '../')
import logging
import io

# 3rd Party Packages
import numpy as np
import scipy.ndimage
import matplotlib.pyplot as plt
from matplotlib import colors, ticker
from matplotlib.ticker import NullFormatter
from PyQt5.QtGui import QImage
from PyQt5.QtWidgets import QApplication

# Local Packages
import modules.options
import modules.utils as utils
import modules.datahelper as datahelper
from modules.variables import InputVariables, OutputVariables
from plotting.modules.plotstyles import PlotStyles, StyleType
import plotting.modules.colormaps


_log = logging.getLogger(__name__)


def run_plotting_loop(vars_to_plot, options, savenameend='', savefig=False, savedata=False):
    """
    Creates PDF Plots of each variable in vars_to_plot

    The x-axis variable is pulled as var_to_scan from the saved options file

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * options (Options): Object containing user options
    * savenameend (str): Name to end file save name with (Optional)
    * savefig (bool): Automatically save the plot if True (Optional)
    * savedata (bool): Automatically save the data if True (Optional)
    """

    def on_press(event):
        if event.key == 'x':  # flip x-axis limits
            plt.xlim(plt.xlim()[::-1])
            plt.gcf().canvas.draw()

        if event.key == 'y':  # flip y-axis limits
            plt.ylim(plt.ylim()[::-1])
            plt.gcf().canvas.draw()

        if event.key == "ctrl+c":  # copy figure to clipboard
            save_format = plt.rcParams['savefig.format']
            plt.rcParams.update({'savefig.format': 'png'})
            with io.BytesIO() as buffer:
                plt.gcf().savefig(buffer)
                QApplication.clipboard().setImage(QImage.fromData(buffer.getvalue()))
                plt.rcParams.update({'savefig.format': save_format})

    def get_base_data(name):
        """Base data is unaltered by any scan factors (multipliers)"""
        obj = None
        if hasattr(output_vars, name):
            obj = getattr(output_vars, name)
        elif hasattr(input_vars, name):
            obj = getattr(input_vars, name)
        elif hasattr(controls, name):
            obj = getattr(controls, name)

        return obj

    def get_var_to_plot_data(name):
        """Get dictionaries of data corresponding to each rho point"""
        obj = None
        if hasattr(output_vars, name):
            obj = output_vars_rho
        elif hasattr(input_vars, name):
            obj = input_vars_rho

        return obj

    def get_vmax_vmin():
        """Get the maximum and minimum values to display"""
        Zmax, Zmin = np.inf, -np.inf
        if 'gma' in var_to_plot:
            Zmax, Zmin = 5e5, -5e5
            Zmax, Zmin = 5e7, -5e7
        elif 'omg' in var_to_plot:
            Zmax, Zmin = 5e5, -5e5
            Zmax, Zmin = 3e6, -3e6
        elif 'xteETGM' == var_to_plot:
            Zmax, Zmin = 2e1, -2e1
            if controls.etgm_sum_modes.values:
                (Zmax, Zmin) = (1e6, -1e6) if var_to_scan != 'time' else (1e6, -1e6)
            else:
                (Zmax, Zmin) = (1e6, -1e6) if var_to_scan != 'time' else (2e1, -2e1)
        elif 'xte2ETGM' == var_to_plot:
            Zmax, Zmin = 1e2, -1e2
            if controls.etgm_sum_modes.values:
                (Zmax, Zmin) = (1e8, -1e8) if var_to_scan != 'time' else (1e8, -1e8)
            else:
                (Zmax, Zmin) = (1e8, -1e8) if var_to_scan != 'time' else (1e1, -1e1)
        elif 'xte' in var_to_plot:
            Zmax, Zmin = 5e2, -5e2
        elif 'xti' in var_to_plot:
            Zmax, Zmin = 5e2, -5e2
        elif 'xdi' in var_to_plot:
            Zmax, Zmin = 5e2, -5e2
        elif 'omegad_gaveETGM' in var_to_plot:
            Zmax, Zmin = 1e7, -1e7
        elif 'gave' in var_to_plot:
            Zmax, Zmin = 5, 0
        elif 'etae' == var_to_plot:
            Zmax, Zmin = 20, -20
        elif 'kyrhosETGM' == var_to_plot:
            Zmax, Zmin = 100, 0
        elif 'kyrhosMTM' == var_to_plot:
            Zmax, Zmin = 5, 0
        elif 'shat' in var_to_plot or 'shear' in var_to_plot:
            Zmax, Zmin = 20, -20
        elif 'omegadETGM' in var_to_plot:
            Zmax, Zmin = 3e6, -3e6
        elif 'omegateETGM' in var_to_plot:
            Zmax, Zmin = 1e7, -1e7
        elif 'omegasETGM' in var_to_plot:
            Zmax, Zmin = 1e7, -1e7
        elif 'omegasetaETGM' in var_to_plot:
            Zmax, Zmin = 1e7, -1e7
        elif 'omegadiffETGM' in var_to_plot:
            Zmax, Zmin = 1e6, -1e6
        elif 'gammadiffETGM' in var_to_plot:
            Zmax, Zmin = 1e6, -1e6
        return Zmax, Zmin

    def get_ylabel():
        """Get the ylabel of the plot"""
        ylabel = xbase.label
        if var_to_scan == 'gne' and options.use_gneabs:
            ylabel = f'$|${ylabel}$|$'

        if 'time' not in var_to_scan and 'kyrho' not in var_to_scan:  # Most yvariables will be plotted as multipliers
            ylabel = f'{ylabel} (multipliers)'
        elif xbase.units_label:
            ylabel = f'{ylabel} ({xbase.units_label})'
        return ylabel

    def get_title():
        """Get the plot title"""
        title = f'{ybase.label} ({ybase.units_label})' if ybase.units_label else f'{ybase.label}'
        if options.use_gnezero:
            title = fr'{title} [$g_\mathrm{{ne}} = 0$]'
        if options.use_gtezero:
            title = fr'{title} [$g_\mathrm{{Te}} = 0$]'
        if options.use_gneabs:
            title = fr'{title} with $|g_\mathrm{{ne}}|$'
        if controls.etgm_exbs.values and var_to_plot != 'wexbs':
            title = fr'{title} $[\omega_{{E \times\! B}}\,\,\mathrm{{on}}]$'
        if controls.etgm_sum_modes.values and 'ETGM' in var_to_plot and '\chi' in ybase.label:
            title = fr'$_{{^\sum}}${title}'
        return title

    def get_smoothing_sigma():
        """Get the sigma used for Gaussian signal smoothing"""
        medium_smoothing = ['xteETGM', 'xte2ETGM', 'xdiETGM']
        light_smoothing = [
            'omgETGM', 'kyrhoeETGM', 'kyrhosETGM',
            'omegateETGM', 'walfvenunit', 'omegadETGM',
            'omegasETGM', 'omegasetaETGM', 'omegadiffETGM',
            'gammadiffETGM'
        ]

        sigma = (0, 0)
        if options.var_to_scan == 'time':
            sigma = (0, 0)
        elif var_to_plot in medium_smoothing:
            sigma = (2, 2)
        elif var_to_plot in light_smoothing:
            sigma = (1, 1)
        return sigma

    def get_contour_levels():
        """Get the displayed contour levels"""
        contour_level_count = 20  # Average number of contours to show (actual number will vary)
        loc = ticker.MaxNLocator(contour_level_count + 1, min_n_ticks=contour_level_count - 4)
        lvls = loc.tick_values(Z.min(), Z.max())

        ## Code for logarithmic scales
        # loc = ticker.LogLocator(base=10)
        # lvls = loc.tick_values(max(np.power(10, np.floor(np.log10(Z.min()) - 1)), 1), np.power(10, np.ceil(np.log10(Z.max()) + 1)))
        # args_both['locator'] = ticker.LogLocator(base=10)
        # locator=ticker.LogLocator()

        ## Condition for nonlinear contour levels
        # lowest_bucket = Z[Z <= lvls[1]].size / Z.size
        # if (var_to_plot == options.var_to_scan or var_to_plot == 'shat_gxi') and lowest_bucket > 0.4:
        #     power = 1.2
        #     lvls = np.round(lvls**power / (lvls**power)[-1] * lvls[-1], 3)

        return lvls

    def get_savename():
        """Get the name to save the file as"""
        savename = adjustment_name
        if options.use_gnezero:
            savename = f'{savename}_gne0'
        if options.use_gtezero:
            savename = f'{savename}_gte0'
        if options.use_gneabs:
            savename = f'{savename}_gneabs'
        if controls.etgm_exbs.values:
            savename = f'{savename}_exbs1'
        if controls.etgm_sum_modes.values and 'xte' in var_to_plot:
            savename = f'{savename}_sum'
        if savenameend:
            savename = f'{savename}_{savenameend}'
        return f'{savename}_{var_to_plot}'

    # Plotting loop initialization
    colormaps = plotting.modules.colormaps.get_colormaps()
    var_to_scan = options.var_to_scan
    adjustment_name = options.adjustment_name or options.var_to_scan

    input_vars_rho, output_vars_rho, controls_rho = datahelper.get_all_rho_data(options)
    input_vars, output_vars, controls = datahelper.get_data_objects(options)

    rho_strs = input_vars_rho.keys()
    xbase = get_base_data(var_to_scan)

    x = np.array(list(output_vars_rho.keys()), dtype=float)
    y = options.scan_range
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)

    for var_to_plot in vars_to_plot:

        fig, ax = plt.gcf(), plt.gca()  # update figure variables for current plot

        if var_to_plot == 'var_to_scan':
            var_to_plot = var_to_scan

        if var_to_plot == 'time':
            _log.warning(f'\n\tNothing to plot when {var_to_plot} is time, continuing to next variable...')
            continue  # nothing to plot (this can happen when using 'var_to_scan')

        ybase = get_base_data(var_to_plot)
        ydata = get_var_to_plot_data(var_to_plot)

        if not isinstance(ybase.values, np.ndarray):
            # Contours need to be np.ndarray. This can happen when using the
            # same variable plotting list for many different scan types, so
            # no exception is raised
            _log.warning(f'\n\t{var_to_plot} is not an np.ndarray, continuing to next variable...')
            continue

        print(f'- {options.scan_num}, {var_to_plot}')

        for i, rho_str in enumerate(rho_strs):
            Z[:, i] = getattr(ydata[rho_str], var_to_plot).values

        if np.isnan(Z).any():
            # This can happen due to calculation errors or or when loading
            # data using old variable definitions, so a warning is thrown but
            # no exception is raised
            _log.warning(f'\n\tnan values found in {var_to_plot}, continuing to next variable...')
            continue

        if not savefig:  # Connect key-press handler when not autosaving figures
            fig.canvas.mpl_connect('key_press_event', on_press)

        for s in ax.spines.values():
            s.set_zorder(10)  # Put axes frame on top of everything

        # Apply smoothing and restore original extrema
        Zmax, Zmin = Z.max(), Z.min()
        Z = scipy.ndimage.gaussian_filter(Z, sigma=get_smoothing_sigma())
        Z = np.minimum(np.maximum(Z, Zmin), Zmax)

        # Default arguments for filled contours, line contours, and both types
        args_fill = {'zorder': -3}
        args_line = {'zorder': -2, 'linewidths': 0.4}
        args_both = {}

        # Max and min values for Z
        Zmax, Zmin = get_vmax_vmin()

        # Clamp Z between Zmin and Zmax
        Z = np.minimum(np.maximum(Z, Zmin), Zmax)

        # Set colormaps
        if Z.min() >= 0:
            args_both['norm'] = None
            args_fill['cmap'] = colormaps['magma_positive']
            args_line['cmap'] = colormaps['magma_positive_lines']
        elif Z.max() <= 0:
            args_both['norm'] = None
            args_fill['cmap'] = colormaps['magma_negative']
            args_line['cmap'] = colormaps['magma_negative_lines']
        else:
            args_both['norm'] = colors.CenteredNorm()
            args_fill['cmap'] = colormaps['magma_both']
            args_line['cmap'] = colormaps['magma_both_lines']

        # Set colorbar extensions
        if Z.max() >= Zmax and Z.min() <= Zmin:
            args_both['extend'] = 'both'
        elif Z.max() >= Zmax:
            args_both['extend'] = 'max'
        elif Z.min() <= Zmin:
            args_both['extend'] = 'min'
        else:
            args_both['extend'] = 'neither'

        # Plot contour fills
        cf = plt.contourf(X, Y, Z, levels=get_contour_levels(), **args_fill, **args_both)

        # Remove lowest level when it is also a boundary, so an extra contour line isn't drawn
        if cf.levels[0] == Z.min():
            if args_both['extend'] != 'min' and args_both['extend'] != 'both':
                cf.levels = cf.levels[1:]

        # Colorbar for filled contour colormap
        cb = plt.colorbar(spacing='proportional', pad=0.02, aspect=30, fraction=0.1, drawedges=True)
        cb.ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        cb.ax.tick_params(size=0, labelsize=plt.rcParams['ytick.labelsize'] - 0.5)

        # An upper colorbar extension moves the offset text down, so we raise it back up
        if 'extend' in args_both:
            if args_both['extend'] == 'max' or args_both['extend'] == 'both':
                cb.ax.yaxis.OFFSETTEXTPAD += 5

        # Plot contour lines
        cl = plt.contour(X, Y, Z, levels=cf.levels, **args_line, **args_both)

        # Override the linestyles based on the levels.
        for line, lvl in zip(cl.collections, cl.levels):
            if lvl < 0:
                line.set_linestyle('--')
            elif lvl == 0:
                line.set_linestyle('-')
            else:
                line.set_linestyle('-')

        ax.yaxis.set_minor_formatter(NullFormatter())

        plt.xlabel(r'$\rho$')
        plt.ylabel(get_ylabel())
        plt.title(get_title())

        if savedata or savefig:
            savedir_base = f'{utils.get_plotting_contours_path()}\\{options.runid}'
            utils.create_directory(savedir_base)

        if savedata:
            savedir = f'{savedir_base}\\data'
            utils.create_directory(savedir)
            _save_to_csv(X, Y, Z, f'{savedir}\\{get_savename()}')

        if savefig:
            savedir = f'{savedir_base}\\figures'
            utils.create_directory(savedir)
            fig.savefig(f'{savedir}\\{get_savename()}')
        else:
            plt.show()

        fig.clear()

    plt.close(fig)


def _save_to_csv(X, Y, Z, savename):
    """
    Save plotted data to a CSV

    Data is saved into three CSV's, one for each dimension X, Y, and Z.  When
    using normalized data for X and Y, the CSV's could be simplified and
    saved as single arrays.  However, the full matrices are being saved
    anyways if generalized values for X and Y are ever used.

    Raises:
    * FileNotFoundError: If the file cannot be found after saving it
    """

    prec, col_pad = 4, 6
    col_len = prec + col_pad
    fmt_str = f'%{col_len}.{prec}e'

    vars_to_save = ['X', 'Y', 'Z']

    # Crop X, Y data if values are uniform per row or column, respectively
    if (X.max(0) == X.min(0)).all():
        X = X[0, :]
    if (Y.max(1) == Y.min(1)).all():
        Y = Y[:, 0]

    for var in vars_to_save:
        var_savename = f'{savename}_{var}.csv'
        var_data = eval(var)
        np.savetxt(var_savename, var_data, fmt=fmt_str, delimiter=',')
        if not utils.check_exists(var_savename):
            raise FileNotFoundError(f'Failed to save {var_savename}\n'
                                    '\tMake sure Python has file writing permissions')


def _verify_vars_to_plot(vars_to_plot):
    '''
    Variables are verified to be members of OutputVariables before the plotting loop begins

    Parameters:
    * vars_to_plot (list):  List of output variables to plot

    Raises:
    * NameError: If the variable to plot is not found in OutputVariables
    '''

    output_vars = OutputVariables()
    input_vars = InputVariables()
    for var_to_plot in vars_to_plot:
        if var_to_plot == 'var_to_scan':
            continue
        if not hasattr(output_vars, var_to_plot) and not hasattr(input_vars, var_to_plot):
            raise NameError(f'{var_to_plot} not found in Variables classes')


def main(vars_to_plot, scan_data, savenameend='', savefig=False, savedata=False):
    '''
    Verifies vars_to_plot, then runs the plotting loop for each var_to_plot, runid, and scan_num

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * scan_data (dict): Dictionary of runid to list of scan numbers
    * savenameend (str): Name to end file save name with
    * savefig (bool): Automatically save the plot if True (Optional)
    * savedata (bool): Automatically save the data if True (Optional)
    '''

    utils.init_logging()
    _verify_vars_to_plot(vars_to_plot)
    options = modules.options.Options()

    if savefig:
        print(f'Files will be saved in:\n\t{utils.get_plotting_contours_path()}')

    plt.figure()  # Instantiate the figure

    for runid, scan_nums in scan_data.items():
        for scan_num in scan_nums:
            options.load(runid, scan_num)
            if options.var_to_scan:
                print(f'\nInitializing data for {runid}, scan {scan_num}, {options.var_to_scan}...')
                run_plotting_loop(vars_to_plot, options, savenameend, savefig, savedata)
            else:
                _log.error(f'\n\tNo variable scan detected for {runid}, scan {scan_num}\n')


# Run this file directly to plot scanned variable profiles from previously created scanned data
if __name__ == '__main__':
    scan_data = {}

    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.SINGLE1B,
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
    })

    # Text to append to the end of the generated save name (Optional)
    savenameend = ''

    """
    Variables to Plot (List):
    * Uncomment the line you wish to use
    """

    # vars_to_plot = ['var_to_scan']
    vars_to_plot = ['gmanormMTM', 'omgnormMTM', 'xteMTM']
    # vars_to_plot = ['gmaETGM', 'omgETGM', 'xteETGM', 'xte2ETGM']
    # vars_to_plot = OutputVariables().get_all_output_vars()
    # vars_to_plot = OutputVariables().get_etgm_vars()
    # vars_to_plot = OutputVariables().get_mtm_vars()

    """
    Scan Data:
    * Uncomment the line you wish to use
        - keys (str): The runid of the scan
        - values (list of int): The scan_numbers to plot from
    """

    # scan_data['120968A02'] = [12]
    scan_data['129016A04'] = [16]
    # scan_data['138536A01'] = [i for i in range(1716, 1738 + 1)]
    # scan_data['138536A01'] = [i for i in [*range(1716, 1738 + 1), *range(1756, 1763 + 1)]]

    """
    Plotting Options:
    * savefig: show figure when True, autosave figure without showing it when False
    * savedata: autosave data into CSVs when True
    """
    main(vars_to_plot, scan_data, savenameend=savenameend, savefig=False, savedata=False)
