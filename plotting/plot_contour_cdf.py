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


class PlotOptions:
    """ Store options for the contour plot

    Initialization parameters determine plot settings and must be specified as
    keyword arguments.
    """

    def __init__(self, **kwargs):
        self.ymin: float | None = None
        self.ymax: float | None = None
        self.xmin: float | None = None
        self.xmax: float | None = None
        self.raw: bool | None = None

        self._set_kwargs(kwargs)

    def _set_kwargs(self, kwargs):
        """Set member values using keyword arguments"""
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                _log.error(f'\n\t"{key}" is not a valid parameter for PlotOptions')


def run_plotting_loop(vars_to_plot, options, plot_options=None, savenameend='', savefig=False, savedata=False):
    """
    Creates PDF Plots of each variable in vars_to_plot

    The x-axis variable is pulled as var_to_scan from the saved options file

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * options (Options): Object containing user options
    * plot_options (PlotOptions): Object containing options for the contour plot (Optional)
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
        # obj = None
        # if hasattr(output_vars, name):
        #     obj = getattr(output_vars, name)
        # elif hasattr(mmm_vars, name):
        #     obj = getattr(mmm_vars, name)
        # elif hasattr(controls, name):
        #     obj = getattr(controls, name)

        # return obj
        ...

    def get_var_to_plot_data(name):
        """Get dictionaries of data corresponding to each rho point"""
        obj = None
        if hasattr(output_vars, name):
            obj = output_vars_rho
        elif hasattr(mmm_vars, name):
            obj = mmm_vars_rho

        return obj

    def get_ylabel():
        """Get the ylabel of the plot"""
        ylabel = f'{ybase.label} ({ybase.units_label})' if ybase.units_label else f'{ybase.label}'
        if var_to_scan == 'gne' and options.use_gneabs:
            ylabel = f'$|${ylabel}$|$'

        if 'time' not in var_to_scan and 'kyrho' not in var_to_scan:  # Most yvariables will be plotted as multipliers
            ylabel = f'{ylabel}'
        elif xbase.units_label:
            ylabel = f'{ylabel} ({xbase.units_label})'
        return ylabel

    def get_title():
        """Get the plot title"""
        title = f'{zbase.label} ({zbase.units_label})' if zbase.units_label else f'{zbase.label}'
        if options.use_gnezero:
            title = fr'{title} [$g_\mathrm{{ne}} = 0$]'
        if options.use_gtezero:
            title = fr'{title} [$g_\mathrm{{Te}} = 0$]'
        if options.use_gneabs:
            title = fr'{title} with $|g_\mathrm{{ne}}|$'
        # if controls.etgm_exbs.values and var_to_plot != 'wexb':
        #     title = fr'{title} $[\omega_{{E \times\! B}}\,\,\mathrm{{on}}]$'
        # if controls.etgm_sum_modes.values and 'ETGM' in var_to_plot and '\chi' in zbase.label:
        #     title = fr'$_{{^\sum}}${title}'
        return title

    def get_smoothing_sigma():
        """Get the sigma used for Gaussian signal smoothing"""
        # medium_smoothing = ['xteETGM', 'xte2ETGM', 'xdiETGM']
        # light_smoothing = [
        #     'omgETGM', 'kyrhoeETGM', 'kyrhosETGM',
        #     'omegateETGM', 'walfvenunit', 'omegadETGM',
        #     'omegasETGM', 'omegasetaETGM', 'omegadiffETGM',
        #     'gammadiffETGM'
        # ]

        # sigma = (0, 0)
        # if options.var_to_scan == 'time':
        #     sigma = (0, 0)
        # elif var_to_plot in medium_smoothing:
        #     sigma = (2, 2)
        # elif var_to_plot in light_smoothing:
        #     sigma = (1, 1)
        # return sigma
        ...

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
    
    adjustment_name = options.adjustment_name or options.var_to_scan

    # mmm_vars_rho, output_vars_rho, controls_rho = datahelper.get_all_rho_data(options)

    if not plot_options.raw:
        mmm_vars, cdf_vars, __ = datahelper.initialize_variables(options)
    else:
        __, cdf_vars, mmm_vars = datahelper.initialize_variables(options)

    # rho_strs = mmm_vars_rho.keys()
    
    if not plot_options.raw:
        x = mmm_vars.rho.values[:, 0]
    else:
        x = mmm_vars.xb.values[:, 0]
    y = mmm_vars.time.values

    x0 = np.zeros_like(x)
    y0 = np.zeros_like(y)

    x0[:] = x[:]
    y0[:] = y[:]

    # Apply boundary limits to x, y variables
    if plot_options is not None:
        if plot_options.xmin is not None:
            x = x[x0 >= plot_options.xmin]
            # y = y[x0 >= plot_options.xmin]
        if plot_options.xmax is not None:
            x = x[x0 <= plot_options.xmax]
            # y = y[x0 <= plot_options.xmax]
        if plot_options.ymin is not None:
            y = y[y0 >= plot_options.ymin]
            # x = x[y0 >= plot_options.ymin]
        if plot_options.ymax is not None:
            y = y[y0 <= plot_options.ymax]
            # x = x[y0 <= plot_options.ymax]

    X, Y = np.meshgrid(x, y)
    X0, Y0 = np.meshgrid(x0, y0)
    Z = np.zeros_like(X)

    for var_to_plot in vars_to_plot:

        var_to_scan = var_to_plot

        if not plot_options.raw:
            xbase = mmm_vars.rho
        else:
            xbase = mmm_vars.xb

        ybase = mmm_vars.time
        zbase = getattr(mmm_vars, var_to_plot)

        fig, ax = plt.gcf(), plt.gca()  # update figure variables for current plot

        if var_to_plot == 'var_to_scan':
            var_to_plot = var_to_scan

        if var_to_plot == 'time':
            _log.warning(f'\n\tNothing to plot when {var_to_plot} is time, continuing to next variable...')
            continue  # nothing to plot (this can happen when using 'var_to_scan')

        # ybase = get_base_data(var_to_plot)
        ydata = mmm_vars

        # if not isinstance(ybase.values, np.ndarray):
        #     # Contours need to be np.ndarray. This can happen when using the
        #     # same variable plotting list for many different scan types, so
        #     # no exception is raised
        #     _log.warning(f'\n\t{var_to_plot} is not an np.ndarray, continuing to next variable...')
        #     continue

        print(f'- {var_to_plot}')

        # for i, rho_str in enumerate(rho_strs):
        Z = getattr(ydata, var_to_plot).values.T

        if np.isnan(Z).any():
            # This can happen due to calculation errors or or when loading
            # data using old variable definitions, so a warning is thrown but
            # no exception is raised
            _log.warning(f'\n\tnan values found in {var_to_plot}, continuing to next variable...')
            continue

        if plot_options is not None:
            zshape = X.shape
            if plot_options.xmin is not None:
                Z = Z[X0 >= plot_options.xmin]
            if plot_options.xmax is not None:
                Z = Z[X0 <= plot_options.xmax]
            if plot_options.ymin is not None:
                Z = Z[Y0 >= plot_options.ymin]
            if plot_options.ymax is not None:
                Z = Z[Y0 <= plot_options.ymax]
            Z = Z.reshape(zshape)

        if not savefig:  # Connect key-press handler when not autosaving figures
            fig.canvas.mpl_connect('key_press_event', on_press)

        for s in ax.spines.values():
            s.set_zorder(10)  # Put axes frame on top of everything

        # Apply smoothing and restore original extrema
        Zmax, Zmin = Z.max(), Z.min()
        # Z = scipy.ndimage.gaussian_filter(Z, sigma=get_smoothing_sigma())
        Z = np.minimum(np.maximum(Z, Zmin), Zmax)

        # Default arguments for filled contours, line contours, and both types
        args_fill = {'zorder': -3}
        args_line = {'zorder': -2, 'linewidths': 0.4}
        args_both = {}

        # Max and min values for Z
        Zmax = zbase.contour_max
        Zmin = zbase.contour_min

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

        plt.xlabel(r'$\hat{\rho}$')
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


def main(vars_to_plot, scan_data, plot_options=None, savenameend='', savefig=False, savedata=False):
    '''
    Verifies vars_to_plot, then runs the plotting loop for each var_to_plot, runid, and scan_num

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * scan_data (dict): Dictionary of runid to list of scan numbers
    * plot_options (PlotOptions): Object containing options for the contour plot (Optional)
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

    for runid in scan_data:
        options.set(runid=runid)
        print(f'\nInitializing data for {runid}...')
        run_plotting_loop(vars_to_plot, options, plot_options, savenameend, savefig, savedata)



# Run this file directly to plot scanned variable profiles from previously created scanned data
if __name__ == '__main__':
    scan_data = []

    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.AIP2,
    )

    plot_options = PlotOptions(
        # xmin=0.0,
        # xmax=0.8,
        # ymin=1.5,
        # ymax=0.17,
        raw=1,
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
    })

    """
    Variables to Plot (List):
    * Uncomment the line you wish to use
    """

    # vars_to_plot = ['var_to_scan']
    vars_to_plot = ['ti', 'te', 'ne', 'ni']
    # vars_to_plot = ['gne', 'gte',]
    # vars_to_plot = ['gne', 'gte',]
    # vars_to_plot = ['gmaW20ii', 'gmaW20ee', 'gmaW20ie']
    # vars_to_plot = ['gmaW20i', 'gmaW20e']
    vars_to_plot = ['xkemmm', 'xkiw20', 'xdew20', 'xkew20', 'xkemtm']
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

    # scan_data.append('138536A01')
    # scan_data.append('129041A10')
    # scan_data.append('129016T17')

    # vars_to_plot = ['walltime', ]
    # vars_to_plot = ['wexbsv2', ]
    # vars_to_plot = ['xkiw20', 'xkdw20', 'xkew20', ]
    # vars_to_plot = ['xkemmm07', 'condepr', 'conde']
    # vars_to_plot = ['xkimmm07', 'condipr', 'condi']
    # vars_to_plot = ['xkdmmm07', 'condipr', 'condi']
    # vars_to_plot = ['te', 'ti',]
    
    vars_to_plot = ['xkimmm', 'xkdmmm', 'xkemmm', 'xkidrbm', 'xkhdrbm', 'xkedrbm', 'gammadbm', 'omegadbm', 'kyrsdbm']
    # scan_data.append('129016Q39') # DRIBM Crash
    # scan_data.append('129016Q71') # DRIBM success
    # scan_data.append('129016Q51') # All Disabled
    # scan_data.append('129016W47') # MMM 9.0.6 DBM Failure, shat_e
    

    # vars_to_plot = ['te',]
    # scan_data.append('129016A03')
    # scan_data.append('129016Q57')

    # vars_to_plot = ['xkemmm', 'xkeetg', 'kyrsetg', 'gammaetg', 'omegaetg']
    # scan_data.append('129016Q58') # ETGM 
    # scan_data.append('129016Q62') # ETGM with overly large diffusivity

    # vars_to_plot = ['xkemtm', 'gammamtm', 'omegamtm', 'kyrsmtm', 'xdbmtm']
    # scan_data.append('129016Q79') # Custom MMM library 

    # vars_to_plot = ['xkiw20', 'xdew20', 'xkew20', 'xkzw20', 'xppw20', 'xptw20', ]
    vars_to_plot = ['xkiw20', 'xkew20','xkemtm', ]
    vars_to_plot = ['gti', 'gte','ti','te', ]
    # vars_to_plot = ['wexbsmod', 'wexbsv2', 'te', 'ti', 'wexb',]
    # vars_to_plot = ['vphimmm', 'gamma1w20', 'gamma2w20', 'omega1w20', 'omega2w20',]
    vars_to_plot = ['xkemmm', 'xkimmm', 'xkdmmm', 'xkzmmm',  'xppmmm',  'xptmmm', ]
    # vars_to_plot = ['vphimmm', 'kyrsdbm', 'kyrsetg', 'kyrsmtm',  'xkeetg',  'xptmmm', ]
    vars_to_plot = [ 'xkemtm', 'xkeetg']
    # vars_to_plot = [ 'gamma1w20', 'xkeetg']
    # scan_data.append('129016Q68') # MMM v8, W20 only
    # scan_data.append('129016Q73') # MMM v8, all models, cal = 1E-6
    # scan_data.append('129016Q57') # W20 only, max xte = 0.01
    # scan_data.append('129016Q92') # MMM9.0.6, ETGM and MTM only
    # scan_data.append('129016Q93') # MMM9.0.6, W20
    # scan_data.append('129016Q94') # MMM9.0.6, W20 with +- diffusivity
    # scan_data.append('129016Q99') # MMM9.0.7 tshare
    # scan_data.append('129016Q50') # MMM disabled
    # scan_data.append('129016W46') # w20, etgm, mtm
    # scan_data.append('129016W53') # MMM 9.0.6, (W20: te, ti, pphi)
    # scan_data.append('129016W56') # MMM 9.0.6, (W20, 0 pinch: te, ti, pphi)
    # scan_data.append('129016W55') # MMM 9.0.6, (W20: te)

    # scan_data.append('129016Z11') # MMM 9.0.7  W20, ETGM, MTM crash; +- Chi
    # scan_data.append('129016Z13') # MMM 9-8  W20, MTM crash; + Chi
    # scan_data.append('129016Z29') # MMM 8.2.1
    # scan_data.append('129016Z33') # MMM 9.0.7
    # scan_data.append('129016Z36') # MMM 9.0.7 + ETGM

    # scan_data.append('120982W30') # 8.2.1
    # scan_data.append('120982W31') # MMM 9.0.7 + ETGM
    # scan_data.append('120982W32') # MMM 9.0.7
    scan_data.append('129017W02') # MMM 8.2.1
    # scan_data.append('129017W03') # MMM 9.0.7
    

    # vars_to_plot = ['xkidrbm', 'xkddrbm', 'xkedrbm', 'xkhdrbm']
    # scan_data.append('129016Q71') # DBM, 0.001 cal

    """
    Plotting Options:
    * savefig: show figure when True, autosave figure without showing it when False
    * savedata: autosave data into CSVs when True
    """
    main(vars_to_plot, scan_data, plot_options, savenameend='NEW', savefig=0, savedata=False)
