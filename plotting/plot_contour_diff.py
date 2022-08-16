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


class ContourPlotData:
    def __init__(self, runid, scannum, var_to_plot):
        self.runid: str = runid
        self.scannum: int = scannum
        self.var_to_plot: str = var_to_plot
        self.options = modules.options.Options().load(runid, scannum)


class PlotOptions:
    def __init__(self, title='', savenameend='', savefig=False, savedata=False, difftype='diff'):
        self.title: str = title
        self.savefig: bool = savefig
        self.savedata: bool = savedata
        self.savenameend: str = savenameend

        difftypes = ['diff', 'absdiff', 'absratio', 'ratio']
        if difftype not in difftypes:
            raise ValueError(
                f'difftype \'{difftype}\' not in list of accepted types:'
                f'\n\t{difftypes}'
            )

        self.difftype: str = difftype


def plot_contour_difference(contour_list, plot_options):
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

    def get_base_data(name, output_vars, input_vars, controls):
        """Base data is unaltered by any scan factors (multipliers)"""
        obj = None
        if hasattr(output_vars, name):
            obj = getattr(output_vars, name)
        elif hasattr(input_vars, name):
            obj = getattr(input_vars, name)
        elif hasattr(controls, name):
            obj = getattr(controls, name)

        return obj

    def get_var_to_plot_data(name, output_vars_rho, input_vars_rho):
        """Get dictionaries of data corresponding to each rho point"""
        obj = None
        if hasattr(output_vars, name):
            obj = output_vars_rho
        elif hasattr(input_vars, name):
            obj = input_vars_rho

        return obj

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

    def get_title(plot_options):
        """Get the plot title"""
        if plot_options.title:
            return plot_options.title

        label = ''
        label1 = ybase1.label
        label2 = ybase2.label

        if label1 == label2:  # Same variables
            label = f'{label1}:'

            if plot_options.difftype == 'diff':
                return f'{label} New - Old'
            elif plot_options.difftype == 'absdiff':
                return fr'{label} $|$New - Old$|$'
            elif plot_options.difftype == 'ratio':
                return f'{label} (New - Old) / (New + Old)'
            elif plot_options.difftype == 'absratio':
                return fr'{label} $|$New - Old$|$ / $|$New + Old$|$'

        else:  # Not same variables
            if plot_options.difftype == 'diff':
                return f'{ybase2.label} - {ybase1.label}'
            elif plot_options.difftype == 'absdiff':
                return fr'$|${ybase2.label} - {ybase1.label}$|$'
            elif plot_options.difftype == 'ratio':
                return f'({ybase2.label} - {ybase1.label}) / ({ybase2.label} + {ybase1.label})'
            elif plot_options.difftype == 'absratio':
                return fr'$|${ybase2.label} - {ybase1.label}$|$ / $|${ybase2.label} + {ybase1.label}$|$'

    def get_smoothing_sigma():
        # """Get the sigma used for Gaussian signal smoothing"""
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
        if controls.etgm_sum_modes.values and 'xte' in var_to_plot1:
            savename = f'{savename}_sum'
        if plot_options.savenameend:
            savename = f'{savename}_{plot_options.savenameend}'
        return f'{savename}_{var_to_plot1}'


    # Plotting loop initialization
    colormaps = plotting.modules.colormaps.get_colormaps()
    options = contour_list[0].options
    options2 = contour_list[1].options
    # print('options check: ', options.scan_num == options2.scan_num)
    var_to_scan = options.var_to_scan
    adjustment_name = options.adjustment_name or options.var_to_scan

    input_vars_rho, output_vars_rho, controls_rho = datahelper.get_all_rho_data(options)
    input_vars_rho2, output_vars_rho2, controls_rho2 = datahelper.get_all_rho_data(options2)
    input_vars, output_vars, controls = datahelper.get_data_objects(options)
    input_vars2, output_vars2, controls2 = datahelper.get_data_objects(options2)

    rho_strs = input_vars_rho.keys()

    if (rho_strs != input_vars_rho2.keys()):
        raise ValueError('Rho values must match from each data set')

    xbase = get_base_data(var_to_scan, output_vars, input_vars, controls)

    x = np.array(list(output_vars_rho.keys()), dtype=float)
    y = options.scan_range
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)
    Z1 = np.zeros_like(X)
    Z2 = np.zeros_like(X)

    fig, ax = plt.gcf(), plt.gca()  # update figure variables for current plot

    # if var_to_plot == 'var_to_scan':
    #     var_to_plot = var_to_scan

    var_to_plot1 = contour_list[0].var_to_plot
    var_to_plot2 = contour_list[1].var_to_plot

    ybase1 = get_base_data(var_to_plot1, output_vars, input_vars, controls)
    ydata1 = get_var_to_plot_data(var_to_plot1, output_vars_rho, input_vars_rho)
    ybase2 = get_base_data(var_to_plot2, output_vars2, input_vars2, controls2)
    ydata2 = get_var_to_plot_data(var_to_plot2, output_vars_rho2, input_vars_rho2)

    # print('base value check: ', (ybase1.values == ybase2.values).all())

    if not isinstance(ybase1.values, np.ndarray) or not isinstance(ybase2.values, np.ndarray):
        # Contours need to be np.ndarray. This can happen when using the
        # same variable plotting list for many different scan types, so
        # no exception is raised
        raise TypeError(f'Either {var_to_plot1} or {var_to_plot2} was not an np.ndarray')

    for i, rho_str in enumerate(rho_strs):
        Z1[:, i] = getattr(ydata1[rho_str], var_to_plot1).values
        Z2[:, i] = getattr(ydata2[rho_str], var_to_plot2).values

    if (Z1 == Z2).all():
        print(f'Identical data sets found for {var_to_plot1}, {var_to_plot2}')
        return

    # Take difference
    Zdiff = Z2 - Z1
    Zsum = np.absolute(Z1 + Z2)
    Zsum[Zsum == 0] = 1

    if plot_options.difftype == 'diff':
        Z = Zdiff
    elif plot_options.difftype == 'absdiff':
        Z = np.absolute(Zdiff)
    elif plot_options.difftype == 'ratio':
        Z = Zdiff / Zsum
    elif plot_options.difftype == 'absratio':
        Z = np.absolute(Zdiff / Zsum)

    if np.isnan(Z).any():
        # This can happen due to calculation errors or or when loading
        # data using old variable definitions, so a warning is thrown but
        # no exception is raised
        raise ValueError(f'nan values found in the computed error between {var_to_plot1}, {var_to_plot2}')

    if not plot_options.savefig:  # Connect key-press handler when not autosaving figures
        fig.canvas.mpl_connect('key_press_event', on_press)

    for s in ax.spines.values():
        s.set_zorder(10)  # Put axes frame on top of everything

    # Apply smoothing and restore original extrema
    Zmax, Zmin = Z.max(), Z.min()
    # Z = scipy.ndimage.gaussian_filter(Z, sigma=get_smoothing_sigma())
    # Z = np.minimum(np.maximum(Z, Zmin), Zmax)

    # Default arguments for filled contours, line contours, and both types
    args_fill = {'zorder': -3}
    args_line = {'zorder': -2, 'linewidths': 0.4}
    args_both = {}

    # Max and min values for Z
    # print(ybase.name, ybase.contour_min, ybase.contour_max, )
    # Zmax = ybase.contour_max
    # Zmin = ybase.contour_min
    # Zmax, Zmin = get_vmax_vmin()

    # Clamp Z between Zmin and Zmax
    # Z = np.minimum(np.maximum(Z, Zmin), Zmax)

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
    plt.title(get_title(plot_options))

    if plot_options.savedata or plot_options.savefig:
        savedir_base = f'{utils.get_plotting_contours_path()}\\{options.runid}'
        utils.create_directory(savedir_base)

    if plot_options.savedata:
        savedir = f'{savedir_base}\\data'
        utils.create_directory(savedir)
        _save_to_csv(X, Y, Z, f'{savedir}\\{get_savename()}')

    if plot_options.savefig:
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


def main(vars_to_plot, scan_data, plot_options):
    '''
    Verifies vars_to_plot, then runs the plotting loop for each var_to_plot, runid, and scan_num

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * scan_data (dict): Dictionary of runid to list of scan numbers
    * plot_options (PlotOptions): Class of plotting options
    '''

    utils.init_logging()
    _verify_vars_to_plot(vars_to_plot)

    if plot_options.savefig:
        print(f'Files will be saved in:\n\t{utils.get_plotting_contours_path()}')

    plt.figure()  # Instantiate the figure

    for var_to_plot in vars_to_plot:
        contour_list = []
        for runid, scan_nums in scan_data.items():
            for scan_num in scan_nums:
                # for var_to_plot in vars_to_plot:
                contour_list.append(ContourPlotData(runid, scan_num, var_to_plot))
                    # if len(contour_list) > 1:
                    #     break
                if len(contour_list) > 1:
                    break
            if len(contour_list) > 1:
                break  

        if len(contour_list) != 2:
            raise ValueError('Need two sets of data to take a percent difference')
        for c in contour_list:
            if not c.options.var_to_scan:
                raise ValueError(f'No variable scan detected for {c.runid}, {c.scan_num}, {c.var_to_plot}')
            if c.var_to_plot == 'time':
                raise ValueError(f'Nothing to plot when {var_to_plot} is time...')
        if contour_list[0].options.var_to_scan != contour_list[1].options.var_to_scan:
            raise ValueError(f'Variable scan is not the same for each data set')            

        plot_contour_difference(contour_list, plot_options)


# Run this file directly to plot scanned variable profiles from previously created scanned data
if __name__ == '__main__':
    scan_data = {}

    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.AIP2,
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
    })

    plot_options = PlotOptions(
        title='',
        savenameend='',
        savefig=False,
        savedata=False,
        difftype='absdiff'
    )

    """
    Variables to Plot (List):
    * Uncomment the line you wish to use
    """

    vars_to_plot = ['xte', 'xti', 'xdi', 'xdz', 'xvt', 'xvp', 'vcz', 'vcp', 'vct']
    # vars_to_plot = ['vci', 'vch', 'vce', 'vcz', 'vct', 'vcp']
    """
    Scan Data:
    * Uncomment the line you wish to use
        - keys (str): The runid of the scan
        - values (list of int): The scan_numbers to plot from
    """

    # scan_data['118341T54'] = [2, 11]  # rmaj/rmaj0
    # scan_data['118341T54'] = [11, 15]  # zz**2/zz**2
    # scan_data['118341T54'] = [15, 16]  # zckb, zcmp in Csound
    # scan_data['118341T54'] = [16, 19]  # 2*zcmu0*zckb in beta, betae
    # scan_data['118341T54'] = [19, 21]  # te * sqrt(te) in vei
    # scan_data['118341T54'] = [21, 27]  # veih fix
    # scan_data['118341T54'] = [2, 27]  # rmaj/rmaj0 to veih fix

    # scan_data['118341T54'] = [27, 54]  # making internal vars into parameters
    # scan_data['118341T54'] = [54, 198]  # min iterations = 1
    # scan_data['118341T54'] = [198, 210]  # requiring an extra convergence iteration
    # scan_data['118341T54'] = [210, 229]  # q2 = q**2
    # scan_data['118341T54'] = [229, 267]  # az_1
    # scan_data['118341T54'] = [267, 287]  # 1E-4 replace 1/10000 
    # scan_data['118341T54'] = [287, 302]  # multiple wde_i definitions
    # scan_data['118341T54'] = [302, 360]  # Temp disable MMM8.1, fixed veih, wde_old, omg_guess
    # scan_data['118341T54'] = [371, 379]  # 361 More Points: 
    # scan_data['118341T54'] = [406, 425]  # Throwing out modes < gma_max / 100  
    # scan_data['118341T54'] = [425, 432]  # vei factor 900 updated to use physical constants
    # scan_data['118341T54'] = [432, 440]  # vei from MMM (including Zeff)
    # scan_data['118341T54'] = [440, 442]  # csound0 = csound
    # scan_data['118341T54'] = [442, 450]  # ion correlation loop
    # scan_data['118341T54'] = [450, 453]  # zflz_i definition
    # scan_data['118341T54'] = [453, 454]  # G_ave_e in diffusivity calculation
    # scan_data['118341T54'] = [454, 463]  # Re-enabled MMM8.1 variables Curr and ALS
    # scan_data['118341T54'] = [463, 469]  # Removed density normalization of 1E-19 
    # scan_data['118341T54'] = [469, 490]  # Mrat using ai instead of az
    # scan_data['118341T54'] = [498, 510]  # Starting with no optimizations

    # scan_data['121123K55'] = [11000, 11002]  # Starting with no optimizations
    # scan_data['120982A09'] = [11000, 11002]  # Starting with no optimizations
    # scan_data['118341T54'] = [11000, 11002]  # Starting with no optimizations
    # scan_data['85126T02']  = [11000, 11002]  # Starting with no optimizations

    scan_data['138536A01']  = [323, 326]  # Starting with no optimizations

    vars_to_plot = [ 'xteMTM', 'xteDBM', 'xtiDBM', 'xteETG', 'xteETGM', 'xte2ETGM', 'xteW20', 'xtiW20', 'xdeW20', 'xdz', 'xvt', 'xvp', 'vcz', 'vcp']
    vars_to_plot = [ 'xtiW20','xteW20', 'xdeW20', 'xdz', 'xvt', 'xvp', 'vcz', 'vcp']
    vars_to_plot = ['gmaEPM']
    # vars_to_plot = ['gmaDBM', 'omgDBM', , 'xtiDBM', 'xdiDBM']

    # scan_data['138536A01'] = [100, 109]  # zcf definition simplification
    # scan_data['138536A01'] = [108, 116]  # *** removing unused nmodes from dribm changes xteMTM output even when dribm isn't enabled!!!
    # scan_data['138536A01'] = [116, 128]  # Fixing redefinition of vthe that was affecting MTM when ETG was enabled
    # scan_data['138536A01'] = [128, 143]  # Updating diffusivity bounds in input file
    # scan_data['138536A01'] = [143, 149]  # kyrhos subroutine (MTM is changing again for no reason!!)
    # scan_data['138536A01'] = [149, 163]  # DRIBM model enabled
    # scan_data['138536A01'] = [163, 175]  # Renaming w20 variables, input rmaj0 instead of eps0
    # scan_data['138536A01'] = [175, 176]  # Reducing kyrhos scan count from 100 to 50 speed scans up
    # scan_data['138536A01'] = [176, 177]  # passing in me_mi instead of calculating it from me_mp / ai
    # scan_data['138536A01'] = [177, 179]  # renaming xdi to xde (breaks comparison)
    # scan_data['138536A01'] = [179, 183]  # If (z_rmin(jz) < 1E-4 * rmaj0) Cycle
    # scan_data['138536A01'] = [183, 222]  # MTM reduced to 100 exp scans, exp switch added
    # scan_data['138536A01'] = [222, 228]  # xteDBM max bound set to 100 from 1000 (DRIBM kyrhos = 0.1 to 1)
    # scan_data['138536A01'] = [228, 231]  # sending unit variables to DRIBM
    # scan_data['138536A01'] = [231, 255]  # zlf_factor added to DRIBM
    # scan_data['138536A01'] = [255, 263]  # testmmm update and prior to removing smoothing from q
    # scan_data['138536A01'] = [157, 158]  #

    # vars_to_plot = ['gmaDBM','gmaETGM']
    # scan_data['118341T54'] = [536, 531]  #
    # plot_options.title, vars_to_plot = r'$|\chi_{\mathrm{i, opt}} - \chi_{\mathrm{i, def}}|$ (m$^2$/s) ', ['xtiW20']
    # plot_options.title, vars_to_plot = r'$|\chi_{\mathrm{e, opt}} - \chi_{\mathrm{e, def}}|$ (m$^2$/s) ', ['xteW20']
    # plot_options.title, vars_to_plot = r'$|\chi_{\mathrm{n, opt}} - \chi_{\mathrm{n, def}}|$ (m$^2$/s) ', ['xdeW20']

    # 498 = all optimizations disabled
    # 499 = min iterations 10
    # 500 = min iterations 1
    # 501 = min iterations 1, step2 +1 = Negligible
    # 502 = min iterations 1, step3 +1 = Noticeable
    # 503 = min iterations 1, step2 +1, step3 -1 = Negligible
    # 504 = min iterations 1, step2 +1, step3 +0 = Negligible

    # scan_data['120968A02'] = [5,6]
    # scan_data['138536A01'] = [i for i in range(1716, 1738 + 1)]
    # scan_data['138536A01'] = [i for i in [*range(1716, 1738 + 1), *range(1756, 1763 + 1)]]

    """
    Plotting Options:
    * savefig: show figure when True, autosave figure without showing it when False
    * savedata: autosave data into CSVs when True
    """

    main(vars_to_plot, scan_data, plot_options)
