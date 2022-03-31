#!/usr/bin/python3

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


def run_plotting_loop(vars_to_plot, options, savenameend='', autosave=False):
    '''
    Creates PDF Plots of each variable in vars_to_plot

    The x-axis variable is pulled as var_to_scan from the saved options file

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * options (Options): Object containing user options
    * savenameend (str): Name to end file save name with (Optional)
    * autosave (bool): Automatically save the plot if True (Optional)
    '''

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
        obj = None
        if hasattr(base_output_vars, name):
            obj = getattr(base_output_vars, name)
        elif hasattr(base_input_vars, name):
            obj = getattr(base_input_vars, name)
        elif hasattr(base_input_controls, name):
            obj = getattr(base_input_controls, name)

        return obj

    def get_var_to_plot_data(name):
        obj = None
        if hasattr(base_output_vars, name):
            obj = output_vars_dict
        elif hasattr(base_input_vars, name):
            obj = input_vars_dict

        return obj

    colormaps = plotting.modules.colormaps.get_colormaps()
    var_to_scan = options.var_to_scan
    adjustment_name = options.adjustment_name or options.var_to_scan

    input_vars_dict, output_vars_dict, input_controls = datahelper.get_all_rho_data(options)
    base_input_vars, base_output_vars, base_input_controls = datahelper.get_data_objects(options)

    xbase = get_base_data(var_to_scan)

    x = np.array(list(output_vars_dict.keys()), dtype=float)
    y = options.scan_range
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)

    levels = 20

    rho_strs = input_vars_dict.keys()

    for var_to_plot in vars_to_plot:

        fig = plt.gcf()
        ax = plt.gca()

        if var_to_plot == 'var_to_scan':
            var_to_plot = var_to_scan

        if var_to_plot == 'time':
            continue  # nothing to plot

        # if var_to_plot == 'shear':
        #     var_to_plot = 'shat_gxi'

        print(options.scan_num, var_to_plot)

        ybase = get_base_data(var_to_plot)
        ydata = get_var_to_plot_data(var_to_plot)

        if isinstance(ybase.values, float):
            continue

        ylabel = xbase.label
        if var_to_scan == 'gne' and options.use_gneabs:
            ylabel = f'$|${ylabel}$|$'

        if var_to_scan not in ['time', 'etgm_kyrhos_min']:
            ylabel = f'{ylabel} (multipliers)'
        elif xbase.units_label:
            ylabel = f'{ylabel} ({xbase.units_label})'

        title = f'{ybase.label} ({ybase.units_label})' if ybase.units_label else f'{ybase.label}'
        if options.use_gnezero:
            title = fr'{title} [$g_\mathrm{{ne}} = 0$]'
        if options.use_gtezero:
            title = fr'{title} [$g_\mathrm{{Te}} = 0$]'
        if options.use_gneabs:
            title = fr'{title} with $|g_\mathrm{{ne}}|$'
        if base_input_controls.etgm_exbs.values and var_to_plot != 'wexbs':
            title = fr'{title} $[\omega_{{E \times\! B}}\,\,\mathrm{{on}}]$'
        if base_input_controls.etgm_sum_modes.values and 'xte' in var_to_plot:
            title = fr'{title} [Sum]'

        for i, rho_str in enumerate(rho_strs):
            Z[:, i] = getattr(ydata[rho_str], var_to_plot).values

        # if var_to_plot == 'wexbs':
        #     Z = np.absolute(Z)

        # # TODO: Hack to calculate omegasETGM
        # if var_to_plot == 'omegasETGM':
        #     GNE = np.zeros_like(Z)

        #     for i, rho_str in enumerate(rho_strs):
        #         Z[:, i] = getattr(output_vars_dict[rho_str], 'omegadETGM').values
        #         GNE[:, i] = getattr(input_vars_dict[rho_str], 'gne').values

        #     Z *= 0.5 * GNE

        # # TODO: Hack to calculate omegasetaETGM
        # if var_to_plot == 'omegasetaETGM':
        #     ETAE = np.zeros_like(Z)
        #     GNE = np.zeros_like(Z)

        #     for i, rho_str in enumerate(rho_strs):
        #         Z[:, i] = getattr(output_vars_dict[rho_str], 'omegadETGM').values
        #         GNE[:, i] = getattr(input_vars_dict[rho_str], 'gne').values
        #         ETAE[:, i] = getattr(input_vars_dict[rho_str], 'etae').values

        #     Z *= 0.5 * GNE * (1 + ETAE)

        # # TODO: Hack to calculate omegateETGM
        # if var_to_plot == 'omegateETGM':
        #     GTE = np.zeros_like(Z)

        #     for i, rho_str in enumerate(rho_strs):
        #         Z[:, i] = getattr(output_vars_dict[rho_str], 'omegadETGM').values
        #         GTE[:, i] = getattr(input_vars_dict[rho_str], 'gte').values

        #     Z *= 0.5 * GTE

        # # TODO: Hack to calculate omegadiffETGM
        # if var_to_plot == 'omegadiffETGM':
        #     OMG = np.zeros_like(Z)
        #     OMGD = np.zeros_like(Z)
        #     GAVE = np.zeros_like(Z)

        #     for i, rho_str in enumerate(rho_strs):
        #         OMG[:, i] = getattr(output_vars_dict[rho_str], 'omgETGM').values
        #         OMGD[:, i] = getattr(output_vars_dict[rho_str], 'omegadETGM').values
        #         GAVE[:, i] = getattr(output_vars_dict[rho_str], 'gaveETGM').values

        #     Z = (OMG - 5 / 3 * OMGD * GAVE)

        # # TODO: Hack to calculate gammadiffETGM
        # if var_to_plot == 'gammadiffETGM':
        #     GMA = np.zeros_like(Z)
        #     OMGD = np.zeros_like(Z)
        #     GAVE = np.zeros_like(Z)

        #     for i, rho_str in enumerate(rho_strs):
        #         GMA[:, i] = getattr(output_vars_dict[rho_str], 'gmaETGM').values
        #         OMGD[:, i] = getattr(output_vars_dict[rho_str], 'omegadETGM').values
        #         GAVE[:, i] = getattr(output_vars_dict[rho_str], 'gaveETGM').values

        #     Z = (GMA - 5 / 3 * OMGD * GAVE)

        if np.isnan(Z).any():
            print(f'nan values found in {var_to_plot}, continuing to next variable...')
            continue

        fig.canvas.mpl_connect('key_press_event', on_press)
        for s in ax.spines.values():
            s.set_zorder(10)  # Put axes frame on top of everything

        # Smoothing weights
        medium_smoothing = ['xteETGM', 'xte2ETGM', 'xdiETGM']
        light_smoothing = [
            'omgETGM', 'kyrhoeETGM', 'kyrhosETGM', 'omegateETGM', 'walfvenunit',
            'omegadETGM', 'omegasETGM', 'omegasetaETGM', 'omegadiffETGM', 'gammadiffETGM'
        ]

        sigma = (0, 0)
        if options.var_to_scan == 'time':
            sigma = (0, 0)
        elif var_to_plot in medium_smoothing:
            sigma = (2, 2)
        elif var_to_plot in light_smoothing:
            sigma = (1, 1)

        # Apply smoothing and restore original extrema
        Zmax, Zmin = Z.max(), Z.min()
        Z = scipy.ndimage.gaussian_filter(Z, sigma=sigma)
        Z = np.minimum(np.maximum(Z, Zmin), Zmax)

        # Default values
        vmax, vmin = np.inf, -np.inf
        args_fill = {'zorder': -3}
        # args_line_base = {'zorder': -3, 'linewidths': 0.1}
        args_line = {'zorder': -2, 'linewidths': 0.4}
        args_both = {}

        # Max and min values per variable
        if 'gma' in var_to_plot:
            vmax, vmin = 5e5, -5e5
            vmax, vmin = 5e7, -5e7
        elif 'omg' in var_to_plot:
            vmax, vmin = 5e5, -5e5
            vmax, vmin = 3e6, -3e6
        elif 'xteETGM' == var_to_plot:
            vmax, vmin = 2e1, -2e1
            if base_input_controls.etgm_sum_modes.values:
                (vmax, vmin) = (1e6, -1e6) if var_to_scan != 'time' else (1e6, -1e6)
                # vmax, vmin = 1e4, -1e4
            else:
                (vmax, vmin) = (1e6, -1e6) if var_to_scan != 'time' else (2e1, -2e1)
        elif 'xte2ETGM' == var_to_plot:
            vmax, vmin = 1e2, -1e2
            if base_input_controls.etgm_sum_modes.values:
                (vmax, vmin) = (1e8, -1e8) if var_to_scan != 'time' else (1e8, -1e8)
                # vmax, vmin = 1e5, -1e5
            else:
                (vmax, vmin) = (1e8, -1e8) if var_to_scan != 'time' else (2e1, -2e1)
        elif 'xte' in var_to_plot:
            vmax, vmin = 5e2, -5e2
        elif 'xti' in var_to_plot:
            vmax, vmin = 5e2, -5e2
        elif 'xdi' in var_to_plot:
            vmax, vmin = 5e2, -5e2
        elif 'gave' in var_to_plot:
            vmax, vmin = 5, 0
        elif 'etae' == var_to_plot:
            vmax, vmin = 20, -20
        # elif 'kyrhosETGM' == var_to_plot:
        #     vmax, vmin = 50, 0
        elif 'shat' in var_to_plot or 'shear' in var_to_plot:
            vmax, vmin = 20, -20
        elif 'omegadETGM' in var_to_plot:
            vmax, vmin = 1e7, -1e7
        elif 'omegateETGM' in var_to_plot:
            vmax, vmin = 1e7, -1e7
        elif 'omegasETGM' in var_to_plot:
            vmax, vmin = 1e7, -1e7
        elif 'omegasetaETGM' in var_to_plot:
            vmax, vmin = 1e7, -1e7
        elif 'omegadiffETGM' in var_to_plot:
            vmax, vmin = 1e6, -1e6
        elif 'gammadiffETGM' in var_to_plot:
            vmax, vmin = 1e6, -1e6

        # Clamp Z between vmin and vmax
        Z = np.minimum(np.maximum(Z, vmin), vmax)

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
        if Z.max() >= vmax and Z.min() <= vmin:
            args_both['extend'] = 'both'
        elif Z.max() >= vmax:
            args_both['extend'] = 'max'
        elif Z.min() <= vmin:
            args_both['extend'] = 'min'
        else:
            args_both['extend'] = 'neither'

        # Set number of contour levels
        loc = ticker.MaxNLocator(levels + 1, min_n_ticks=levels - 4)
        lvls = loc.tick_values(Z.min(), Z.max())

        # if var_to_plot == 'gmaETGM':
        #     print(lvls)
        #     lvls = np.sort(np.append(lvls, 1e5))
        #     print(lvls)

        # lev_exp = np.arange(np.floor(np.log2(Z.min())-1),
        #                    np.ceil(np.log2(Z.max())+1))
        # lvls = np.power(2, lev_exp)
        # args_both['norm'] = colors.LogNorm(base=2)
        # args_both['norm'] = None
        # print(max(np.power(10, np.floor(np.log10(Z.min()) - 1)), 1))
        # print(np.power(10, np.ceil(np.log10(Z.max()) + 1)))
        # loc = ticker.LogLocator(base=10)
        # lvls = loc.tick_values(max(np.power(10, np.floor(np.log10(Z.min()) - 1)), 1), np.power(10, np.ceil(np.log10(Z.max()) + 1)))
        # args_both['locator'] = ticker.LogLocator(base=10)
        # locator=ticker.LogLocator()

        # lowest_bucket = Z[Z <= lvls[1]].size / Z.size

        # Condition for nonlinear contour levels
        # if (var_to_plot == options.var_to_scan or var_to_plot == 'shat_gxi') and lowest_bucket > 0.4:
        #     power = 1.2
        #     lvls = np.round(lvls**power / (lvls**power)[-1] * lvls[-1], 3)

        # Plot contour fills
        cf = plt.contourf(X, Y, Z, levels=lvls, **args_fill, **args_both)

        # Remove lowest level when it is also a boundary, so an extra line isn't drawn
        if cf.levels[0] == Z.min():
            if args_both['extend'] != 'min' and args_both['extend'] != 'both':
                cf.levels = cf.levels[1:]

        # Colorbar for filled contour colormap
        cb = plt.colorbar(spacing='proportional', pad=0.02, aspect=30, fraction=0.1, drawedges=True)
        cb.ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        cb.ax.tick_params(size=0, labelsize=6.5)

        # An upper colorbar extension moves the offset text down, so we raise it back up
        if 'extend' in args_both:
            if args_both['extend'] == 'max' or args_both['extend'] == 'both':
                cb.ax.yaxis.OFFSETTEXTPAD += 5

        # Plot contour lines
        # plt.contour(X, Y, Z, levels=cf.levels, **args_line_base, **args_both)
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
        plt.ylabel(ylabel)
        plt.title(title)

        if autosave:
            savename = adjustment_name
            if options.use_gnezero:
                savename = f'{savename}_gne0'
            if options.use_gtezero:
                savename = f'{savename}_gte0'
            if options.use_gneabs:
                savename = f'{savename}_gneabs'
            if base_input_controls.etgm_exbs.values:
                savename = f'{savename}_exbs1'
            if base_input_controls.etgm_sum_modes.values and 'xte' in var_to_plot:
                savename = f'{savename}_sum'
            if savenameend:
                savename = f'{savename}_{savenameend}'

            savename = f'{savename}_{var_to_plot}'

            fig.savefig(f'{utils.get_plotting_contours_path()}\\{options.runid}\\{savename}')
        else:
            plt.show()

        fig.clear()

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
    input_vars = InputVariables()
    for var_to_plot in vars_to_plot:
        if var_to_plot == 'var_to_scan':
            continue
        if not hasattr(output_vars, var_to_plot) and not hasattr(input_vars, var_to_plot):
            raise NameError(f'{var_to_plot} not found in Variables classes')


def main(vars_to_plot, scan_data, savenameend='', autosave=False):
    '''
    Verifies vars_to_plot, then runs the plotting loop for each var_to_plot, runid, and scan_num

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * scan_data (dict): Dictionary of runid to list of scan numbers
    * savenameend (str): Name to end file save name with
    * autosave (bool): Automatically save the plot if True (Optional)
    '''

    utils.init_logging()
    verify_vars_to_plot(vars_to_plot)
    options = modules.options.Options()

    if autosave:
        print(f'Files will be saved in:\n\t{utils.get_plotting_contours_path()}')

    plt.figure()

    for runid, scan_nums in scan_data.items():
        for scan_num in scan_nums:
            print(f'Initializing data for {runid}, scan {scan_num}...')
            options.load(runid, scan_num)
            utils.create_directory(f'{utils.get_plotting_contours_path()}\\{runid}')
            utils.clear_temp_folder(options)
            if options.var_to_scan:
                run_plotting_loop(vars_to_plot, options, savenameend, autosave)
                utils.clear_temp_folder(options)
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
        'savefig.format': 'pdf',
        'pdf.compression': 6,
    })

    savenameend = ''

    '''
    Input options:
    * vars_to_plot (list): List of output variables to plot

    Examples:
    * vars_to_plot = ['gmaETGM']
    * vars_to_plot = ['xteMTM', 'xteETGM', 'xteETG', 'gmaMTM', 'omgMTM', 'dbsqprf']
    * vars_to_plot = OutputVariables().get_all_output_vars()
    * vars_to_plot = OutputVariables().get_etgm_vars()
    '''
    # vars_to_plot = ['gmaETGM', 'omgETGM', 'xteETGM', 'xte2ETGM', 'kyrhoeETGM', 'kyrhosETGM', 'gave', 'var_to_scan']
    # vars_to_plot = ['gmaETGM', 'omgETGM', 'xteETGM']
    # vars_to_plot = ['var_to_scan', 'gmaETGM', 'lareunit', 'alphamhdunit', 'xteETGM', 'xte2ETGM', 'gaveETGM']
    # vars_to_plot = ['var_to_scan']

    '''
    Scan Data:
    * Uncomment the lines you wish to use
        - keys (str): The runid of the scan
        - values (list of int): The scan_numbers to plot from
    '''
    """exbs off"""
    vars_to_plot = [
        'gmaETGM', 'omgETGM', 'xteETGM', 'xte2ETGM',
        'xteETG', 'walfvenunit', 'phi2ETGM', 'Apara2ETGM', 'satETGM',
        'gaveETGM', 'kyrhosETGM', 'kyrhoeETGM', 'kpara2ETGM', 'fleETGM', 'omegateETGM',
        'omegadETGM', 'omegasETGM', 'omegasetaETGM', 'omegadiffETGM', 'gammadiffETGM',
        'gne', 'gte', 'shat_gxi', 'etae', 'betaeunit', 'wexbs', 'bunit', 'te', 'ne', 'q'
    ]
    """SUMMED MODES"""
    # scan_data['121123K55'] = [4]
    # scan_data['120968A02'] = [4]
    # scan_data['120982A09'] = [4]
    # scan_data['129016A04'] = [4]
    # scan_data['129017A04'] = [4]
    # scan_data['129018A02'] = [4]
    # scan_data['129019A02'] = [4]
    # scan_data['129020A02'] = [4]
    # scan_data['129041A10'] = [4]
    # scan_data['138536A01'] = [1475]  # 1236 (7), 1256 (8), 1269 (9)
    # scan_data['141007A10'] = [4]
    # scan_data['141031A01'] = [4]
    # scan_data['141032A01'] = [4]
    # scan_data['141040A01'] = [4]
    # scan_data['141716A80'] = [4]
    # scan_data['132017T01'] = [4]
    # scan_data['141552A01'] = [4]
    """MAX MODE"""
    # scan_data['121123K55'] = [5]
    # scan_data['120968A02'] = [5]
    # scan_data['120982A09'] = [5]
    # scan_data['129016A04'] = [5]
    # scan_data['129017A04'] = [5]
    # scan_data['129018A02'] = [5]
    # scan_data['129019A02'] = [5]
    # scan_data['129020A02'] = [5]
    # scan_data['129041A10'] = [5]
    # scan_data['138536A01'] = [1476]  # 1236 (7), 1256 (8), 1269 (9)
    # scan_data['141007A10'] = [5]
    # scan_data['141031A01'] = [5]
    # scan_data['141032A01'] = [5]
    # scan_data['141040A01'] = [5]
    # scan_data['141716A80'] = [5]
    # scan_data['132017T01'] = [5]
    # scan_data['141552A01'] = [5]

    """exbs on"""
    # vars_to_plot = ['gmaETGM', 'xteETGM', 'xte2ETGM']
    # scan_data['121123K55'] = [2]
    # scan_data['120968A02'] = [2]
    # scan_data['120982A09'] = [2]
    # scan_data['129016A04'] = [2]
    # scan_data['129017A04'] = [2]
    # scan_data['129018A02'] = [2]
    # scan_data['129019A02'] = [2]
    # scan_data['129020A02'] = [2]
    # scan_data['129041A10'] = [2]
    # scan_data['138536A01'] = [1237]
    # scan_data['141007A10'] = [2]
    # scan_data['141031A01'] = [2]
    # scan_data['141032A01'] = [2]
    # scan_data['141040A01'] = [2]
    # scan_data['141716A80'] = [2]
    # scan_data['132017T01'] = [2]
    # scan_data['141552A01'] = [2]

    # vars_to_plot = [
    #     'gmaETGM', 'omgETGM', 'xteETGM', 'xte2ETGM', 'xteETG', 'kyrhoeETGM', 'kyrhosETGM', 'gaveETGM', 'kpara2ETGM', 'fleETGM',
    #     'omegadETGM', 'var_to_scan', #'lareunit', 'alphamhdunit'
    # ]

    # vars_to_plot = [
    #     'gmaETGM', 'omgETGM', 'xteETGM', 'xte2ETGM',  # 'xteETG',
    #     'gaveETGM', 'alphaETGM', 'kyrhosETGM', 'kyrhoeETGM', 'kpara2ETGM',  # 'fleETGM',
    #     'omegadETGM', 'omegasETGM', 'omegadiffETGM', 'gammadiffETGM', 'omegasetaETGM',
    # ]
    # vars_to_plot = ['gmaETGM']
    # scan_data['138536A01'] = [i for i in range(433, 450)]

    # scan_data['TEST'] = [554]
    # scan_data['138536A01'] = [1779, 1780]
    scan_data['138536A01'] = [1815]
    # scan_data['138536A01'] = [i for i in range(1716, 1738 + 1)]
    # scan_data['138536A01'] = [i for i in range(1749, 1750 + 1)]
    # scan_data['138536A01'] = [i for i in range(1756, 1763 + 1)]
    # scan_data['138536A01'] = [i for i in range(1779, 1780 + 1)]
    # scan_data['138536A01'] = [i for i in [*range(1716, 1738 + 1), *range(1749, 1750 + 1), *range(1756, 1763 + 1), *range(1779, 1780 + 1)]]

    # scan_data['138536A01'] = [i for i in range(569, 575)]
    # vars_to_plot = [
    #     'omegadiffETGM', 'gammadiffETGM',
    # ]
    # scan_data['138536A01'] = [1268]  # 1264 = no zamr(3,2), 1265 = no zbmr(3,4), 1266 no zamr or zbmr

    # scan_data['138536A01'], savenameend = [632], 'alpha'
    # scan_data['138536A01'], savenameend = [633], 'gave'
    # scan_data['138536A01'], savenameend = [634], 'kpara'
    # scan_data['138536A01'], savenameend = [635], 'gave'

    main(vars_to_plot, scan_data, savenameend=savenameend, autosave=0)
