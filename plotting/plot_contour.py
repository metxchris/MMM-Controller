#!/usr/bin/python3

# Standard Packages
import sys; sys.path.insert(0, '../')
import logging
import io

# 3rd Party Packages
import numpy as np
import scipy.ndimage
import matplotlib.pyplot as plt
from matplotlib import cm, colors, ticker
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.ticker import StrMethodFormatter, NullFormatter
from PyQt5.QtGui import QImage
from PyQt5.QtWidgets import QApplication

# Local Packages
import settings
import modules.options
import modules.utils as utils
import modules.datahelper as datahelper
import modules.constants as constants
from modules.enums import ScanType, MergeType
from modules.variables import InputVariables, OutputVariables
from plotting.modules.plotstyles import PlotStyles, StyleType
import plotting.modules.colormaps


_log = logging.getLogger(__name__)


def run_plotting_loop(vars_to_plot, options, autosave=False):
    '''
    Creates PDF Plots of each variable in vars_to_plot

    The x-axis variable is pulled as var_to_scan from the saved options file

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * options (Options): Object containing user options
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
    scan_type = options.scan_type

    input_vars_dict, output_vars_dict, input_controls = datahelper.get_all_rho_data(options)
    base_input_vars, base_output_vars, base_input_controls = datahelper.get_data_objects(options)

    xbase = get_base_data(var_to_scan)

    # xvar_data = input_vars_dict[rho_str] if scan_type == ScanType.VARIABLE else input_controls
    # xvar = getattr(xvar_data, var_to_scan)
    # yvar = getattr(output_vars_dict[rho_str], var_to_plot)

    x = np.array(list(output_vars_dict.keys()), dtype=float)
    y = options.scan_range
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)
    G = np.zeros_like(X)

    gamma_min = 1e4

    rho_strs = input_vars_dict.keys()

    for var_to_plot in vars_to_plot:

        levels = 20

        fig = plt.gcf()
        ax = plt.gca()

        if var_to_plot == 'var_to_scan':
            var_to_plot = options.var_to_scan

        if var_to_plot == 'shear':
            var_to_plot = 'shat_gxi'

        print(options.scan_num, var_to_plot)

        for s in ax.spines.values():
            s.set_zorder(10)  # Put axes frame on top of everything

        fig.canvas.mpl_connect('key_press_event', on_press)

        # profile_type = f'{var_to_plot}_{var_to_scan}'
        ybase = get_base_data(var_to_plot)
        ydata = get_var_to_plot_data(var_to_plot)

        if isinstance(ybase.values, float):
            continue

        for i, rho_str in enumerate(rho_strs):
            G[:, i] = getattr(output_vars_dict[rho_str], 'gmaETGM').values
            Z[:, i] = getattr(ydata[rho_str], var_to_plot).values

        if np.isnan(Z).any():
            print(f'nan values found in {var_to_plot}, continuing to next variable...')
            continue

        sigma = 0
        if var_to_plot in ['xteETGM', 'xte2ETGM', 'xdiETGM']:
            sigma = 2
        elif var_to_plot in ['omgETGM', 'kyrhoeETGM', 'kyrhosETGM']:
            sigma = 1

        Zmin = Z.min()
        Zmax = Z.max()
        Z = scipy.ndimage.gaussian_filter(Z, sigma=sigma)
        Z[Z < Zmin] = Zmin
        Z[Z > Zmax] = Zmax

        # np.set_printoptions(threshold=sys.maxsize)
        # print(Z[100, :])

        # if var_to_plot in ['omgETGM', 'xteETGM', 'xte2ETGM', 'kyrhosETGM', 'kyrhoeETGM', 'xdiETGM']:
        #     Z[G < gamma_min] = 0

        norm = None
        # Use gmaETGM as a mask for omgETGM values
        if var_to_plot == 'omgETGM':
            norm = colors.CenteredNorm()
            # levels += 10
            omg_max = 4e5 + 1e-6
            omg_min = 0
            # print(Z.min(), Z.max())

            # Z[G < G.max() / 50] = 0
            # Z[G < gamma_min] = 0
            Z[Z > omg_max] = omg_max
            Z[Z < omg_min] = omg_min
            # print(Z.min())
            cmap_fill = colormaps['magma_both']
            cmap_line = colormaps['binary']
            # divnorm = colors.TwoSlopeNorm(vmin=omg_min, vcenter=0, vmax=omg_max)
            # args_fill = {'norm': divnorm}
            args_line = {'cmap': cmap_line}

            if Z.max() >= omg_max:
                args_fill = {'extend': 'max'}
            # if Z.min() <= omg_min and Z.max() >= omg_max:
            #     args_fill = {'extend': 'both'}
            # elif Z.min() <= omg_min and Z.max() < omg_max:
            #     args_fill = {'extend': 'min'}
            # elif Z.min() > omg_min and Z.max() >= omg_max:
            #     args_fill = {'extend': 'max'}
            # else:
            #     args_fill = {'extend': 'neither'}
            # args_line = {'colors': 'k'}
        # elif Z.min() < 0:
        #     cmap_fill = cmap3
        #     cmap_line = cmap4
        #     divnorm = colors.TwoSlopeNorm(vmin=-5., vcenter=0., vmax=5)
        #     args_fill = {'norm': divnorm}
        #     args_line = {'cmap': cmap_line}

        elif var_to_plot in ['gmaETGM']:
            vmax = 5e5
            vmin = 0
            # Z[Z > vmax] = vmax
            # Z[Z < vmin] = vmin
            cmap_fill = colormaps['magma_both']
            cmap_line = colormaps['binary']
            args_fill = {'extend': 'min'}
            args_fill = {}
            args_line = {'cmap': cmap_line}
            args_line = {'colors': 'k'}
            norm = colors.CenteredNorm()

        elif var_to_plot in ['xteETGM', 'xte2ETGM']:
            norm = colors.CenteredNorm()
            cmap_fill = colormaps['magma_both']
            cmap_line = colormaps['binary']
            args_fill = {'extend': 'neither'}
            args_line = {'cmap': cmap_line}
            # args_line = {'colors': 'k'}

        elif var_to_plot in ['kyrhoeETGM', 'kyrhosETGM']:
            # norm = colors.CenteredNorm()
            # Z[G < gamma_min] = 0
            cmap_fill = colormaps['magma_positive']
            cmap_line = colormaps['binary']
            args_fill = {'extend': 'neither'}
            args_line = {'cmap': cmap_line}

        elif var_to_plot == 'gave':
            # norm = colors.CenteredNorm()
            gave_max = 4
            gave_min = -1e-6
            Z[Z < gave_min] = gave_min
            Z[Z > gave_max] = gave_max

            cmap_fill = colormaps['magma_positive']
            cmap_line = colormaps['binary']

            if Z.max() >= gave_max and Z.min() <= gave_min:
                args_fill = {'extend': 'both'}
            elif Z.max() >= gave_max:
                args_fill = {'extend': 'max'}
            elif Z.min() <= gave_min:
                args_fill = {'extend': 'min'}
            else:
                args_fill = {'extend': 'neither'}

            # if Z.max() >= gave_max:
            #     args_fill = {'extend': 'max'}
            # else:
            #     args_fill = {'extend': 'neither'}

        elif var_to_plot == options.var_to_scan:
            cmap_line = colormaps['binary']
            args_fill = {'extend': 'neither'}
            args_line = {'cmap': cmap_line}

            if Z.min() < 0:
                norm = colors.CenteredNorm()
                cmap_fill = colormaps['magma_both']
                # args_fill = {'extend': 'min'}
            else:
                cmap_fill = colormaps['magma_positive']

        else:
            norm = colors.CenteredNorm()
            cmap_fill = colormaps['magma_both']
            cmap_line = colormaps['binary']
            args_fill = {'extend': 'neither'}
            args_line = {'cmap': cmap_line}

            # if Z.min() < 0:
            #     args_fill = {'extend': 'min'}


            # norm = colors.CenteredNorm()
            # cmap_fill = cmap3
            # cmap_line = cmap2
            # args_fill = {'extend': 'neither'}
            # args_line = {'cmap': cmap_line}

        # if var_to_plot == 'gte':
        #     norm = colors.LogNorm(vmin=Z.min(), vmax=Z.max())
        #     norm = colors.PowerNorm(1)

        # # Plot filled contour areas
        # # define the colormap
        # cmap = plt.get_cmap('PuOr')

        # # extract all colors from the .jet map
        # cmaplist = [cmap(i) for i in range(cmap.N)]
        # # create the new map
        # cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
        #
        # cmap_fill = cmap3
        # Z = norm(Z)

        # norm=colors.BoundaryNorm(np.linspace(1, Z.max(), levels), ncolors=cmap_fill.N)
        # norm=colors.PowerNorm(0.5)
        # norm2=colors.BoundaryNorm(np.linspace(1, Z.max(), levels), ncolors=cmap_line.N)

        cf = plt.contourf(X, Y, Z, cmap=cmap_fill, levels=levels, zorder=3, **args_fill, norm=norm)
        # print(cmap_fill)
        # print(cf.levels)

        lowest_bucket = Z[Z <= cf.levels[1]].size / Z.size
        # print(lowest_bucket)
        if (var_to_plot == options.var_to_scan or var_to_plot == 'shat_gxi') and lowest_bucket > 0.3:
            power = 1.2
            cf.levels = cf.levels**power / (cf.levels**power)[-1] * cf.levels[-1]
            cf.levels = np.round(cf.levels, 1)
            cf = plt.contourf(X, Y, Z, cmap=cmap_fill, levels=cf.levels, zorder=3, **args_fill, norm=norm)
            # print(cf.levels)


        # lowest_bucket = Z[Z <= cf.levels[1]].size / Z.size
        # print(lowest_bucket)
        # if lowest_bucket > 0.2 and cf.levels[0] >= 0:
        #     cf.levels = np.hstack((cf.levels[0], (cf.levels[0] + cf.levels[1]) / 2, cf.levels[1:]))
        #     cf = plt.contourf(X, Y, Z, cmap=cmap_fill, levels=cf.levels, zorder=3, **args_fill, norm=norm)

        if cf.levels[0] == Z.min():
            cf.levels = cf.levels[1:]

        # Colorbar for filled contour colormap
        cb = plt.colorbar(spacing='proportional', pad=0.02, aspect=30, fraction=0.1, drawedges=True)
        # plt.clim(0, 100)
        cb.ax.ticklabel_format(axis="y", style="sci", scilimits=(-1, 2))
        cb.ax.tick_params(size=0, labelsize=6.5)


        if var_to_plot == 'gmaETGM':
            gamma_min = cb.ax.get_yticks().min()

        if 'extend' in args_fill:
            if args_fill['extend'] == 'max' or args_fill['extend'] == 'both':
                cb.ax.yaxis.OFFSETTEXTPAD += 5
        # cb.ax.yaxis.set_major_locator(plt.MaxNLocator(8))  # number of tick labels

        # Plot contour lines
        CS = plt.contour(X, Y, Z, levels=cf.levels, linewidths=0.5, alpha=0.3, zorder=4, colors='k', **args_fill)#, norm=norm2)
        # print(CS.levels)


        # print(matplotlib.ticker.LinearLocator())

        # if Z.min() < 0:
        #     # for line in CS.collections:
        #     #     print(line.get_linestyle())
        #     #     if line.get_linestyle() == [(None, None)]:
        #     #         print("Solid Line")
        #     #     else:
        #     #         line.set_linestyle('--')
        #     #         line.set_color('red')
        # if 'colors' in args_line:
        #     plt.clabel(CS, levels=[0], fontsize=plt.rcParams["legend.fontsize"])
        #     CS = plt.contour(X, Y, Z, cmap=cmap_line, levels=[0], linewidths=0.5, alpha=0.8, zorder=4)
        #     plt.clabel(CS, levels=[0], fontsize=plt.rcParams["legend.fontsize"])

        # Base value line
        # plt.hlines(1, 0, 1, color="#333", lw=0.5, alpha=0.2, ls='--', zorder=2)

        # ax.set(yscale='log')
        # ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.1f}'))
        ax.yaxis.set_minor_formatter(NullFormatter())

        plt.xlabel(r'$\rho$')
        plt.ylabel(f'{xbase.label} (multipliers)')

        title = f'{ybase.label} ({ybase.units_label})' if ybase.units_label else f'{ybase.label}'

        if options.use_gnezero:
            title = fr'{title} [$g_\mathrm{{ne}} = 0$]'

        if options.use_gtezero:
            title = fr'{title} [$g_\mathrm{{Te}} = 0$]'

        if options.use_gneabs:
            title = fr'{title} with $|g_\mathrm{{ne}}|$'

        # if base_input_controls.etgm_diffusivity_type.values == 1:
        #     title = fr'{title} [Alternate]'

        # if base_input_controls.etgm_kyrhoe_scan.values > 0:
        #     title = fr'{title} with $k_y\rho_\mathrm{{e}}$ scan'

        # title = fr'{title} with $4\times g_\mathrm{{Te}}$'

        plt.title(title)

        if autosave:
            savename = f'{var_to_scan}_{var_to_plot}'
            if options.use_gnezero:
                savename = f'{savename}_gne0'
            if options.use_gtezero:
                savename = f'{savename}_gte0'
            if options.use_gneabs:
                savename = f'{savename}_gneabs'
            # if base_input_controls.etgm_kyrhoe_scan.values == 0:
            #     savename = f'{savename}_nokyrhoescan'
            # if base_input_controls.etgm_diffusivity_type.values == 1:
            #     savename = f'{savename}_xte2'

            fig.savefig(f'{utils.get_plotting_contours_path()}\\{savename}')
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
            raise NameError(f'{var_to_plot} not found in OutputVariables class')


def main(vars_to_plot, scan_data, autosave=False):
    '''
    Verifies vars_to_plot, then runs the plotting loop for each var_to_plot, runid, and scan_num

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * scan_data (dict): Dictionary of runid to list of scan numbers
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
            utils.clear_temp_folder(options)
            if options.var_to_scan:
                run_plotting_loop(vars_to_plot, options, autosave)
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
        'contour.negative_linestyle': 'dotted',
        'pdf.compression': 6,
    })

    '''
    Input options:
    * vars_to_plot (list): List of output variables to plot

    Examples:
    * vars_to_plot = ['gmaETGM']
    * vars_to_plot = ['xteMTM', 'xteETGM', 'xteETG', 'gmaMTM', 'omgMTM', 'dbsqprf']
    * vars_to_plot = OutputVariables().get_all_output_vars()
    * vars_to_plot = OutputVariables().get_etgm_vars()
    '''
    vars_to_plot = ['gmaETGM', 'omgETGM', 'xteETGM', 'xte2ETGM', 'kyrhoeETGM', 'kyrhosETGM', 'gave', 'var_to_scan']
    # vars_to_plot = ['gmaETGM', 'omgETGM', 'xteETGM']
    # vars_to_plot = ['gmaETGM', 'omgETGM',]
    # vars_to_plot = ['gmaETGM', 'kyrhoeETGM', 'kyrhosETGM',]
    # vars_to_plot = ['xteETGM']
    # vars_to_plot = ['gmaETGM', 'xteETGM', 'xte2ETGM']

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
    # scan_data['138536A01'] = [i for i in range(218, 231)]
    # scan_data['138536A01'] = [i for i in range(231, 244)]
    # scan_data['138536A01'] = [i for i in range(265, 290)]
    # scan_data['138536A01'] = [i for i in range(298, 310)]
    scan_data['138536A01'] = [i for i in range(387, 398)]
    scan_data['138536A01'] = [i for i in range(433, 450)]
    # scan_data['138536A01'] = [i for i in range(444, 450)]
    # scan_data['138536A01'] = [406, 407, 408]
    scan_data['138536A01'] = [500]

    main(vars_to_plot, scan_data, autosave=0)
