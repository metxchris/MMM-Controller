#!/usr/bin/python3

# Standard Packages
import logging
import io

# 3rd Party Packages
import numpy as np
from matplotlib import colors, ticker
from matplotlib.ticker import NullFormatter
import matplotlib.pyplot as plt
from PyQt5.QtGui import QImage
from PyQt5.QtWidgets import QApplication

# Local Packages
import modules.utils as utils
import plotting.modules.colormaps
from modules.variables import OutputVariables


_log = logging.getLogger(__name__)


def make_contour_plot(cd):

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

    def set_colorbar_extensions(cd, args_both):
        """Set colorbar extensions"""
        if cd.plot_options.difftype:
            args_both['extend'] = 'both'  # Gets white to show up for zero abs diff
            return

        Zmax, Zmin = cd.zvar.contour_max, cd.zvar.contour_min
        if cd.Z.max() >= Zmax and cd.Z.min() <= Zmin:
            args_both['extend'] = 'both'
        elif cd.Z.max() >= Zmax:
            args_both['extend'] = 'max'
        elif cd.Z.min() <= Zmin:
            args_both['extend'] = 'min'
        else:
            args_both['extend'] = 'neither'

    def set_colormaps(cd, args_both, args_fill, args_line):
        """Set colormaps for contour fills and lines"""
        args_fill['zorder'] = -3
        args_line['zorder'] = -2
        args_line['linewidths'] = 0.4
        colormaps = plotting.modules.colormaps.get_colormaps()

        if hasattr(cd, 'isidentical') and cd.isidentical or (cd.Z == 0).all():
            args_both['norm'] = None
            args_fill['cmap'] = colormaps['white']
            args_line['cmap'] = colormaps['white']
        elif cd.Z.min() >= 0:
            args_both['norm'] = None
            args_fill['cmap'] = colormaps['magma_positive']
            args_line['cmap'] = colormaps['magma_positive_lines']
        elif cd.Z.max() <= 0:
            args_both['norm'] = None
            args_fill['cmap'] = colormaps['magma_negative']
            args_line['cmap'] = colormaps['magma_negative_lines']
        else:
            args_both['norm'] = colors.CenteredNorm()
            args_fill['cmap'] = colormaps['magma_both']
            args_line['cmap'] = colormaps['magma_both_lines']

    def set_contour_levels(cd, args_fill, args_both, use_log=0):
        """Set the displayed contour levels"""

        # Average number of contours to show (actual number will vary)
        contour_level_count = 20
        Z = cd.Z

        if cd.zvar.is_int and Z.max() - Z.min() <= contour_level_count:
            Zmax = int(Z.max())
            Zmin = int(Z.min())
            if Zmax == Zmin:
                Zmax += 1
                Zmin -= 1
            Zcount = max(min(Zmax - Zmin + 1, contour_level_count), 2)
            args_fill['levels'] = np.linspace(Zmin, Zmax, Zcount)
        elif hasattr(cd, 'isidentical') and cd.isidentical or (cd.Z == 0).all():
            loc = ticker.MaxNLocator(1, min_n_ticks=1)
            args_fill['levels'] = loc.tick_values(Z.min(), Z.max())
        elif cd.plot_options.difftype:
            loc = ticker.MaxNLocator(contour_level_count + 1, min_n_ticks=contour_level_count - 4)
            args_fill['levels'] = loc.tick_values(Z.min(), Z.max())
        elif not use_log:
            unique_levels = len(np.unique(np.round(Z, 3)))
            contour_level_count = min(contour_level_count, unique_levels)
            if contour_level_count == unique_levels:
                min_n_ticks = contour_level_count + 0
            else:
                min_n_ticks = contour_level_count - 4
            min_n_ticks = max(min_n_ticks, 2)
            loc = ticker.MaxNLocator(contour_level_count + 0, min_n_ticks=min_n_ticks)
            args_fill['levels'] = loc.tick_values(Z.min(), Z.max())
        else:
            loc = ticker.LogLocator(numticks=contour_level_count + 0, base=10)  # set max ticks
            tick_min = np.power(10, np.floor(np.log10(Z.min())))
            tick_max = np.power(10, np.ceil(np.log10(Z.max())))
            levels = loc.tick_values(tick_min, tick_max)
            levels = levels[levels >= tick_min]
            levels = levels[levels <= tick_max]
            args_fill['levels'] = levels
            args_both['locator'] = loc

    def plot_contours():
        # Plot contour fills
        cf = plt.contourf(cd.X, cd.Y, cd.Z, **args_fill, **args_both)

        # Remove lowest level when it is also a boundary, so an extra contour line isn't drawn
        if cf.levels[0] == cd.Z.min():
            if args_both['extend'] != 'min' and args_both['extend'] != 'both':
                cf.levels = cf.levels[1:]

        # Colorbar for filled contour colormap
        cb = plt.colorbar(spacing='uniform', pad=0.02, aspect=30, fraction=0.1, drawedges=True)
        cb.ax.tick_params(size=0, labelsize=plt.rcParams['ytick.labelsize'] - 0.5)

        if cd.use_log():
            ticks = [utils.get_power_10(t) for t in cf.levels]
            cb.ax.axes.set_yticklabels(ticks)
        else:
            scilimits = (-1, 2) if min(cf.levels) < 0 else (-2, 2)
            cb.ax.ticklabel_format(axis="y", style="sci", scilimits=scilimits)

        # An upper colorbar extension moves the offset text down, so we raise it back up
        if 'extend' in args_both:
            if args_both['extend'] == 'max' or args_both['extend'] == 'both':
                cb.ax.yaxis.OFFSETTEXTPAD += 5

        # Plot contour lines
        cl = plt.contour(cd.X, cd.Y, cd.Z, levels=cf.levels, **args_line, **args_both)

        # Override the linestyles based on the levels.
        for line, lvl in zip(cl.collections, cl.levels):
            if lvl < 0:
                line.set_linestyle('--')
            elif lvl == 0:
                line.set_linestyle('-')
            else:
                line.set_linestyle('-')

    def init_figure(plt, plot_options):
        """Initializes the figure"""

        fig, ax = plt.gcf(), plt.gca()  # update figure variables for current plot

        for s in ax.spines.values():
            s.set_zorder(10)  # Put axes frame on top of everything

        if plot_options.showfig:  # Connect key-press handler when not autosaving figures
            fig.canvas.mpl_connect('key_press_event', on_press)

        return fig, ax

    def write_wexb_text(fig, cd):
        """Writes whether wexb is enabled or disabled to the plot"""
        if cd.controls is None or not hasattr(OutputVariables(cd.options), cd.var_to_plot) or 'MTM' in cd.var_to_plot:
            return  # Don't print wexb status for non-output variables or MTM output variables

        if cd.options.wexb_factor or cd.controls.etgm_exbs.values:
            wexb_status = 'on'
            text_x, text_y = 0.8, 0.035
        else:
            wexb_status = 'off'
            text_x, text_y = 0.798, 0.035

        fig.text(text_x, text_y,
                 rf'$\omega_{{E \times\! B}}\,\,\mathrm{{{wexb_status}}}$',
                 size=8, color='#777')

    def write_zero_difference(fig, cd):
        """Writes text when two variables in difference plot are identical"""
        fig.text(0.525, 0.5, rf'Zero Difference', ha='center', size=10, color='#777')

    if cd.var_to_plot == 'nR8TOMSQZ':
        print('\t R8TOMSQZ Calls Per Point:', round(np.average(cd.Z[:, 1:]), 3))
    if cd.var_to_plot == 'nCubic':
        print('\t Cubic Calls Per Point:', round(np.average(cd.Z[:, 1:]), 3))
    if cd.var_to_plot == 'nWarning':
        print('\t Total Warnings:', int(np.sum(cd.Z)))
    if cd.var_to_plot == 'nError':
        print('\t Total Errors:', int(np.sum(cd.Z)))

    options = cd.options
    plot_options = cd.plot_options

    args_both, args_fill, args_line = {}, {}, {}
    set_colormaps(cd, args_both, args_fill, args_line)
    set_colorbar_extensions(cd, args_both)
    set_contour_levels(cd, args_fill, args_both)

    # Init figure and create contours
    fig, ax = init_figure(plt, plot_options)
    plot_contours()

    ax.yaxis.set_minor_formatter(NullFormatter())  # Forget why we want this

    # Set labels and other text
    plt.xlabel(cd.get_xlabel())
    plt.ylabel(cd.get_ylabel())
    plt.title(cd.get_title())

    if hasattr(cd, 'controls'):
        write_wexb_text(fig, cd)

    if hasattr(cd, 'isidentical') and cd.isidentical:
        write_zero_difference(fig, cd)

    # Save handling
    if plot_options.savedata or plot_options.savefig:
        savedir_base = f'{utils.get_plotting_contours_path()}\\{options.runid}'
        utils.create_directory(savedir_base)

        if plot_options.savedata:
            savedir = f'{savedir_base}\\data'
            utils.create_directory(savedir)
            savename = cd.get_savename()
            cd.save_to_csv(f'{savedir}\\{savename}')

        if plot_options.savefig:
            savedir = f'{savedir_base}\\figures'
            utils.create_directory(savedir)
            savename = cd.get_savename()
            fig.savefig(f'{savedir}\\{savename}')

    if plot_options.showfig:
        plt.show()

    fig.clear()
    plt.close(fig)  # might need to move outside of plotting loop


def run_plotting_loop(cd, vars_to_plot):
    """
    Creates PDF Plots of each variable in vars_to_plot

    The x-axis variable is pulled as var_to_scan from the saved options file

    Parameters:
    * cd (ContourData): Initialized contour data
    * vars_to_plot (list): List of output variables to plot
    """

    # Replace instances of 'var_to_scan' with actual scanned variable name
    vars_to_plot = [v if v != 'var_to_scan' else cd.options.var_to_scan for v in vars_to_plot]

    # Loop through all variables to plot
    for var_to_plot in vars_to_plot:

        print(f'- {cd.options.scan_num or cd.options.runid}, {var_to_plot}')

        if var_to_plot == 'time':
            _log.warning(f'\n\tNothing to plot when var_to_plot is time, continuing to next variable...')
            continue  # nothing to plot (this can happen when using 'var_to_scan')

        is_error = cd.set_z(var_to_plot)
        if is_error:
            continue

        make_contour_plot(cd)
