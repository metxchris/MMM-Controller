#!/usr/bin/python3

"""Creates contour plots of differences in data from different variable scans
"""

# Standard Packages
import sys; sys.path.insert(0, '../')
import logging

# 3rd Party Packages
import matplotlib.pyplot as plt

# Local Packages
import modules.options
import modules.utils as utils
from plotting.modules.plotstyles import PlotStyles, StyleType
import plotting.modules.contourdata as contourdata
import plotting.modules.makecontourplot as makecontourplot


_log = logging.getLogger(__name__)


class ContourPlotData:
    def __init__(self, runid, scannum, var_to_plot):
        self.runid: str = runid
        self.scannum: int = scannum
        self.var_to_plot: str = var_to_plot
        self.options = modules.options.Options().load(runid, scannum)


def plot_contour_difference(contour_list, plot_options):
    """
    Creates difference plots of variables specified in the contour data

    Parameters:
    * contour_list (list[ContourPlotData]): List of data to be plotted
    * plot_options (PlotOptions): Options for the contour plot
    """

    print(f'- {contour_list[0].var_to_plot}')

    cd = contourdata.ContourDataDiff(contour_list[0].options, contour_list[1].options, plot_options)

    is_error = cd.set_z(contour_list[0].var_to_plot)
    if is_error:
        return

    # cd.Z[cd.Z > 2.0] = 2.0
    # cd.Z[cd.Z < 0.5] = 0.5
    # cd.Z -= 1

    makecontourplot.make_contour_plot(cd)


def main(vars_to_plot, scan_data, plot_options):
    '''
    Verifies vars_to_plot, then plots a contour difference of all variables
    for each runid.

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * scan_data (dict): Dictionary of runid to list of scan numbers
    * plot_options (PlotOptions): Class of plotting options
    '''

    utils.init_logging()
    contourdata.verify_vars_to_plot(vars_to_plot)
    contour_list = []

    if plot_options.savefig:
        print(f'Figures will be saved in:\n\t{utils.get_plotting_contours_path()}')

    for runid, scan_nums in scan_data.items():
        print('\nInitializing data for', runid, scan_nums)

        if len(scan_nums) != 2:
            raise ValueError('Two scan numbers required per runid')

        for var_to_plot in vars_to_plot:
            contour_list.clear()
            for scan_num in scan_nums:
                contour_list.append(ContourPlotData(runid, scan_num, var_to_plot))

            if contour_list[0].var_to_plot == 'time' or contour_list[1].var_to_plot == 'time':
                _log.warning(
                    f'\n\tNothing to plot when var_to_plot is time'
                    f'\n\t{contour_list[0].scannum}: {contour_list[0].var_to_plot}'
                    f'\n\t{contour_list[0].scannum}: {contour_list[0].var_to_plot}'
                )
                continue

            if not contour_list[0].options.var_to_scan or not contour_list[1].options.var_to_scan:
                raise ValueError(
                    f'Variable scan is missing'
                    f'\n\t{contour_list[0].scannum}: {contour_list[0].options.var_to_scan}'
                    f'\n\t{contour_list[1].scannum}: {contour_list[1].options.var_to_scan}'
                )

            if contour_list[0].options.var_to_scan != contour_list[1].options.var_to_scan:
                raise ValueError(
                    f'Variable scan differs between each data set'
                    f'\n\t{contour_list[0].scannum}: {contour_list[0].options.var_to_scan}'
                    f'\n\t{contour_list[1].scannum}: {contour_list[1].options.var_to_scan}'
                )

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

    plot_options = contourdata.PlotOptions(
        xmin=0.2,
        xmax=0.8,
        ymin=0.2,
        ymax=0.6,
        showfig=1,
        savefig=0,
        savedata=0,
        saveappend='',
        difftype='pureratio',  # 'diff', 'ratio', 'absdiff', 'absratio', 'pureratio'
    )

    """
    Variables to Plot (List):
    * Uncomment the line you wish to use
    """

    vars_to_plot = ['xte', 'xti', 'xdi', 'xdz', 'xvt', 'xvp', 'vcz', 'vcp', 'vct']
    vars_to_plot = ['xteW20', 'xtiW20', 'xdiW20']
    # vars_to_plot = ['vci', 'vch', 'vce', 'vcz', 'vct', 'vcp']
    """
    Scan Data:
    * Uncomment the line you wish to use
        - keys (str): The runid of the scan
        - values (list of int): The scan_numbers to plot from
    """

    # vars_to_plot = ['fti', 'fte', 'fde', 'xdz', 'xvt', 'xvp']
    vars_to_plot = ['gmaEPM', 'xteEPM', 'xtiEPM', 'xdeEPM',]
    # vars_to_plot = ['gaveW20i', 'gaveW20e', 'kyrhosW20i', 'kyrhosW20e', 'kparaW20i', 'kparaW20e']

    # vars_to_plot = ['xkemtm', 'wexb', 'xtiEPM', 'xdeEPM',]
    # vars_to_plot = ['vcz', 'vct', 'vcp', 'xdeEPM',]
    # scan_data['138536A01'] = [1335, 1336]
    # scan_data['85610T01'] = [18, 19]

    vars_to_plot = ['xteETGM']
    # scan_data['120968A02'] = [25001, 25002]  # ETGM kyrhos optimization vs no optimization
    # scan_data['120982A09'] = [25001, 25002]  # ETGM kyrhos optimization vs no optimization
    # scan_data['129016A04'] = [25001, 25002]  # ETGM kyrhos optimization vs no optimization
    # scan_data['138536A01'] = [25001, 25002]  # ETGM kyrhos optimization vs no optimization
    # scan_data['101381T31'] = [25001, 25002]  # ETGM kyrhos optimization vs no optimization
    # scan_data['85126T02'] = [25001, 25002]  # ETGM kyrhos optimization vs no optimization
    # scan_data['15334T03'] = [25001, 25002]  # ETGM kyrhos optimization vs no optimization
    # scan_data['08505Z06'] = [25001, 25002]  # ETGM kyrhos optimization vs no optimization
    # scan_data['18696R06'] = [25001, 25002]  # ETGM kyrhos optimization vs no optimization
    # scan_data['84599T01'] = [25001, 25002]  # ETGM kyrhos optimization vs no optimization

    vars_to_plot = ['xtiW20'] # 39 = def, 40 = +, 41 = -, 42 = abs +, 43 = abs -
    # scan_data['120968A02'] = [25003, 25002]  # 7x7 no optimization vs 1000x1
    # scan_data['120982A09'] = [25003, 25002]  # 7x7 no optimization vs 1000x1
    # scan_data['129016A04'] = [25003, 25002]  # 7x7 no optimization vs 1000x1
    # scan_data['138536A01'] = [25004, 25002]  # 7x7 no optimization vs 1000x1
    # scan_data['101381T31'] = [25003, 25002]  # 7x7 no optimization vs 1000x1
    # scan_data['85126T02'] = [25003, 25002]   # 7x7 no optimization vs 1000x1
    # scan_data['15334T03'] = [25003, 25002]   # 7x7 no optimization vs 1000x1
    # scan_data['08505Z06'] = [25003, 25002]   # 7x7 no optimization vs 1000x1
    # scan_data['18696R06'] = [25003, 25002]   # 7x7 no optimization vs 1000x1
    # scan_data['84599T01'] = [25003, 25002]     # 7x7 no optimization vs 1000x1

    vars_to_plot = ['fte']
    # scan_data['120968A02'] = [25041, 25015]  # X vs 1000x1 default xte
    # scan_data['120982A09'] = [25017, 25015]    # X vs 1000x1 default xte
    # scan_data['129016A04'] = [25039, 25038]  # X vs 1000x1 default xte
    # scan_data['138536A01'] = [25025, 25013]  # X vs 1000x1 default xte
    # scan_data['138536A01'] = [25060, 25061]  # X vs 1000x1 default xte
    # scan_data['101381T31'] = [25011, 25012]  # X vs 1000x1 default xte
    # scan_data['85126T02'] = [25010, 25009]   # X vs 1000x1 default xte
    # scan_data['15334T03'] = [25008, 25007]   # X vs 1000x1 default xte
    # scan_data['08505Z06'] = [25008, 25007]   # X vs 1000x1 default xte
    # scan_data['18696R06'] = [25008, 25007]   # X vs 1000x1 default xte
    # scan_data['84599T01'] = [25008, 25007]   # X vs 1000x1 default xte

    vars_to_plot = ['nR8TOMSQZ', 'xtiW20', 'xteW20']
    # vars_to_plot = ['gmaW20i', 'gmaW20e']
    # vars_to_plot = ['xtiW20']
    # scan_data['120968A02'] = [26205, 26204]  # X vs 1000x1 default xte
    # scan_data['120982A09'] = [26101, 26103]    # X vs 1000x1 default xte
    # scan_data['129016A04'] = [26101, 26103]  # X vs 1000x1 default xte
    # scan_data['138536A01'] = [26101, 26103]  # X vs 1000x1 default xte
    # scan_data['101381T31'] = [26101, 26102]  # X vs 1000x1 default xte
    # scan_data['132498J05'] = [26101, 26201]  # X vs 1000x1 default xte
    # scan_data['85126T02']  = [26101, 26102]   # X vs 1000x1 default xte
    # scan_data['90328T01']  = [26101, 26102]   # X vs 1000x1 default xte
    # scan_data['15334T03']  = [26101, 26102]   # X vs 1000x1 default xte
    # scan_data['18399T05']  = [26101, 26102]   # X vs 1000x1 default xte
    # scan_data['24899R05']  = [26207, 26206]   # X vs 1000x1 default xte
    # scan_data['29976U69']  = [26101, 26102]   # X vs 1000x1 default xte
    # scan_data['38265R80']  = [26101, 26102]   # X vs 1000x1 default xte
    # scan_data['80200A13']  = [26101, 26102]   # X vs 1000x1 default xte
    # scan_data['84599T01']  = [26101, 26102]   # X vs 1000x1 default xte
    # scan_data['87261T01']  = [26101, 26102]   # X vs 1000x1 default xte


    ## Akima derivative vs normal derivative
    vars_to_plot = ['gte']
    # scan_data['138536A01'] = [26205, 26204]  # Traditional vs Akima derivative
    # scan_data['129016A04'] = [26202, 26201]  # Traditional vs Akima derivative

    # vars_to_plot = ['gmaDBM','gmaETGM']
    # scan_data['118341T54'] = [536, 531]  #
    # plot_options.title, vars_to_plot = r'$|\chi_{\mathrm{i, opt}} - \chi_{\mathrm{i, def}}|$ (m$^2$/s) ', ['xtiW20']
    # plot_options.title, vars_to_plot = r'$|\chi_{\mathrm{e, opt}} - \chi_{\mathrm{e, def}}|$ (m$^2$/s) ', ['xteW20']
    # plot_options.title, vars_to_plot = r'$|\chi_{\mathrm{n, opt}} - \chi_{\mathrm{n, def}}|$ (m$^2$/s) ', ['xdeW20']

    # scan_data['120968A02'] = [5,6]
    # scan_data['138536A01'] = [i for i in range(1716, 1738 + 1)]
    # scan_data['138536A01'] = [i for i in [*range(1716, 1738 + 1), *range(1756, 1763 + 1)]]

    main(vars_to_plot, scan_data, plot_options)
