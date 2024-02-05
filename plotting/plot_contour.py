#!/usr/bin/python3

"""Creates contour plots of all data from a variable scan

This module is designed to automatically save contour plots and their
associated data from a variable scan.
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


def main(vars_to_plot, scan_data, plot_options):
    '''
    Verifies vars_to_plot, then runs the plotting loop for each var_to_plot, runid, and scan_num

    Parameters:
    * vars_to_plot (list): List of output variables to plot
    * scan_data (dict): Dictionary of runid to list of scan numbers
    * plot_options (PlotOptions): Object containing options for the contour plot (Optional)
    '''

    utils.init_logging()
    contourdata.verify_vars_to_plot(vars_to_plot)
    options = modules.options.Options()

    if plot_options.savefig:
        print(f'Files will be saved in:\n\t{utils.get_plotting_contours_path()}')

    plt.figure()  # Instantiate the figure

    for runid, scan_nums in scan_data.items():
        for scan_num in scan_nums:
            options.load(runid, scan_num)
            if options.var_to_scan:
                print(f'\nInitializing data for {runid}, scan {scan_num}, {options.var_to_scan}...')
                cd = contourdata.ContourDataMMM(options, plot_options)
                makecontourplot.run_plotting_loop(cd, vars_to_plot)
            else:
                _log.error(f'\n\tNo variable scan detected for {runid}, scan {scan_num}\n')


# Run this file directly to plot scanned variable profiles from previously created scanned data
if __name__ == '__main__':
    scan_data = {}

    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.AIP2,
    )

    plot_options = contourdata.PlotOptions(
        # xmin=0.1,
        # xmax=0.9,
        # ymin=0.5,
        # ymax=0.6,
        savefig=0,
        savedata=0,
        smoothing=0,
        showfig=1,
        saveappend='',
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
    })

    """
    Variables to Plot (List):
    * Uncomment the line you wish to use
    """
    # vars_to_plot = ['var_to_scan']
    # vars_to_plot = ['xti', 'xte', 'xde', 'xdz', 'xvt', 'xvp', 'vcz', 'vcp', 'vct']
    # vars_to_plot = ['xti', 'xte', 'xde']
    # vars_to_plot = ['xteDBM', 'xtiDBM', 'xdiDBM', 'gmaDBM', 'omgDBM', 'kyrhosDBM', 'phi2DBM', 'Apara2DBM', 'gaveDBM', 'satDBM']
    # vars_to_plot = ['gmaW20i', 'gmaW20e']
    # vars_to_plot = OutputVariables().get_all_output_vars()
    # vars_to_plot = OutputVariables().get_etgm_vars()

    """
    Scan Data:
    * Uncomment the line you wish to use
        - keys (str): The runid of the scan
        - values (list of int): The scan_numbers to plot from
    """

    vars_to_plot = ['satETGM']
    ## 20001: 7x7 w/ optimization
    ## 20002: 7x7 no opt
    ## 20003: 1000x1
    # scan_data['120968A02'] = [25044];
    # scan_data['120982A09'] = [25018];  
    # scan_data['129016A04'] = [25040];  
    # scan_data['138536A01'] = [25010];
    # scan_data['101381T31'] = [25008];
    # scan_data['85126T02'] = [25010];
    # scan_data['15334T03'] = [25008];
    # scan_data['08505Z06'] = [25008];
    # scan_data['18696R06'] = [25008]; 
    # scan_data['84599T01'] = [25008]; 

    # scan_data['129016A04'] = [25029]; # kyrhos scan

    # vars_to_plot = ['satETGM','satDBM',] # 20001 - 7x7 scans w/ opt
    vars_to_plot = ['xteMTM', 'xtiW20', 'xteW20'] # 20001 - 7x7 scans w/ opt
    vars_to_plot = ['gmaW20i', 'gmaW20e'] # 20001 - 7x7 scans w/ opt
    # vars_to_plot = ['nR8TOMSQZ'] # 20001 - 7x7 scans w/ opt
    # vars_to_plot = ['xteETGM'] # 20001 - 7x7 scans w/ opt
    # scan_data['120968A02'] = [26207];  # 7x7 scans, w/ optimization (29.3 scans)
    # scan_data['120982A09'] = [26103];  # 7x7 scans, w/ optimization (32.2 scans)
    # scan_data['129016A04'] = [26103];  # 7x7 scans, w/ optimization (28.5 scans)
    # scan_data['138536A01'] = [26202];  # 7x7 scans, w/ optimization (29.2 scans)
    # scan_data['101381T31'] = [26102];  # 7x7 scans, w/ optimization (34.0 scans)
    # scan_data['132498J05'] = [26101];  # 7x7 scans, w/ optimization (34.0 scans)
    # scan_data['85126T02']  = [26201];  # 7x7 scans, w/ optimization  (30.1 scans)
    # scan_data['90328T01']  = [26101];  # 7x7 scans, w/ optimization  (32.4 scans)
    # scan_data['15334T03']  = [26102];  # 7x7 scans, w/ optimization  (26.3 scans)
    # scan_data['18399T05']  = [26100];  # 7x7 scans, w/ optimization  (21.7 scans)
    # scan_data['24899R05']  = [26206];  # 7x7 scans, w/ optimization  (31.9 scans)
    # scan_data['29976U69']  = [26100];  # 7x7 scans, w/ optimization  (30.1 scans)
    # scan_data['38265R80']  = [26100];  # 7x7 scans, w/ optimization  (32.4 scans)
    # scan_data['80200A13']  = [26101];  # 7x7 scans, w/ optimization  (26.3 scans)
    # scan_data['84599T01']  = [26100];  # 7x7 scans, w/ optimization  (21.7 scans)
    # scan_data['87261T01']  = [26100];  # 7x7 scans, w/ optimization  (31.9 scans)

    vars_to_plot = ['gte'] # 20001 - 7x7 scans w/ opt
    scan_data['15334T03'] = [26224];  # 7x7 scans, w/ optimization (34.0 scans)

    ## Examples using ranges
    # scan_data['138536A01'] = [i for i in range(1716, 1738 + 1)]
    # scan_data['138536A01'] = [i for i in [*range(1716, 1738 + 1), *range(1756, 1763 + 1)]]

    ## Run program
    main(vars_to_plot, scan_data, plot_options)
