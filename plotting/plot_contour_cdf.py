#!/usr/bin/python3

"""Creates contour plots of all data from a variable scan

This module is designed to automatically save contour plots and their
associated data from a CDF.

TODO:
* This module could also be adapted to plot variables directly from a CDF
"""

# Standard Packages
import sys; sys.path.insert(0, '../')
import logging

# 3rd Party Packages
import matplotlib.pyplot as plt
import numpy as np

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
    * plot_options (PlotOptions): Object containing options for the contour plot
    '''

    utils.init_logging()
    contourdata.verify_vars_to_plot(vars_to_plot)
    options = modules.options.Options()
    options.set(
        adjustment_name='time',
        scan_range=np.linspace(
            start=plot_options.time_min,
            stop=plot_options.time_max,
            num=plot_options.time_count,
        ),
    )

    if plot_options.savefig:
        print(f'Files will be saved in:\n\t{utils.get_plotting_contours_path()}')

    plt.figure()  # Instantiate the figure

    for runid in scan_data:
        options.set(runid=runid)
        print(f'\nInitializing data for {runid}...')
        cd = contourdata.ContourDataCDF(options, plot_options)
        makecontourplot.run_plotting_loop(cd, vars_to_plot)


# Run this file directly to plot scanned variable profiles from previously created scanned data
if __name__ == '__main__':
    scan_data = []

    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.AIP2,
    )

    plot_options = contourdata.PlotOptions(
        # xmin=0.05,
        # xmax=0.95,
        # ymin=0.5,
        # ymax=0.55,
        time_min=0,
        time_max=1,
        time_count=100,
        raw=0,
        savefig=0,
        savedata=0,
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
    vars_to_plot = ['ti', 'te', 'ne', 'ni']
    # vars_to_plot = ['gne', 'gte',]
    vars_to_plot = ['fte', 'gte',]
    # vars_to_plot = ['gmaW20ii', 'gmaW20ee', 'gmaW20ie']
    # vars_to_plot = ['gmaW20i', 'gmaW20e']
    # vars_to_plot = ['xkemmm', 'xkiw20', 'xdew20', 'xkew20', 'xkemtm']
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
    # vars_to_plot = ['xkiw20', 'xkew20','xkemtm', ]
    # vars_to_plot = ['gti', 'gte','ti','te', ]
    # vars_to_plot = ['kyrsetg', 'xke', 'xkemmm', 'gte', 'te', ]
    # vars_to_plot = ['agxi_1', 'agxi_1b','agxi2_1',]
    vars_to_plot = ['ti', 'te', ]
    # vars_to_plot = ['xkemmm', 'gte', 'te', ]
    # vars_to_plot = ['tz', 'fki', 'xki', 'xkimmm', 'gti', 'ti', ]
    # vars_to_plot = ['wexbsmod', 'wexbsv2', 'te', 'ti', 'wexb',]
    # vars_to_plot = ['vphimmm', 'gamma1w20', 'gamma2w20', 'omega1w20', 'omega2w20',]
    # vars_to_plot = ['xkemmm', 'xkimmm', 'xkdmmm', 'xkzmmm',  'xppmmm',  'xptmmm', ]
    # vars_to_plot = ['xkemmm', 'xkemtm', 'xkeetgm', 'xkew20', 'xkedrbm']
    # vars_to_plot = ['xkdmmm', 'xkemtm', 'xkeetgm']
    # vars_to_plot = ['vphimmm', 'kyrsdbm', 'kyrsetg', 'kyrsmtm',  'xkeetg',  'xptmmm', ]
    # vars_to_plot = [ 'te', 'xkeetg']
    # vars_to_plot = [ 'gamma1w20', 'xkeetg']
    # vars_to_plot = [ 'omega']
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
    # scan_data.append('129017W02') # MMM 8.2.1
    # scan_data.append('129017W03') # MMM 9.0.7
    # scan_data.append('129016W02') # MMM 9.0.7 W20+MTM Failure
    # scan_data.append('129016W03') # MMM 9.0.7 ETGM +chi
    # scan_data.append('129016W11') # MMM 9.0.7 ETGM +-chi
    # scan_data.append('129016W15') # MMM 9.0.7 ETGM +-chi, 0 PT
    # scan_data.append('129016W12') # MMM 9.0.7 Default = pm chi
    # scan_data.append('129016W16') # MMM 9.0.7 pphi
    # scan_data.append('129016W21') # MMM 9.0.9 g>=1
    # scan_data.append('129016W23') # MMM 9.0.9 g>=1
    # scan_data.append('129016W24') # Default

    # scan_data.append('147634T61') # Default
    # scan_data.append('113944W20') # Failed EAST run
    # scan_data.append('113944W15') # Good old EAST run

    # scan_data.append('129016W12') # Testing

    vars_to_plot = ['vndnc_e', 'dvtor_dr', 'dwtor_dr']
    vars_to_plot = ['wexb',]
    # scan_data.append('129016A04')  # Testing
    scan_data.append('85122L01')  # Testing

    # vars_to_plot = ['xkidrbm', 'xkddrbm', 'xkedrbm', 'xkhdrbm']
    # scan_data.append('129016Q71') # DBM, 0.001 cal

    main(vars_to_plot, scan_data, plot_options)
