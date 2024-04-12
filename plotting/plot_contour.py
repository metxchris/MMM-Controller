#!/usr/bin/python3

"""Creates contour plots of all data from a variable scan

This module is designed to automatically save contour plots and their
associated data from a variable scan.
"""

# Standard Packages
import sys; sys.path.insert(0, '../')
import logging
import warnings

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

    for i, (runid, scan_nums) in enumerate(scan_data.items()):
        for scan_num in scan_nums:
            options.load(runid, scan_num)
            if options.var_to_scan:
                print(f'\nInitializing data for {runid}, scan {scan_num}, {options.var_to_scan} ({i + 1}/{len(scan_data)})...')
                cd = contourdata.ContourDataMMM(options, plot_options)
                makecontourplot.run_plotting_loop(cd, vars_to_plot)
            else:
                _log.error(f'\n\tNo variable scan detected for {runid}, scan {scan_num}\n')


# Run this file directly to plot scanned variable profiles from previously created scanned data
if __name__ == '__main__':
    warnings.filterwarnings("ignore")

    scan_data = {}
    vars_to_plot = []

    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.AIP2,
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
        'figure.dpi': 150,
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

    # vars_to_plot = ['satETGM']
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
    # vars_to_plot = ['xteMTM', 'xtiW20', 'xteW20'] # 20001 - 7x7 scans w/ opt
    # vars_to_plot = ['gmaW20i', 'gmaW20e'] # 20001 - 7x7 scans w/ opt
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

    # vars_to_plot = ['xtiW20'] # 20001 - 7x7 scans w/ opt
    # scan_data['118341T54'] = [26301];  # 7x7 scans, w/ optimization (34.0 scans)

    ## Examples using ranges
    # scan_data['138536A01'] = [i for i in range(1716, 1738 + 1)]
    # scan_data['138536A01'] = [i for i in [*range(1716, 1738 + 1), *range(1756, 1763 + 1)]]

    # scan_data['141716A80']  = [26202];  # 1 mode, scan both ways

    # vars_to_plot = ['gxi','gaveW20i', 'gaveW20e']
    # vars_to_plot = ['fti', 'fte', 'fde', 'fdz', 'fvp', 'fvt']
    # vars_to_plot = ['xtiW20', 'xteW20', 'xdeW20', 'xdz', 'xvt', 'xvp']
    # vars_to_plot += ['gti', 'gte', 'gne', 'gnz', 'gvt', 'gvp']
    # vars_to_plot += ['nR8TOMSQZ', 'nWarning', 'nError']
    # vars_to_plot = ['gma0W20', 'gmagW20', 'gmaW20']
    # vars_to_plot = ['gaveETGM', 'gmaW20e', 'gmaW20']
    # vars_to_plot = ['nR8TOMSQZ']
    # vars_to_plot += ['tvt', 'xvt', 'vct', 'gvtor', 'fvt']
    # vars_to_plot = ['gaveW20i', 'gaveW20e']
    # num = 26325  # 26322, replacement rule on initial solutions
    # num = 26324  # 26314, only ion direction guesses
    # num = 26323  # 26314, but no mode replacement
    # num = 26322  # 26314, but converge in both directions
    # num = 26314  # 26302, only conv initial modes + gave kpara fix
    # num = 26311  # 26302, with max mode from match array fixed
    # num = 26307  # 26305, but coded without overwriting ion vars
    # num = 26305  # 26303, but searching in both directions
    # num = 26303  # 26302, but using single most unstable mode
    # num = 26302  # ky_i = ky_e, wexb_i = wexb_e
    # num = 26301  # W20 9.0.10
    
    num = 27083  # #114 w20 no negative chi
    # num = 27073  # #114 w20 only
    # num = 27070  # #114 Channel swap, var rename
    # num = 27030  # Major w20 fixes, + alp min zepsqrt
    # num = 27027  # Major w20 fixes, + suppressed kap1, ne/nh Curr, guess fix
    # num = 27026  # Major w20 fixes, + removing kap1 max
    # num = 27025  # Major w20 fixes, + zflh, zflz, geometry
    # num = 27024  # Major w20 fixes, + epsilon Gave
    # num = 27023  # Major w20 fixes, Fixing XI(3) term
    # num = 27021  # Major w20 fixes
    # num = 27020  # #103 (cleanup of w20 calculations)
    # num = 27015  # #102, QZ zepsqrt, Gave min = zepslon
    # num = 27014  # Commit 102, Gave min = zepslon
    # num = 27013  # Commit 101
    # num = 27011  # pre w20 fixes, Commit 102
    # num = 26386  # convstrat 6 rolling guess denominator
    # num = 26385  # convstrat 5 cleanup, max 99, stuck 30
    # num = 26382  # convstrat 4 counting opp gma, omgr
    # num = 26381  # convstrat 3 counting mode stuck
    # num = 26380  # convstrat 2 based on rel dif > 1E-1
    # num = 26376  # convergence strategy OG
    # num = 26375  # NEW Baseline (stricter duplicate modes)
    # num = 26371  # conv fix 6 OPT#7 (e, i mode replacement)
    # num = 26370  # conv fix 6 OPT#7 (no mode replacement)
    # num = 26357  # conv fix 6 OPT#7
    # num = 26356  # conv fix 6 OPT#6
    # num = 26355  # conv fix 6 OPT#5
    # num = 26354  # conv fix 6 OPT#4
    # num = 26352  # 9.0.10 OPT
    # num = 26351  # 9.0.10 NO OPT
    # num = 26348  # conv fix 6 OPT#3
    # num = 26347  # conv fix 6 OPT#2
    # num = 26343  # i = e OPT
    # num = 26342  # i = e NO OPT
    # num = 26341  # conv fix 6 OPT
    # num = 26340  # conv fix 6 NO OPT
    # num = 61
    vars_to_plot = ['nR8TOMSQZ', 'nWarning', 'nError']
    vars_to_plot = ['nCubic']
    # vars_to_plot = ['gmaW20', 'gma0W20', 'gmagW20', ]
    # vars_to_plot += ['omgW20', 'omg0W20', 'omggW20', ]
    # vars_to_plot = ['gmagW20', 'gma0W20']
    # vars_to_plot = ['gmaDBM']
    # vars_to_plot = ['fde', 'xde', 'vde']
    # vars_to_plot = ['fvt', 'xvt', 'vvt']
    # vars_to_plot += ['fte', 'xte', 'vte']
    # vars_to_plot += ['fti', 'xti', 'vti']
    # vars_to_plot = ['fdz', 'xdz', 'vdz']
    # vars_to_plot = ['dvtor_dr', 'dwtor_dr']
    # vars_to_plot = ['fti', 'fte', 'fde']
    # vars_to_plot = ['te', 'ti']
    # vars_to_plot = ['gne', 'fde', 'xde']
    # vars_to_plot = ['fti', 'fte', 'fne', 'fnz', 'fvt', 'fvp']
    # vars_to_plot = ['fti', 'fte', 'fde', 'fdz', 'fvt', 'fvp']
    # vars_to_plot = ['xti', 'xte', 'xne', 'xnz', 'xvt', 'xvp']
    # vars_to_plot = ['fne', 'xne', 'vne']
    # vars_to_plot = ['gmaETGM', 'phi2ETGM', 'Apara2ETGM', 'xteETGM', 'xte2ETGM']
    # vars_to_plot = ['gmaW20', 'xtiW20', 'xteW20', 'xneW20', 'xnz', 'xvt', 'xvp']
    # vars_to_plot = ['gmaETGM', 'satETGM', 'phiETGM', 'AparaETGM', 'xteETGM', 'xte2ETGM',]
    # vars_to_plot = ['xteETGM', 'xte2ETGM',]
    vars_to_plot = ['fte', 'xte', 'xteETGM',]
    # vars_to_plot = ['kyrhosMTM',]

    scan_data['179415P02'] = [27117]
    # scan_data['138536A01'] = [10027]

    # vars_to_plot = ['gxi',]
    # scan_data['120968A02'] = [42]

    # 48: no split, no convert
    # 49: split, no convert
    # 56: no split, convert
    # 57: split, convert
    # 58: no split, no convert, no neg
    # 63: solvers
    # 64: No solvers
    # 73: zflh no ah
    # 74: zflh with ah
    # 79, 77: no limit small grad
    # 80, 78: limit small grad
    # scan_data['85122L01'] = [11]
    # scan_data['129016A04'] = [127]
    # scan_data['138536A01'] = [84]
    # scan_data['183743H01'] = [8]
    # scan_data['120968A02'] = [55]
    # scan_data['90949R01'] = [8]

    plot_options = contourdata.PlotOptions(
        # xmin=0.60,
        # xmax=0.75,
        ymin=0.3,
        # ymax=500,
        showfig=1,
        savefig=0,
        savedata=0,
        smoothing=0,
        saveappend='',
    )

    ## Run program
    main(vars_to_plot, scan_data, plot_options)
