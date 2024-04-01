#!/usr/bin/python3

"""Creates contour plots of differences in data from different variable scans
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

    for i, (runid, scan_nums) in enumerate(scan_data.items()):
        print(f'\nInitializing data for {runid} {scan_nums} ({i + 1}/{len(scan_data)})')

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
    })

    """
    Variables to Plot (List):
    * Uncomment the line you wish to use
    """

    # vars_to_plot = ['xte', 'xti', 'xdi', 'xdz', 'xvt', 'xvp', 'vcz', 'vcp', 'vct']
    # vars_to_plot = ['xteW20', 'xtiW20', 'xdiW20']
    # vars_to_plot = ['vci', 'vch', 'vce', 'vcz', 'vct', 'vcp']
    """
    Scan Data:
    * Uncomment the line you wish to use
        - keys (str): The runid of the scan
        - values (list of int): The scan_numbers to plot from
    """

    # vars_to_plot = ['gmaEPM', 'xteEPM', 'xtiEPM', 'xdeEPM',]
    # vars_to_plot = ['gaveW20i', 'gaveW20e', 'kyrhosW20i', 'kyrhosW20e', 'kparaW20i', 'kparaW20e']
    # vars_to_plot = ['xkemtm', 'wexb', 'xtiEPM', 'xdeEPM',]
    # vars_to_plot = ['vcz', 'vct', 'vcp', 'xdeEPM',]
    # scan_data['138536A01'] = [1335, 1336]
    # scan_data['85610T01'] = [18, 19]

    # vars_to_plot = ['xteETGM']
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

    # vars_to_plot = ['xtiW20'] # 39 = def, 40 = +, 41 = -, 42 = abs +, 43 = abs -
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

    # vars_to_plot = ['fte']
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

    # vars_to_plot = ['nR8TOMSQZ', 'xtiW20', 'xteW20']
    # vars_to_plot = ['gmaW20i', 'gmaW20e']
    # vars_to_plot = ['xteETGM']
    # scan_data['120968A02'] = [26205, 26204]  # X vs 1000x1 default xte
    # scan_data['120982A09'] = [26101, 26103]    # X vs 1000x1 default xte
    # scan_data['129016A04'] = [26101, 26103]  # X vs 1000x1 default xte
    # scan_data['138536A01'] = [26215, 26212]  # X vs 1000x1 default xte
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

    # scan_data['141716A80'] = [26202, 26302]

    # vars_to_plot += ['xtiW20', 'xteW20', 'xdeW20', 'xdz', 'xvt', 'xvp']
    # vars_to_plot += ['nR8TOMSQZ']
    # vars_to_plot = ['gmaW20e', 'gmaW20i', 'kparaW20i', 'kparaW20e', 'gaveW20i', 'gaveW20e', 'omgW20i', 'omgW20e']
    # vars_to_plot += ['nR8TOMSQZ', 'nWarning', 'nError']
    # vars_to_plot += ['gmaW20','gma0W20', 'gmagW20']
    # vars_to_plot = ['xtiW20', 'fte', 'fvt']

    n_new = 27041  # #108 Limiting diff due to small grad
    # n_new = 27040  # #108 Default
    # n_new = 27032  # Major w20 fixes, + R_curv_min 0.1
    # n_new = 27031  # Major w20 fixes, + R_curv_min 0.01
    # n_new = 27030  # Major w20 fixes, + alp min zepsqrt
    # n_new = 27027  # Major w20 fixes, + suppressed kap1, ne/nh Curr, guess fix
    # n_new = 27026  # Major w20 fixes, + removing kap1 max
    # n_new = 27025  # Major w20 fixes, + zflh, zflz, geometry
    # n_new = 27024  # w20 Major fixes, + epsilon Gave
    # n_new = 27023  # w20 Major fixes, fixing XI(3) term
    # n_new = 27022  # w20 Major fixes, no kap1 max
    # n_new = 27021  # w20 Major fixes
    # n_new = 27015  # #102, QZ zepsqrt, Gave min = zepslon
    # n_new = 27014  # #102, Gave min = zepsqrt
    # n_new = 27011  # #102
    # n_new = 26387  # convstrat 6 rel_diff 1E-3
    # n_new = 26386  # convstrat 6 rolling guess denominator
    # n_new = 26385  # convstrat 5 cleanup, max 99, stuck 30
    # n_new = 26382  # convstrat 4 counting opp gma, omgr
    # n_new = 26381  # convstrat 3 counting mode stuck
    # n_new = 26380  # convstrat 2 based on rel dif > 1E-1
    # n_new = 26376  # New convergence strategy OG
    # n_new = 26375  # NEW Baseline (stricter duplicate modes)
    # n_new = 26371  # conv fix Opt #7  (elc + ion mode replacement)
    # n_new = 26357  # conv fix Opt #7  (i > 4 .AND. im_omg_max > 0)
    # n_new = 26356  # conv fix Opt #6  (i > 2 .AND. im_omg_max > 0)
    # n_new = 26355  # conv fix Opt #5  NSTX OPT
    # n_new = 26354  # conv fix Opt #4
    # n_new = 26352  # 9.0.10 OPT
    # n_new = 26348  # conv fix Opt #3
    # n_new = 26347  # conv fix Opt #2
    # n_new = 26343  # i=e (OPT)
    # n_new = 26341  # conv fix 6 (OPT)
    # n_new = 26340  # conv fix 6 (NO OPT)
    # n_new = 26325  # conv fix 5 (NO OPT)
    # n_new = 26324  # conv fix 4 (NO OPT)
    # n_new = 26323  # conv fix 3 (NO OPT)
    # n_new = 26322  # conv fix 2 (NO OPT)
    # n_new = 26321  # conv fix 2 (BASELINE OPT)
    # n_new = 26314  # conv fix 1

    n_old = 27110  # #122 v9.1.0 w20 
    # n_old = 27040  # #108 Default
    # n_old = 27027  # Major w20 fixes, + suppressed kap1, ne/nh Curr, guess fix
    # n_old = 27025  # Major w20 fixes, + zflh, zflz, geometry
    # n_old = 27020  # #103 (cleanup of w20 equations)
    # n_old = 27014  # #102, Gave min = zepslon
    # n_old = 27013  # #101
    # n_old = 27011  # #102
    # n_old = 26385  # convstrat 5
    # n_old = 26375  # NEW Baseline (stricter duplicate modes)
    # n_old = 26370  # conv fix opt #7 (no replacement)
    # n_old = 26358  # conv fix opt #7 (no elc freq guess)
    # n_old = 26355  # conv fix Opt #5  NSTX OPT
    # n_old = 26348  # conv fix Opt #3
    # n_old = 26352  # 9.0.10 OPT
    # n_old = 26351  # 9.0.10 NO OPT
    # n_old = 26343  # i=e OPT
    # n_old = 26342  # i=e NO OPT
    # n_old = 26340  # conv fix 6 NO OPT
    # n_old = 26302  # 9.0.10 with i = e


    # scan_data['80200A13'] = [225, 224] #127 Base
    # scan_data['120968A02'] = [47, 45] #16 Base
    # scan_data['150840T02'] = [13, 26380]
    # scan_data['138536A01'] = [35, 36]
    # scan_data['138536A01'] = [37, 38]
    # 11 default no wexb
    # 12 default with wexb
    # scan_data['138536A01'] = [10, 11] # 11 default no wexb
    
    # 14: zflh change
    # 14: rmaj0 to rmaj
    # scan_data['138536A01'] = [15, 14] # 12 default with wexb

    ## OLD
    # vars_to_plot = ['fti', 'fte', 'fde', 'fdz', 'fvt', 'fvp']
    # vars_to_plot = ['xti', 'xte', 'xde', 'xdz', 'xvt', 'xvp']

    # vars_to_plot = ['fde', 'xde', 'vde']
    # vars_to_plot = ['fdz', 'xdz', 'vdz']
    # vars_to_plot = ['fti', 'xti', 'vti']
    # vars_to_plot = ['fte', 'xte', 'vte']
    # vars_to_plot = ['fvt', 'xvt', 'vvt']
    vars_to_plot = ['fti', 'fte', 'fne', 'fnz', 'fvt', 'fvp']
    # vars_to_plot = ['xti', 'xte', 'xne', 'xnz', 'xvt', 'xvp']
    # vars_to_plot = ['fne', 'xne',]
    # vars_to_plot = ['gmaW20', 'xtiW20', 'xteW20', 'xneW20', 'xnz', 'xvt', 'xvp']
    # vars_to_plot = ['gmaETGM', 'phiETGM', 'AparaETGM', 'xteETGM', 'xte2ETGM']
    # scan_data['138536A01'] = [118, 116]
    
    # vars_to_plot = ['gxi',]
    scan_data['129016A04'] = [2, 27117]

    # scan_data['90949R01'] = [11, 27064]


    ## Akima derivative vs normal derivative
    # vars_to_plot = ['gte']
    # scan_data['138536A01'] = [26205, 26204]  # Traditional vs Akima derivative
    # scan_data['129016A04'] = [26202, 26201]  # Traditional vs Akima derivative

    # scan_data['120968A02'] = [5,6]
    # scan_data['138536A01'] = [i for i in range(1716, 1738 + 1)]
    # scan_data['138536A01'] = [i for i in [*range(1716, 1738 + 1), *range(1756, 1763 + 1)]]

    # vars_to_plot = ['xtiW20', 'xteW20', 'xvt',]
    # vars_to_plot = ['fti', 'fte', 'fvt',]
    # scan_data['147634T61'] = [45, 26314]  # 26314

    plot_options = contourdata.PlotOptions(
        # xmin=0.2,
        # xmax=0.8,
        # ymin=0.2,
        # ymax=0.6,
        zmax_diff=100,
        zmin_diff=-100,
        showfig=1,
        savefig=0,
        savedata=0,
        plotidentical=0,
        saveappend='',
        difftype='diff',  # 'diff', 'rel', 'absdiff', 'absrel', 'ratio'
    )

    main(vars_to_plot, scan_data, plot_options)
