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
    # OLD
    # vars_to_plot += ['fti', 'fte', 'fde', 'fdz', 'fvt', 'fvp']
    # vars_to_plot += ['xti', 'xte', 'xde', 'xdz', 'xvp', 'xvt']

    # vars_to_plot += ['xtiW20', 'xteW20', 'xdeW20', 'xdz', 'xvt', 'xvp']
    vars_to_plot += ['fti', 'fte', 'fne', 'fnz', 'fvt', 'fvp']
    # vars_to_plot += ['xti', 'xte', 'xne', 'xnz', 'xvp', 'xvt']
    # vars_to_plot += ['fte']
    # vars_to_plot += ['fti', 'fte', 'fne']
    # vars_to_plot += ['nR8TOMSQZ']
    # vars_to_plot += ['nCubic']
    # vars_to_plot = ['gmaW20e', 'gmaW20i', 'kparaW20i', 'kparaW20e', 'gaveW20i', 'gaveW20e', 'omgW20i', 'omgW20e']
    vars_to_plot += ['nR8TOMSQZ', 'nWarning', 'nError']
    # vars_to_plot += ['gmaW20','gma0W20', 'gmagW20']
    # vars_to_plot = ['xtiW20', 'fte', 'fvt']
    
    n_new = 27124  # #126 9.1.1 w20 zggev
    # n_new = 27123  # #126 9.1.1 mtm default
    # n_new = 27122  # #126 9.1.1 etgm default
    # n_new = 27121  # #126 9.1.1 dbm default
    # n_new = 27120  # #126 9.1.1 w20
    # n_new = 27101  # #117 dbm shat_e
    # n_new = 27100  # #117 dbm default vei (no etanc)
    # n_new = 27099  # #117 dbm default xne, dbm nh_ne, ti/te no max
    # n_new = 27098  # #117 dbm nh_ne, ti/te no max
    # n_new = 27097  # #117 dbm ti/te max 3
    # n_new = 27095  # #117 dbm ti/te no max
    # n_new = 27094  # #117 dbm ti/te in alphaMHD
    # n_new = 27092  # #117 w20 derived gxi
    # n_new = 27091  # #117 w20 kparaN in disp change
    # n_new = 27084  # #115 w20 eigensolver fix
    # n_new = 27083  # #114 w20 no negative chi
    # n_new = 27065  # #111 Newton-Rahpson
    # n_new = 27064  # #111 w20 min gma diffusivity zepsqrt
    # n_new = 27063  # #111 kap1 no max
    # n_new = 27062  # #111 w20 velocity min zepslon
    # n_new = 27061  # #111 min eps zepslon
    # n_new = 27060  # #111 min grad zepslon + limit small grad
    # n_new = 27059  # #111 min grad zepslon
    # n_new = 27058  # #111 w20 no max ti/te
    # n_new = 27057  # #111 W20 XI(4) Term
    # n_new = 27056  # #111 W20 max ti/te = 10
    # n_new = 27055  # #111 W20 FIH min zepslon
    # n_new = 27053  # #111 DBM Default no sum n + T
    # n_new = 27052  # #111 DBM Default ti/te normal, max 5
    # n_new = 27050  # #109 MTM 400
    # n_new = 27049  # #109 MTM 100
    # n_new = 27047  # #109 MTM 200
    # n_new = 27045  # #108 increase pol to 1E1
    # n_new = 27041  # #108 Limiting diff due to small grad
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

    n_old = 27120  # #126 9.1.1 w20
    # n_old = 27117  # #122 v9.1.0 mtm
    # n_old = 27116  # #122 v9.1.0 etgm horton 2
    # n_old = 27115  # #122 v9.1.0 etgm horton 1
    # n_old = 27114  # #122 v9.1.0 etgm alternate
    # n_old = 27113  # #122 v9.1.0 etgm default
    # n_old = 27112  # #122 v9.1.0 dbm alternate
    # n_old = 27111  # #122 v9.1.0 dbm default
    # n_old = 27110  # #122 v9.1.0 w20 
    # n_old = 27099  # #117 default xne, dbm nh_ne, ti/te no max
    # n_old = 27095  # #117 dbm ti/te no max
    # n_old = 27094  # #117 dbm ti/te in alphaMHD
    # n_old = 27093  # #117 dbm gxi calculated
    # n_old = 27091  # #117 w20 kparaN in disp change
    # n_old = 27073  # #114 w20 only
    # n_old = 27064  # #111 w20 min gma diffusivity zepsqrt
    # n_old = 27062  # #111 w20 velocity min zepslon
    # n_old = 27061  # #111 min eps zepslon
    # n_old = 27060  # #111 min grad zepslon + limit small grad
    # n_old = 27059  # #111 min grad zepslon
    # n_old = 27056  # #111 max ti/te = 10
    # n_old = 27055  # #111 FIH min zepslon
    # n_old = 27054  # #111 W20 Default
    # n_old = 27051  # #111 DBM Default ti/te = 1
    # n_old = 27046  # #109 MTM 2000
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


    ## NSTU
    # scan_data['121123K55'] = [n_new, n_old]
    # scan_data['202946A02'] = [n_new, n_old]
    # scan_data['203531A08'] = [n_new, n_old]
    # scan_data['203592A02'] = [n_new, n_old]
    # scan_data['204179A01'] = [n_new, n_old]
    # scan_data['204198A03'] = [n_new, n_old]
    # scan_data['204201A01'] = [n_new, n_old]
    # scan_data['204509A02'] = [n_new, n_old]
    # scan_data['204511A03'] = [n_new, n_old]
    # scan_data['204519A03'] = [n_new, n_old]
    # scan_data['204556A01'] = [n_new, n_old]
    # scan_data['204651A04'] = [n_new, n_old]
    # scan_data['204665A01'] = [n_new, n_old]
    # scan_data['204963A08'] = [n_new, n_old]
    # scan_data['204980A05'] = [n_new, n_old]
    # scan_data['205042A02'] = [n_new, n_old]

    ##

    ## NSTX
    scan_data['120968A02'] = [n_new, n_old]
    scan_data['120982A09'] = [n_new, n_old]
    scan_data['129016A04'] = [n_new, n_old]
    scan_data['129017A04'] = [n_new, n_old]
    scan_data['129018A02'] = [n_new, n_old]
    scan_data['129019A02'] = [n_new, n_old]
    scan_data['129020A02'] = [n_new, n_old]
    scan_data['129041A10'] = [n_new, n_old]
    scan_data['133964Z02'] = [n_new, n_old]
    scan_data['134020N01'] = [n_new, n_old]
    scan_data['138536A01'] = [n_new, n_old]
    scan_data['141007A10'] = [n_new, n_old]
    scan_data['141031A01'] = [n_new, n_old]
    scan_data['141032A01'] = [n_new, n_old]
    scan_data['141040A01'] = [n_new, n_old]
    scan_data['141716A80'] = [n_new, n_old]
    scan_data['204100J02'] = [n_new, n_old]
    scan_data['204202Z02'] = [n_new, n_old]

    ## D3D
    scan_data['101381A01'] = [n_new, n_old]
    scan_data['101391A07'] = [n_new, n_old]
    scan_data['141069A01'] = [n_new, n_old]
    scan_data['142111A03'] = [n_new, n_old]
    scan_data['144226A01'] = [n_new, n_old]
    scan_data['147634A02'] = [n_new, n_old]
    scan_data['150840T02'] = [n_new, n_old]
    scan_data['175275K01'] = [n_new, n_old]
    scan_data['175288A01'] = [n_new, n_old]
    scan_data['176523L01'] = [n_new, n_old]
    scan_data['179415P02'] = [n_new, n_old]
    scan_data['183743H01'] = [n_new, n_old]
    scan_data['184822M01'] = [n_new, n_old]

    # ## EAST
    scan_data['85122L01'] = [n_new, n_old]
    scan_data['85606W02'] = [n_new, n_old]
    scan_data['85610W01'] = [n_new, n_old]
    scan_data['90328W02'] = [n_new, n_old]
    scan_data['90949R01'] = [n_new, n_old]
    scan_data['100131N01'] = [n_new, n_old]
    scan_data['100137N01'] = [n_new, n_old]
    scan_data['100205N01'] = [n_new, n_old]
    scan_data['100206N01'] = [n_new, n_old]
    scan_data['101085N02'] = [n_new, n_old]
    scan_data['102054N16'] = [n_new, n_old]
    scan_data['113944B01'] = [n_new, n_old]
    scan_data['128474X01'] = [n_new, n_old]

    # ## KSTR
    scan_data['16295A00'] = [n_new, n_old]
    scan_data['16325A00'] = [n_new, n_old]
    scan_data['17231P02'] = [n_new, n_old]
    scan_data['18399P01'] = [n_new, n_old]
    scan_data['18402H01'] = [n_new, n_old]
    scan_data['18476H01'] = [n_new, n_old]
    scan_data['18477P02'] = [n_new, n_old]
    scan_data['18492D01'] = [n_new, n_old]
    scan_data['18499D01'] = [n_new, n_old]
    scan_data['18602P03'] = [n_new, n_old]
    scan_data['22663J01'] = [n_new, n_old]
    scan_data['22773J12'] = [n_new, n_old]
    scan_data['22937T01'] = [n_new, n_old]
    scan_data['25768R01'] = [n_new, n_old]
    scan_data['26607A01'] = [n_new, n_old]

    ## MAST
    scan_data['08505Z06'] = [n_new, n_old]
    scan_data['22341P01'] = [n_new, n_old]
    scan_data['29976P01'] = [n_new, n_old]
    scan_data['45083P01'] = [n_new, n_old]
    scan_data['45163P01'] = [n_new, n_old]
    scan_data['45238P01'] = [n_new, n_old]
    scan_data['47014P01'] = [n_new, n_old]

    ## ITER
    scan_data['20102A12']  = [n_new, n_old]
    scan_data['38265A28']  = [n_new, n_old]
    scan_data['38275A19']  = [n_new, n_old]
    scan_data['38285A42']  = [n_new, n_old]
    scan_data['38530A80']  = [n_new, n_old]
    scan_data['50000A10']  = [n_new, n_old]
    scan_data['59100A05']  = [n_new, n_old]
    scan_data['80200A13']  = [n_new, n_old]

    ## JET
    scan_data['84599T01'] = [n_new, n_old]
    scan_data['86911T01'] = [n_new, n_old]
    scan_data['87215T01'] = [n_new, n_old]
    scan_data['87261T01'] = [n_new, n_old]

    # Break scan data up into smaller dictionaries
    cdfList = list(scan_data.keys())
    keyCount = len(cdfList)
    idx1 = int(0.25 * keyCount)
    idx2 = int(0.50 * keyCount)
    idx3 = int(0.75 * keyCount)

    # cdfList = cdfList[:idx2]
    # cdfList = cdfList[idx2:]

    idx1 = int(0.33333 * keyCount) + 1
    idx2 = int(0.66666 * keyCount) + 1

    # cdfList = cdfList[:idx1]
    # cdfList = cdfList[idx1:idx2]
    cdfList = cdfList[idx2:]

    scan_data_new = {}

    for key in cdfList:
        scan_data_new[key] = scan_data[key]

    scan_data = scan_data_new

    plot_options = contourdata.PlotOptions(
        # xmin=0.2,
        # xmax=0.8,
        # ymin=0.2,
        # ymax=0.6,
        showfig=0,
        savefig=1,
        savedata=0,
        plotidentical=1,
        saveappend='',
        difftype='diff',  # 'diff', 'rel', 'absdiff', 'absrel', 'ratio'
    )

    main(vars_to_plot, scan_data, plot_options)
