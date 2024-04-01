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

    # OLD
    # vars_to_plot += ['xti', 'xte', 'xde', 'xdz', 'xvp', 'xvt']
    # vars_to_plot += ['fti', 'fte', 'fde', 'fdz', 'fvp', 'fvt']

    # vars_to_plot = ['gxi','gaveW20i', 'gaveW20e']
    # vars_to_plot += ['fti', 'fte', 'fne', 'fnz', 'fvp', 'fvt']
    # vars_to_plot += ['xti', 'xte', 'xne', 'xnz', 'xvp', 'xvt']
    # vars_to_plot = ['xtiW20', 'xteW20', 'xdeW20', 'xdz', 'xvt', 'xvp']
    vars_to_plot += ['fte', 'xte']
    # vars_to_plot += ['fti', 'fte', 'fne']
    # vars_to_plot += ['xti', 'xte', 'xne']
    # vars_to_plot += ['gti', 'gte', 'gne', 'gnz', 'gvt', 'gvp']
    # vars_to_plot += ['nR8TOMSQZ', 'nWarning', 'nError']
    # vars_to_plot = ['gma0W20', 'gmagW20', 'gmaW20']
    # vars_to_plot = ['gaveETGM', 'gmaW20e', 'gmaW20']
    # vars_to_plot += ['nR8TOMSQZ']
    # vars_to_plot = ['fti', 'fte', 'xvt']
    # vars_to_plot = ['gaveW20i', 'gaveW20e']

    num = 27117  # #122 v9.1.0 mtm
    # num = 27116  # #122 v9.1.0 etgm horton 2
    # num = 27115  # #122 v9.1.0 etgm horton 1
    # num = 27114  # #122 v9.1.0 etgm alternate
    # num = 27113  # #122 v9.1.0 etgm default
    # num = 27112  # #122 v9.1.0 dbm alternate
    # num = 27111  # #122 v9.1.0dbm default
    # num = 27110  # #122 v9.1.0 w20 
    # num = 27101  # #117 dbm shat_e
    # num = 27100  # #117 dbm default vei (no etanc)
    # num = 27099  # #117 dbm default xne, dbm nh_ne, ti/te no max
    # num = 27098  # #117 dbm nh_ne, ti/te no max
    # num = 27097  # #117 dbm ti/te max 3
    # num = 27095  # #117 dbm ti/te = z_ti/te
    # num = 27094  # #117 dbm ti/te in alphaMHD
    # num = 27093  # #117 dbm gxi calculated
    # num = 27092  # #117 w20 gxi calculated
    # num = 27091  # #117 w20 kparaN in disp change
    # num = 27084  # #115 w20 eigensolver fix
    # num = 27083  # #114 w20 no negative chi
    # num = 27082  # #114 mtm
    # num = 27081  # #114 etgm horton 2
    # num = 27080  # #114 etgm horton 1
    # num = 27079  # #114 etgm alt
    # num = 27078  # #114 etgm default
    # num = 27077  # #114 dbm alternate
    # num = 27076  # #114 dbm default
    # num = 27075  # #114 dbm default no wexb
    # num = 27074  # #114 w20 only no wexb
    # num = 27073  # #114 w20 only
    # num = 27072  # #114 Channel swap, xne, xnz rename
    # num = 27065  # #111 Newton-Rahpson
    # num = 27064  # #111 min gma diffusivity zepsqrt
    # num = 27063  # #111 kap1 no max
    # num = 27062  # #111 w20 velocity min zepslon
    # num = 27061  # #111 min eps zepslon
    # num = 27060  # #111 min grad zepslon + limit small grad
    # num = 27059  # #111 min grad zepslon
    # num = 27058  # #111 W20 no max ti/te
    # num = 27057  # #111 W20 XI(4) term
    # num = 27056  # #111 W20 max ti/te = 10
    # num = 27055  # #111 W20 FIH min zepslon
    # num = 27054  # #111 W20 default
    # num = 27053  # #111 DBM Default no sum n + T
    # num = 27052  # #111 DBM Default ti/te normal, max 5
    # num = 27051  # #111 DBM Default ti/te = 1
    # num = 27050  # #109 MTM 400
    # num = 27049  # #109 MTM 100
    # num = 27047  # #109 MTM 200
    # num = 27046  # #109 MTM 2000
    # num = 27045  # #108 Limiting diff due to small grad, increase p to 1E1
    # num = 27041  # #108 Limiting diff due to small grad
    # num = 27040  # #108 Default
    # num = 27032  # Major w20 fixes, + R_curv_min = 0.1
    # num = 27031  # Major w20 fixes, + R_curv_min = 0.01
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

    # ## NSTU
    # scan_data['121123K55'] = [num]
    # scan_data['202946A02'] = [num]
    # scan_data['203531A08'] = [num]
    # scan_data['203592A02'] = [num]
    # scan_data['204179A01'] = [num]
    # scan_data['204198A03'] = [num]
    # scan_data['204201A01'] = [num]
    # scan_data['204509A02'] = [num]
    # scan_data['204511A03'] = [num]
    # scan_data['204519A03'] = [num]
    # scan_data['204556A01'] = [num]
    # scan_data['204651A04'] = [num]
    # scan_data['204665A01'] = [num]
    # scan_data['204963A08'] = [num]
    # scan_data['204980A05'] = [num]
    # scan_data['205042A02'] = [num]



    ## NSTX
    scan_data['120968A02'] = [num]
    scan_data['120982A09'] = [num]
    scan_data['129016A04'] = [num]
    scan_data['129017A04'] = [num]
    scan_data['129018A02'] = [num]
    scan_data['129019A02'] = [num]
    scan_data['129020A02'] = [num]
    scan_data['129041A10'] = [num]
    scan_data['133964Z02'] = [num]
    scan_data['134020N01'] = [num]
    scan_data['138536A01'] = [num]
    scan_data['141007A10'] = [num]
    scan_data['141031A01'] = [num]
    scan_data['141032A01'] = [num]
    scan_data['141040A01'] = [num]
    scan_data['141716A80'] = [num]
    scan_data['204100J02'] = [num]
    scan_data['204202Z02'] = [num]

    ## D3D
    scan_data['101381A01'] = [num]
    scan_data['101391A07'] = [num]
    scan_data['141069A01'] = [num]
    scan_data['142111A03'] = [num]
    scan_data['144226A01'] = [num]
    scan_data['147634A02'] = [num]
    scan_data['150840T02'] = [num]
    scan_data['175275K01'] = [num]
    scan_data['175288A01'] = [num]
    scan_data['176523L01'] = [num]
    scan_data['179415P02'] = [num]
    scan_data['183743H01'] = [num]
    scan_data['184822M01'] = [num]

    ## EAST
    scan_data['85122L01']  = [num]
    scan_data['85606W02']  = [num]
    scan_data['85610W01']  = [num]
    scan_data['90328W02']  = [num]
    scan_data['90949R01']  = [num]
    scan_data['100131N01'] = [num]
    scan_data['100137N01'] = [num]
    scan_data['100205N01'] = [num]
    scan_data['100206N01'] = [num]
    scan_data['101085N02'] = [num]
    scan_data['102054N16'] = [num]
    scan_data['113944B01'] = [num]
    scan_data['128474X01'] = [num]

    ## KSTR
    scan_data['16295A00'] = [num]
    scan_data['16325A00'] = [num]
    scan_data['17231P02'] = [num]
    scan_data['18399P01'] = [num]
    scan_data['18402H01'] = [num]
    scan_data['18476H01'] = [num]
    scan_data['18477P02'] = [num]
    scan_data['18492D01'] = [num]
    scan_data['18499D01'] = [num]
    scan_data['18602P03'] = [num]
    scan_data['22663J01'] = [num]
    scan_data['22773J12'] = [num]
    scan_data['22937T01'] = [num]
    scan_data['25768R01'] = [num]
    scan_data['26607A01'] = [num]

    ## MAST
    scan_data['08505Z06'] = [num]
    scan_data['22341P01'] = [num]
    scan_data['29976P01'] = [num]
    scan_data['45083P01'] = [num]
    scan_data['45163P01'] = [num]
    scan_data['45238P01'] = [num]
    scan_data['47014P01'] = [num]

    ## ITER
    scan_data['20102A12']  = [num]
    scan_data['38265A28']  = [num]
    scan_data['38275A19']  = [num]
    scan_data['38285A42']  = [num]
    scan_data['38530A80']  = [num]
    scan_data['50000A10']  = [num]
    scan_data['59100A05']  = [num]
    scan_data['80200A13']  = [num]

    ## JET
    scan_data['84599T01'] = [num]
    scan_data['86911T01'] = [num]
    scan_data['87215T01'] = [num]
    scan_data['87261T01'] = [num]

    cdfList = list(scan_data.keys())
    keyCount = len(cdfList)
    idx1 = int(0.25 * keyCount) + 1
    idx2 = int(0.50 * keyCount) + 1
    idx3 = int(0.75 * keyCount) + 1

    # cdfList = cdfList[:idx1]
    # cdfList = cdfList[idx1:idx2]
    # cdfList = cdfList[idx2:idx3]
    cdfList = cdfList[idx3:]

    scan_data_new = {}

    for key in cdfList:
        scan_data_new[key] = scan_data[key]

    scan_data = scan_data_new

    plot_options = contourdata.PlotOptions(
        # xmin=0.60,
        # xmax=0.75,
        # ymin=0.5,
        # ymax=500,
        showfig=0,
        savefig=1,
        savedata=0,
        smoothing=0,
        saveappend='',
    )

    ## Run program
    main(vars_to_plot, scan_data, plot_options)
