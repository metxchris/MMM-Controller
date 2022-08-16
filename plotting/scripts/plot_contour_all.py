#!/usr/bin/python3

"""Saved contour plots for the ETGM paper"""

# Standard Packages
import sys; sys.path.insert(0, '../'); sys.path.insert(0, '../../')
import logging

# 3rd Party Packages
import matplotlib.pyplot as plt

# Local Packages
import modules.options  # Import here even though it's unused to avoid issues
from modules.variables import OutputVariables
from plotting.modules.plotstyles import PlotStyles, StyleType
from plotting.plot_contour import main, PlotOptions


_log = logging.getLogger(__name__)


# Run this file directly to plot scanned variable profiles from previously created scanned data
if __name__ == '__main__':
    scan_data = {}

    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.AIP,
    )

    plot_options = PlotOptions(
        xmin=0,
        xmax=1,
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
        # 'text.usetex': True,
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
    # vars_to_plot = ['gmaETGM', 'omgETGM', 'xteETGM', 'xte2ETGM', 'kyrhoeETGM', 'kyrhosETGM', 'gave', 'var_to_scan']
    # vars_to_plot = ['gmaETGM', 'omgETGM', 'xteETGM']
    # vars_to_plot = ['var_to_scan', 'gmaETGM', 'lareu', 'alphamhdu', 'xteETGM', 'xte2ETGM', 'gaveETGM']
    # vars_to_plot = ['var_to_scan']

    '''
    Scan Data:
    * Uncomment the lines you wish to use
        - keys (str): The runid of the scan
        - values (list of int): The scan_numbers to plot from
    '''
    vars_to_plot = [
        'gmaETGM', 'omgETGM', 'xteETGM', 'xte2ETGM', 'xteETG',
        'gmaMTM', 'omgMTM', 'xteMTM',
        'gmaW20ii', 'gmaW20ie', 'gmaW20ei', 'gmaW20ee',
        'omgW20ii', 'omgW20ie', 'omgW20ei', 'omgW20ee',
        'xteW20', 'xtiW20', 'xdiW20',
        'xdz', 'xvt', 'xvp',
        'xte', 'xti', 'xdi',
        # 'fte', 'fti', 'fdi', 'fdz',
    ]
    vars_to_plot = [
        'xdz', 'xvt', 'xvp',
        'xte', 'xti', 'xdi',
        'vcz', 'vct', 'vcp',
    ]
    # vars_to_plot = ['gmaW20ii', 'gmaW20ei', 'gmaW20ee', 'gmaW20ie']
    # vars_to_plot = ['omgW20ii', 'omgW20ei', 'omgW20ee', 'omgW20ie']
    # vars_to_plot = ['vcz', 'vct', 'vcp']
    # vars_to_plot = ['xte', 'xti', 'xdi', 'xdz', 'xvt', 'xvp', 'vcz', 'vcp', 'vct']
    vars_to_plot = ['xte', 'xti', 'xdz', 'xvt', 'xvp',]
    vars_to_plot = OutputVariables().get_dbm_vars()

    # n, s = 1, 'OPT2'
    # n, s = 2, 'OPT1'
    n, s = 1000, 'OLD'
    n, s = 1001, 'OLDW'  # WEXB
    n, s = 1002, 'NEW'  # KAPPA
    n, s = 1003, 'NEWW'  # KAPPA, WEXB
    n, s = 1004, 'NEWG'  # GXI
    n, s = 1005, 'NEWGW'  # GXI, WEXB
    n, s = 1006, 'NEWNO'  # KAPPA, NO OPT

    nlist = [6100, 6101]
    slist = ['', '']

    # scan_data['118341T54'] = [350]

    for n, s in zip(nlist, slist):
        # # NSTU
        # scan_data['121123K55'] = [n]
        # NSTX
        scan_data['120968A02'] = [n]
        scan_data['120982A09'] = [n]
        scan_data['129016A04'] = [n]
        scan_data['129017A04'] = [n]
        scan_data['129018A02'] = [n]
        scan_data['129019A02'] = [n]
        scan_data['129020A02'] = [n]
        scan_data['129041A10'] = [n]
        scan_data['138536A01'] = [n]
        scan_data['141007A10'] = [n]
        scan_data['141031A01'] = [n]
        scan_data['141032A01'] = [n]
        scan_data['141040A01'] = [n]
        scan_data['141716A80'] = [n]
        # D3D
        scan_data['98777V06']  = [n]
        scan_data['101381J05'] = [n]
        scan_data['101381T31'] = [n]
        scan_data['101391J08'] = [n]
        scan_data['118341T54'] = [n]
        scan_data['132017T01'] = [n]
        scan_data['132411T02'] = [n]
        scan_data['132498J05'] = [n]
        scan_data['141552A01'] = [n]
        scan_data['150840T02'] = [n]
        scan_data['153283T50'] = [n]
        # # EAST
        # scan_data['85126T02'] = [n]
        # scan_data['85610T01'] = [n]
        # scan_data['85122T04'] = [n]
        # scan_data['80208T04'] = [n]
        # scan_data['90328T01'] = [n]
        # # KSTR
        # scan_data['15334T03'] = [n]
        # scan_data['16297T01'] = [n]
        # scan_data['16299T01'] = [n]
        # scan_data['16325T10'] = [n]
        # scan_data['16901T01'] = [n]
        # scan_data['18399T05'] = [n]
        # scan_data['18400T01'] = [n]
        # scan_data['18402T01'] = [n]
        # scan_data['18404T05'] = [n]
        # scan_data['18476T02'] = [n]
        # scan_data['18477T01'] = [n]
        # scan_data['18492T01'] = [n]
        # scan_data['18495T01'] = [n]
        # scan_data['18499T01'] = [n]
        # scan_data['18602T01'] = [n]
        # # MAST
        # scan_data['08505Z06'] = [n]
        # scan_data['18696B01'] = [n]
        # scan_data['22341P37'] = [n]
        # scan_data['24899R05'] = [n]
        # scan_data['27527M34'] = [n]
        # scan_data['29271A01'] = [n]
        # scan_data['29976U69'] = [n]
        # scan_data['45424H01'] = [n]
        # # ITER
        # scan_data['18696R06'] = [n]
        # scan_data['38265R80'] = [n]
        # scan_data['50000A10'] = [n]
        # scan_data['80200A13'] = [n]
        # scan_data['80300A02'] = [n]

        
        # ???
        # scan_data['16949T02'] = [1]


        main(vars_to_plot, scan_data, plot_options, savenameend=s, savefig=1, savedata=0)
