#!/usr/bin/python3

"""Saved contour plots for the ETGM paper"""

# Standard Packages
import sys; sys.path.insert(0, '../'); sys.path.insert(0, '../../')
import logging

# 3rd Party Packages
import matplotlib.pyplot as plt

# Local Packages
import modules.options  # Import here even though it's unused to avoid issues
from plotting.modules.plotstyles import PlotStyles, StyleType
from plotting.plot_contour import main


_log = logging.getLogger(__name__)


# Run this file directly to plot scanned variable profiles from previously created scanned data
if __name__ == '__main__':
    scan_data = {}

    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.SINGLE1B,
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
    })

    # Text to append to the end of the generated save name
    savenameend = ''

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
    # vars_to_plot = ['var_to_scan', 'gmaETGM', 'lareunit', 'alphamhdunit', 'xteETGM', 'xte2ETGM', 'gaveETGM']
    # vars_to_plot = ['var_to_scan']

    '''
    Scan Data:
    * Uncomment the lines you wish to use
        - keys (str): The runid of the scan
        - values (list of int): The scan_numbers to plot from
    '''
    """exbs off"""
    vars_to_plot = [
        'gmaETGM', 'omgETGM', 'xteETGM', 'xte2ETGM',
        'xteETG', 'walfvenunit', 'phi2ETGM', 'Apara2ETGM', 'satETGM',
        'gaveETGM', 'kyrhosETGM', 'kyrhoeETGM', 'kpara2ETGM', 'fleETGM', 'omegateETGM',
        'omegadETGM', 'omegasETGM', 'omegasetaETGM', 'omegadiffETGM', 'gammadiffETGM',
        'gne', 'gte', 'shat_gxi', 'etae', 'betaeunit', 'wexbs', 'bunit', 'te', 'ne', 'q'
    ]

    """SUMMED MODES"""
    # scan_data['121123K55'] = [4]
    # scan_data['120968A02'] = [4]
    # scan_data['120982A09'] = [4]
    # scan_data['129016A04'] = [4]
    # scan_data['129017A04'] = [4]
    # scan_data['129018A02'] = [4]
    # scan_data['129019A02'] = [4]
    # scan_data['129020A02'] = [4]
    # scan_data['129041A10'] = [4]
    # scan_data['138536A01'] = [1475]  # 1236 (7), 1256 (8), 1269 (9)
    # scan_data['141007A10'] = [4]
    # scan_data['141031A01'] = [4]
    # scan_data['141032A01'] = [4]
    # scan_data['141040A01'] = [4]
    # scan_data['141716A80'] = [4]
    # scan_data['132017T01'] = [4]
    # scan_data['141552A01'] = [4]

    """MAX MODE"""
    # scan_data['121123K55'] = [5]
    # scan_data['120968A02'] = [5]
    # scan_data['120982A09'] = [5]
    # scan_data['129016A04'] = [5]
    # scan_data['129017A04'] = [5]
    # scan_data['129018A02'] = [5]
    # scan_data['129019A02'] = [5]
    # scan_data['129020A02'] = [5]
    # scan_data['129041A10'] = [5]
    # scan_data['138536A01'] = [1476]  # 1236 (7), 1256 (8), 1269 (9)
    # scan_data['141007A10'] = [5]
    # scan_data['141031A01'] = [5]
    # scan_data['141032A01'] = [5]
    # scan_data['141040A01'] = [5]
    # scan_data['141716A80'] = [5]
    # scan_data['132017T01'] = [5]
    # scan_data['141552A01'] = [5]

    """exbs on"""
    # vars_to_plot = ['gmaETGM', 'xteETGM', 'xte2ETGM']
    # scan_data['121123K55'] = [2]
    # scan_data['120968A02'] = [2]
    # scan_data['120982A09'] = [2]
    # scan_data['129016A04'] = [2]
    # scan_data['129017A04'] = [2]
    # scan_data['129018A02'] = [2]
    # scan_data['129019A02'] = [2]
    # scan_data['129020A02'] = [2]
    # scan_data['129041A10'] = [2]
    # scan_data['138536A01'] = [1237]
    # scan_data['141007A10'] = [2]
    # scan_data['141031A01'] = [2]
    # scan_data['141032A01'] = [2]
    # scan_data['141040A01'] = [2]
    # scan_data['141716A80'] = [2]
    # scan_data['132017T01'] = [2]
    # scan_data['141552A01'] = [2]

    # scan_data['138536A01'] = [1779, 1780]
    # scan_data['138536A01'] = [i for i in range(1716, 1738 + 1)]
    # scan_data['138536A01'] = [i for i in range(1749, 1750 + 1)]
    # scan_data['138536A01'] = [i for i in range(1756, 1763 + 1)]
    # scan_data['138536A01'] = [i for i in range(1779, 1780 + 1)]
    # scan_data['138536A01'] = [
    #     i for i in [*range(1716, 1738 + 1), *range(1749, 1750 + 1), *range(1756, 1763 + 1), *range(1779, 1780 + 1), *range(1860, 1861 + 1)]
    # ]

    main(vars_to_plot, scan_data, savenameend=savenameend, savefig=1, savedata=1)
