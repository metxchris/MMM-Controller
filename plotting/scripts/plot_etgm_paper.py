#!/usr/bin/python3

"""General non-automated plotting for the ETGM paper"""


# Standard Packages
import sys; sys.path.insert(0, '../'), sys.path.insert(0, '../../')
import logging
import io

# 3rd Party Packages
import matplotlib.pyplot as plt
import numpy as np
from PyQt5.QtGui import QImage
from PyQt5.QtWidgets import QApplication

# Local Packages
import modules.options
import modules.constants as constants
import modules.datahelper as datahelper
import modules.utils as utils
from plotting.modules.plotstyles import PlotStyles, StyleType
from plotting.plot_variables import FigData, PlotDataCsv, PlotDataCdf, main


_log = logging.getLogger(__name__)


if __name__ == '__main__':
    """Run this module directly to plot variable data stored in CDF or CSV files"""

    utils.init_logging()

    # Initialize visual styles for the plot
    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.SINGLE1B,
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
        # 'text.usetex': True,
    })

    # Define settings for the plot
    fig_data = FigData(
        replace_offset_text=0,
        allow_title_runid=1,
        allow_title_time=1,
        allow_title_factor=True,
        allow_title_rho=True,
        invert_y_axis=False,
        invert_x_axis=False,
        nomralize_y_axis=False,
        nomralize_x_axis=False,
        # title_override=r'using $B_\mathrm{unit} = 2\hat{\rho}/r \cdot \mathtt{TRFMP}_\mathrm{B} \cdot \mathtt{GXI}$',
        # title_override=r'using $B_\mathrm{unit} = \hat{\rho}/(\pi r) \cdot \mathtt{TRFMP}_\mathrm{B} \cdot \mathtt{GXI}$',
        # title_override=r'using $\rho = \sqrt{\chi_T/(\pi B_T)}$',
        # title_override=r'using $\rho = \sqrt{\chi_T/(\pi B_T)}$',
        # title_override=r'using $\rho = \sqrt{\chi_T/(\pi B_0)}$',
        # title_override=r'using $\rho = \sqrt{2\chi_T/B_T}$',
        # title_override=r'using $\rho = \sqrt{2\chi_T/B_0}$',
        # title_override=r'using $\rho \sim \sqrt{\chi_T/B_0}$',
        # title_override=r'using $B_\mathrm{unit} = \hat{\rho}/(\pi r) \cdot \mathtt{TRFMP}_\mathrm{B} \cdot \mathtt{GXI}$',
        # title_override=r'Changed $\rho$',
        # title_override=r'Changed $\rho$, $B_0$',
        # title_override=r'Changed $\rho$, $B_0$, $\overline{G}$, $B_\mathrm{unit}$',
        # title_override=r'$\omega_{E \times B}$ off',
        # title_override=r'$k_y\rho_\mathrm{s}$ Scan Counts',
        # title_override=r'$k_y\rho_\mathrm{e}$ Scan Counts',
        # title_override=r'Max vs Summed Modes',
        # title_override=r'Horton ETG',
        # title_override=r'Saturation Exponent = 2',
        # title_override=r'Electrostatic Case',
        # title_override=r'Max Mode',
        # title_override=r'Summed Modes',
        # title_override=r'Summed Modes ($m$=2)',
        # title_override=r'Max Mode ($m$=2)',
        # title_override=r'Single $k_y\rho_\mathrm{s}$ vs Sum ($m$=1)',
        # title_override=r'Convergence Speed ($m=1$)',
        # title_override=r'zami(4,4) = 0',
        title_override=r' ',
        # ylabel_override=r'$|g_\mathrm{ne}|$',
        ylabel_override=r'',
        xlabel_override=r'',
        summed_modes=0,
    )

    # Define data for the plot
    r = '138536A01'
    # r = '141031A01'
    # r = '129016A04'
    n1 = 95
    n2 = n1 + 1
    fig_data.set(

        # PlotDataCdf('18476T02', yname='ti', zval=0.4, timeplot=1, legend='Prediction'),
        # PlotDataCdf('18476P01', yname='ti', zval=0.4, timeplot=1, legend='Experiment'),
        # xmin=1, xmax=7,
        # PlotDataCdf('18476T02', yname='xki', zval=7, legend='Total'),
        # PlotDataCdf('18476T02', yname='xkimmm07', zval=7, legend='MMM'),
        # xmin=0.01, xmax=0.8, ymax=6,

        PlotDataCsv('138536A01', 291, yname='xteETGM'),
        PlotDataCsv('138536A01', 292, yname='xteETGM'),
        PlotDataCsv('138536A01', 293, yname='xteETGM'),
        PlotDataCsv('138536A01', 294, yname='xteETGM'),
        # PlotDataCdf('16297T01', yname='beta', zval=0.5, timeplot=1),

        # PlotDataCdf('85610T01', yname='curdlh', zval=4, legend='Low Density Discharge'),
        # PlotDataCdf('85126T02', yname='curdlh', zval=4, legend='High Density Discharge'),
        # ymax=1.28, ymin=0,
        # PlotDataCdf('85610T01', yname='curdoh', zval=4, legend='Low Density Discharge'),
        # PlotDataCdf('85126T02', yname='curdoh', zval=4, legend='High Density Discharge'),
        # ymax=2.5, ymin=0,

        # PlotDataCdf('85610T01', yname='curoh', zval=4, legend='Low Density Discharge'),
        # PlotDataCdf('85126T02', yname='curoh', zval=4, legend='High Density Discharge'),
        # title_override=r'using $\pi r^2$'

        # PlotDataCdf('85610T01', yname='curlh', zval=4, timeplot=1, legend='Low Density Discharge'),
        # PlotDataCdf('85126T02', yname='curlh', zval=4, timeplot=1, legend='High Density Discharge'),
        # xmax=6, ymin=0,

        # PlotDataCdf('85610T01', yname='darea', zval=4, legend='Low Density Discharge'),
        # PlotDataCdf('85126T02', yname='darea', zval=4, legend='High Density Discharge'),

        # PlotDataCdf('85126T02', yname='curdoh', zval=4),
        # PlotDataCdf('85610T02', yname='curdoh', zval=4),

        # PlotDataCdf('101381T31', yname='gte', zval=0.490),
        # PlotDataCdf('120968A02', yname='gte', zval=0.490),

        # PlotDataCdf('18492P02', yname='xke', zval=2.175),
        # PlotDataCdf('18492P02', yname='xki', zval=2.175),
        # xmax=0.91,

        # PlotDataCdf('18399P01', yname='te', zval=10),
        # PlotDataCdf('18399T05', yname='te', zval=10),
        # xmax=0.91, ymax=3.5,

        # PlotDataCdf('85126T02', yname='xke', zval=2.175, legend=r'$\mathtt{CONDE:PR,WNC+XKEPALEO}$'),
        # PlotDataCdf('85126T02', yname='xkemmm07', zval=2.175, legend=r'$\mathtt{XKEMMM07}$'),
        # PlotDataCdf('85126T02', yname='conde', zval=2.175, legend=r'$\mathtt{CONDE}$'),
        # PlotDataCdf('85126T02', yname='condepr', zval=2.175, legend=r'$\mathtt{CONDEPR}$'),
        # xmax=0.91, xmin=0, ymax=4, ymin=0,

        # PlotDataCdf('85610T01', yname='xki', zval=2.175, legend=r'$\mathtt{CONDI:PR,WNC}$'),
        # PlotDataCdf('85610T01', yname='xkimmm07', zval=2.175, legend=r'$\mathtt{XKIMMM07}$'),
        # PlotDataCdf('85610T01', yname='condi', zval=2.175, legend=r'$\mathtt{CONDI}$'),
        # PlotDataCdf('85610T01', yname='condipr', zval=2.175, legend=r'$\mathtt{CONDIPR}$'),
        # PlotDataCdf('85610T01', yname='condiwnc', zval=2.175, legend=r'$\mathtt{CONDIWNC}$'),
        # xmax=0.91, ymax=4, ymin=0,

        # PlotDataCsv(r, 1928,  'xte2ETGM', xname='rho', legend=r'min $g_\mathrm{Te} = 0.5$'),
        # PlotDataCsv(r, 1929, 'xte2ETGM', xname='rho', legend=r'min $g_\mathrm{Te} = 1.0$'),
        # PlotDataCsv(r, 1930, 'xte2ETGM', xname='rho', legend=r'min $g_\mathrm{Te} = 2.0$'),

        # PlotDataCsv(r, 1931, 'xte2ETGM', xname='rho', legend=r'min $g_\mathrm{Te} = 0.5$'),
        # PlotDataCsv(r, 1932, 'xte2ETGM', xname='rho', legend=r'min $g_\mathrm{Te} = 1.0$'),
        # PlotDataCsv(r, 1933, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 1933, 'xte2ETGM', xname='rho'),

        # PlotDataCsv(r, 1840, 'kyrhosETGM', xname='rho'),
        # PlotDataCsv(r, 1841, 'kyrhosETGM', xname='rho'),
        # PlotDataCsv(r, 1807, 'gmaETGM', xname='rho'),

        # PlotDataCsv(r, 1924, 'te', xname='rho', rho_value=0.1),

        # PlotDataCsv(r, 1867, 'gtecritETG', xname='rho', legend='with $Z_\mathrm{eff}$'),
        # PlotDataCsv(r, 1869, 'gtecritETG', xname='rho', legend='without $Z_\mathrm{eff}$'),
        # PlotDataCsv(r, 1869, 'gte', xname='rho'),
        # PlotDataCsv(r, 1869, 'gtecritETG', xname='rho'),

        # PlotDataCsv('TEST', 17, 'gmaDBM', xname='rho'),
        # PlotDataCsv('TEST', 17, 'xteDBM', xname='rho'),

        # PlotDataCsv(r, 1664, 'xteETGM', xname='rho'),  # Max, Calibration
        # PlotDataCsv(r, 1664, 'xte2ETGM', xname='rho'),
        # PlotDataCsv(r, 1665, 'xteETGM', xname='rho'),  # Sum, Calibration
        # PlotDataCsv(r, 1665, 'xte2ETGM', xname='rho'),

        # PlotDataCsv(r, 1773, 'gmaETGM', xname='rho', legend=r'Default'),
        # PlotDataCsv(r, 1772, 'gmaETGM', xname='rho', legend=r'$g_\mathrm{Te}=0$'), #gte0
        # PlotDataCsv(r, 1771, 'gmaETGM', xname='rho', legend=r'$g_\mathrm{ne}=0$'), #gne0
        # PlotDataCsv(r, 1774, 'gmaETGM', xname='rho', legend=r'$g_\mathrm{Te}, g_\mathrm{ne}=0$'), #gne0
        # PlotDataCsv(r, 1765, 'gmaETGM', xname='rho', legend=r'$\beta_\mathrm{e}=0$'),

        # PlotDataCdf('120982A09', yname='nuste', xname='rho', zval=0.620),
        # PlotDataCdf('120968A02', yname='nuste', xname='rho', zval=0.560),
        # PlotDataCdf('129041A10', yname='nuste', xname='rho', zval=0.490),

        # PlotDataCdf('138536A01', yname='lb', xname='rho', zval=0.629),
        # PlotDataCdf('138536A01', yname='lbunit', xname='rho', zval=0.629),
        # PlotDataCdf('138536A01', yname='rmaj', xname='rho', zval=0.629),
        # PlotDataCdf('138536A01', yname='gbu', xname='rho', zval=0.629),

        # PlotDataCdf('129016A04', yname='rmin', xname='rho', zval=0.5),
        # PlotDataCdf('129016A04', yname='gbtor', xname='rho', zval=0.5),

        # PlotDataCsv('TEST', 86, 'gmaMTM', xname='rho'),
        # PlotDataCsv('TEST', n1, 'gmaMTM', xname='rho'),
        # PlotDataCsv('TEST', 86, 'xteMTM', xname='rho'),
        # PlotDataCsv('TEST', n1, 'xteMTM', xname='rho'),
        # PlotDataCsv('TEST', 86, 'kyrhosMTM', xname='rho'),
        # PlotDataCsv('TEST', n1, 'kyrhosMTM', xname='rho'),
        # PlotDataCsv('TEST', 86, 'dbsqprf', xname='rho'),
        # PlotDataCsv('TEST', n1, 'dbsqprf', xname='rho'),

        # PlotDataCsv('TEST', 78, 'xteETGM', xname='rho')right,
        # PlotDataCsv('TEST', n1, 'xteETGM', xname='rho'),

        # PlotDataCdf('138536A01', yname='beta', xname='rho', zval=0.629, transp_calcs=True),
        # PlotDataCdf('138536A01', yname='beta', xname='rho', zval=0.629),

        # PlotDataCdf(r, yname='e_r_grp', xname='rho', zval=0.629, transp_calcs=True),
        # PlotDataCdf(r, yname='e_r_grp', xname='rho', zval=0.629),
        # PlotDataCdf(r, yname='e_r_phi', xname='rho', zval=0.629, transp_calcs=True),
        # PlotDataCdf(r, yname='e_r_phi', xname='rho', zval=0.629),
        # PlotDataCdf(r, yname='e_r_tht', xname='rho', zval=0.629, transp_calcs=True),
        # PlotDataCdf(r, yname='e_r_tht', xname='rho', zval=0.629),

        # PlotDataCdf(r, yname='wexb', xname='rho', zval=0.629, transp_calcs=True),
        # PlotDataCdf(r, yname='wexb', xname='rho', zval=0.629, apply_smoothing=0),

        # PlotDataCdf(r, yname='wexbsgrp', xname='rho', zval=0.629, transp_calcs=True),
        # PlotDataCdf(r, yname='wexb', xname='rho', zval=0.629, apply_smoothing=0),
        # PlotDataCdf(r, yname='wexbstht', xname='rho', zval=0.629, transp_calcs=True),
        # PlotDataCdf(r, yname='wexb', xname='rho', zval=0.629, apply_smoothing=0),
        # PlotDataCdf(r, yname='wexbsphi', xname='rho', zval=0.629, transp_calcs=True),
        # PlotDataCdf(r, yname='wexb', xname='rho', zval=0.629, apply_smoothing=0),


        # PlotDataCdf(r, yname='wexbs_unit1', xname='rho', zval=0.629),
        # PlotDataCdf(r, yname='wexbs_unit2', xname='rho', zval=0.629),
        # PlotDataCdf(r, yname='wexbs_unit1_ratio', xname='rho', zval=0.629),
        # PlotDataCdf(r, yname='wexbs_unit2_ratio', xname='rho', zval=0.629),
        # PlotDataCdf(r, yname='bu_btor', xname='rho', zval=0.629),

        # PlotDataCdf(r, yname='vtoravg', xname='rho', zval=0.629, transp_calcs=True),
        # PlotDataCdf(r, yname='vtoravg', xname='rho', zval=0.629),

        # PlotDataCdf(r, yname='vtorx', xname='rho', zval=0.629),
        # PlotDataCdf(r, yname='vtord', xname='rho', zval=0.629),
        # PlotDataCdf(r, yname='vtorh', xname='rho', zval=0.629),

        # PlotDataCdf(r, yname='zeff', xname='rho', zval=0.629, transp_calcs=True),
        # PlotDataCdf(r, yname='zeff', xname='rho', zval=0.629),

        # PlotDataCsv(r, 753, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 754, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 761, 'xte2ETGM', xname='rho', legend=r'$k_y\rho_s$ max = 0.50'),
        # PlotDataCsv(r, 762, 'xte2ETGM', xname='rho', legend=r'$k_y\rho_s$ max = 0.75'),
        # PlotDataCsv(r, 763, 'xte2ETGM', xname='rho', legend=r'$k_y\rho_s$ max = 1.00'),
        # PlotDataCsv(r, 764, 'xte2ETGM', xname='rho', legend=r'$k_y\rho_s$ max = 2.00'),
        # PlotDataCsv(r, 765, 'xte2ETGM', xname='rho', legend=r'$k_y\rho_s$ max = 50'),
        # PlotDataCsv(r, 766, 'xte2ETGM', xname='rho', legend=r'$k_y\rho_s$ max = 75'),
        # PlotDataCsv(r, 774, 'xteETGM', xname='rho', legend=r'$k_y\rho_s$ max = 100'),
        # PlotDataCsv(r, 777, 'xteETGM', xname='rho', legend=r'$k_y\rho_s$ max = 100'),
        # PlotDataCsv(r, 778, 'xteETGM', xname='rho', legend=r'using $\gamma_{max}$'),

        # PlotDataCsv(r, 795, 'xteETGM', xname='rho', legend=r''),
        # PlotDataCsv(r, 815, 'fte', xname='rho', legend=r''),  # GOOD (includes etg)
        # PlotDataCsv(r, 830, 'xteETGM', xname='rho', legend=r''), #OLD with zepsqrt check
        # PlotDataCsv(r, 836, 'xteETGM', xname='rho', legend=r''), #GOOD


        # PlotDataCsv(r, 838, 'kyrhosETGM', xname='rho', legend=r'Previous'),  #GOOD
        # PlotDataCsv(r, 844, 'gmaETGM', xname='rho', legend=r'CL On'),  # Hem Term
        # PlotDataCsv(r, 845, 'kyrhosETGM', xname='rho', legend=r'CL Off'),  # Hem Term

        # PlotDataCsv(r, 838, 'xteETGM', xname='rho', legend=r''),  #GOOD

        # PlotDataCsv(r, 852, 'xte2ETGM', xname='rho', legend=r''),  # kyrhos 0.05, 200 (500)
        # PlotDataCsv(r, 853, 'kyrhosETGM', xname='rho', legend=r''),  # kyrhos 0.05, 400 (1000)
        # PlotDataCsv(r, 854, 'kyrhosETGM', xname='rho', legend=r''),  # kyrhos 1, 400 (1000)
        # PlotDataCsv(r, 855, 'kyrhosETGM', xname='rho', legend=r''),  # kyrhos 0.1, 400 (1000)
        # PlotDataCsv(r, 856, 'kyrhosETGM', xname='rho', legend=r''),  # kyrhos 0.5, 400 (1000)

        # PlotDataCsv(r, 866, 'gte', xname='rho', legend=r''),  # kyrhos 0.4, 100 (2000)

        # PlotDataCsv(r, 867, 'xteETGM', xname='rho', legend=r''),  # kyrhos 0.1, 100 (100)
        # PlotDataCsv(r, 873, 'xteETGM', xname='rho', legend=r''),  # kyrhos 1, 100 (1000)
        # PlotDataCsv(r, 874, 'xteETGM', xname='rho', legend=r''),  # kyrhos 0.5, 100 (100)
        # PlotDataCsv(r, 875, 'kyrhosETGM', xname='rho', legend=r''),  # kyrhos 0.5, 100 (1000)
        # PlotDataCsv(r, 876, 'kyrhosETGM', xname='rho', legend=r''),  # kyrhos 0.5, 100 (10000)

        # PlotDataCsv(r, 884, 'xteETGM', xname='rho', legend=r''),  # kyrhos 0.05, 100 (10000)
        # PlotDataCsv(r, 885, 'xteETGM', xname='rho', legend=r''),  # kyrhos 0.05, 100 (100) linear
        # PlotDataCsv(r, 888, 'xteETGM', xname='rho', legend=r''),  # kyrhos 0.05, 100 (100) nonlinear
        # PlotDataCsv(r, 889, 'xteETGM', xname='rho', legend=r''),  # kyrhos 0.05, 100 (100) nonlinear2

        # PlotDataCsv(r, 890, 'xteETGM', xname='rho', legend=r''),  # kyrhos 0.05, 100 (100) nonlinear2
        # PlotDataCsv(r, 906, 'xteETGM', xname='rho', legend=r''),  # kyrhos 0.05, 100 (100) nonlinear
        # PlotDataCsv(r, 909, 'fte', xname='rho', legend=r'$k_x/k_y = 0.10$'),  # kyrhos 0.05, 100 (100) nonlinear
        # PlotDataCsv(r, 910, 'fte', xname='rho', legend=r'$k_x/k_y = 0.50$'),  # kyrhos 0.05, 100 (100) nonlinear
        # PlotDataCsv(r, 911, 'fte', xname='rho', legend=r'$k_x/k_y = 0.75$'),  # kyrhos 0.05, 100 (100) nonlinear

        # PlotDataCsv(r, 909, 'gmaETGM', xname='rho', legend=r'min $k_y\rho_s = 0.05$'),  # kyrhos 0.05, 100 (100) nonlinear
        # PlotDataCsv(r, 912, 'gmaETGM', xname='rho', legend=r'min $k_y\rho_s = 1$'),  # kyrhos 0.05, 100 (100) nonlinear

        # PlotDataCsv(r, 912, 'gmaETGM', xname='rho', legend=r'all terms'),  # kyrhos 0.05, 100 (100) nonlinear
        # PlotDataCsv(r, 914, 'gmaETGM', xname='rho', legend=r'zami(3,4)=0'),  # kyrhos 0.05, 100 (100) nonlinear

        # PlotDataCsv(r, 912, 'gmaETGM', xname='rho', legend=r'all terms'),
        # PlotDataCsv(r, 920, 'gmaETGM', xname='rho', legend=r'min $k_y\rho_s = 0.05$'),  # kyrhos 0.05, 100 (100) nonlinear
        # PlotDataCsv(r, 919, 'gmaETGM', xname='rho', legend=r'min $k_y\rho_s = 1$'),  # kyrhos 0.05, 100 (100) nonlinear

        # PlotDataCsv(r, 921, 'xteETGM', xname='rho', legend=r'Eq. (1)'),
        # PlotDataCsv(r, 922, 'xteETGM', xname='rho', legend=r'Eq. (2)'),
        # PlotDataCsv(r, 923, 'xteETGM', xname='rho', legend=r'Eq. (3)'),
        # PlotDataCsv(r, 924, 'xteETGM', xname='rho', legend=r'Eq. (4)'),

        # PlotDataCsv(r, 920, 'xteETGM', xname='rho', legend=r'Sum'),
        # PlotDataCsv(r, 925, 'xteETGM', xname='rho', legend=r'No Sum'),

        # PlotDataCsv(r, 926, 'kyrhosETGM', xname='rho', legend=r'No Sum'),
        # PlotDataCsv(r, 927, 'kyrhosETGM', xname='rho', legend=r'No Sum'),
        # PlotDataCsv(r, 928, 'xteETGM', xname='rho', legend=r'No Sum'),

        # PlotDataCsv(r, 929, 'xteETGM', xname='rho', legend=r'No Sum'),
        # PlotDataCsv(r, 935, 'xteETGM', xname='rho', legend=r'No Sum'),
        # PlotDataCsv(r, 939, 'xteETGM', xname='rho', legend=r'No Sum'),
        # PlotDataCsv(r, 945, 'xteETGM', xname='rho', legend=r'original'),
        # PlotDataCsv(r, 954, 'xteETGM', xname='rho', legend=r'$\gamma$'), # kind of matches xte2


        # Saturation Scan with kx/ky = 0.1, kyrhos=1,80 with 1000 scans
        # PlotDataCsv(r, 1040, 'kyrhosETGM', xname='rho', legend=r'$k_y \rightarrow$saturation, 2'), # norm by |phi|**2
        # PlotDataCsv(r, 1045, 'kyrhosETGM', xname='rho', legend=r'$k_y \rightarrow$saturation, 1'), # norm by |phi|
        # PlotDataCsv(r, 1043, 'kyrhosETGM', xname='rho', legend=r'$k_y \rightarrow$$\gamma_\mathrm{sat}$, 2'), # norm by |phi|**2
        # PlotDataCsv(r, 1044, 'kyrhosETGM', xname='rho', legend=r'$k_y \rightarrow$$\gamma_\mathrm{sat}$, 1'), # norm by |phi|

        # PlotDataCsv(r, 1040, 'omegadETGM', xname='rho', legend=r'$k_y \rightarrow$saturation'), # norm by |phi|**2
        # PlotDataCsv(r, 1043, 'omegadETGM', xname='rho', legend=r'$k_y \rightarrow$$\gamma_\mathrm{sat}$, 2.0'), # norm by |phi|**2
        # PlotDataCsv(r, 1048, 'omegadETGM', xname='rho', legend=r'$k_y \rightarrow$$\gamma_\mathrm{sat}$, 1.5'), # norm by |phi|**1.5
        # PlotDataCsv(r, 1044, 'omegadETGM', xname='rho', legend=r'$k_y \rightarrow$$\gamma_\mathrm{sat}$, 1.0'), # norm by |phi|

        # PlotDataCsv(r, 1040, 'xte2ETGM', xname='rho', legend=r'$k_y \rightarrow$saturation'), # norm by |phi|**2
        # PlotDataCsv(r, 1045, 'xte2ETGM', xname='rho', legend=r'$k_y \rightarrow$saturation, 1'), # norm by |phi|
        # PlotDataCsv(r, 1043, 'kyrhosETGM', xname='rho', legend=r'$k_y \rightarrow$$\gamma_\mathrm{s}$ max$(\gamma)$ exp'), # norm by |phi|**2
        # PlotDataCsv(r, 1064, 'kyrhosETGM', xname='rho', legend=r'$k_y \rightarrow$$\gamma_\mathrm{s}$ max$(\gamma)$ lin'), # norm by |phi|**2
        # PlotDataCsv(r, 1061, 'kyrhosETGM', xname='rho', legend=r'$k_y \rightarrow$$\gamma_\mathrm{s}$ max$(\gamma_s)$ exp'), # # max(gammaS) exp
        # PlotDataCsv(r, 1062, 'kyrhosETGM', xname='rho', legend=r'$k_y \rightarrow$$\gamma_\mathrm{s}$ max$(\gamma_s)$ lin'), # # max(gammaS) linear
        # PlotDataCsv(r, 1059, 'xteETGM', xname='rho', legend=r'$k_y \rightarrow$saturation'),

        # PlotDataCsv(r, 1040, 'satETGM', xname='rho', legend=r'$k_y \rightarrow$saturation'), # norm by |phi|**2
        # PlotDataCsv(r, 1064, 'gmaETGM', xname='rho', legend=r'$k_y \rightarrow$$\gamma_\mathrm{s}$ max$(\gamma)$'), # norm by |phi|**2

        # PlotDataCsv(r, 1064, 'gmaETGM', xname='rho', legend=r'$k_y \rightarrow$$\gamma_\mathrm{s}$ max$(\gamma)$'), # norm by |phi|**2
        # PlotDataCsv(r, 1062, 'gmaETGM', xname='rho', legend=r'$k_y \rightarrow$$\gamma_\mathrm{s}$ max$(\gamma_s)$'), # # max(gammaS) linear

        # PlotDataCsv(r, 1066, 'gmaETGM', xname='rho', legend=r'using $\gamma_\mathrm{max}$'),
        # PlotDataCsv(r, 1068, 'gmaETGM', xname='rho', legend=r'using $\gamma_\mathrm{sat}$'),

        # PlotDataCsv(r, 1068, 'gmaETGM', xname='rho', legend=r'using $\gamma_\mathrm{sat}$, max$(\gamma)$'),
        # PlotDataCsv(r, 1069, 'gmaETGM', xname='rho', legend=r'using $\gamma_\mathrm{sat}$, max$(\gamma_s)$'),

        # PlotDataCsv(r, 1068, 'omgETGM', xname='rho', legend=r'max$(\gamma)$ scan'),
        # PlotDataCsv(r, 1073, 'omgETGM', xname='rho', legend=r'saturation scan'),  # gamma_s

        # PlotDataCsv(r, 1085, 'xteETGM', xname='rho'),  # gamma_s
        # PlotDataCsv(r, 1086, 'xte2ETGM', xname='rho'),  # gamma_s

        # xte comparison
        # PlotDataCsv(r, 1115, 'xte2ETGM', xname='rho', legend='base'),  # gamma_s
        # PlotDataCsv(r, 1117, 'xte2ETGM', xname='rho', legend='base'),  # gamma_s

        # xte comparison with 1/r ~ LB
        # PlotDataCsv(r, 1118, 'xteETGM', xname='rho', legend='base'),  # gamma_s
        # PlotDataCsv(r, 1118, 'xte2ETGM', xname='rho', legend='base'),  # gamma_s

        # PlotDataCsv(r, 1134, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 1134, 'xte2ETGM', xname='rho'),

        # PlotDataCsv(r, 1134, 'satETGM', xname='rho'),
        # PlotDataCsv(r, 1130, 'satETGM', xname='rho'),

        # saturation scan xte2 = xte2 (both xte2)
        # PlotDataCsv(r, 1141, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 1141, 'xte2ETGM', xname='rho'),

        # max growth rate scan xte2 = xte2 (both xte2)
        # PlotDataCsv(r, 1146, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 1146, 'xte2ETGM', xname='rho'),
        # PlotDataCsv(r, 1147, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 1147, 'xte2ETGM', xname='rho'),
        # PlotDataCsv(r, 1152, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 1152, 'xte2ETGM', xname='rho'),

        # PlotDataCsv(r, 1177, 'satETGM', xname='rho'),
        # PlotDataCsv(r, 1189, 'satETGM', xname='rho'),
        # PlotDataCsv(r, 1214, 'satETGM', xname='rho'),

        # PlotDataCsv(r, 1207, 'xte2ETGM', xname='rho'),
        # PlotDataCsv(r, 1211, 'xte2ETGM', xname='rho'),

        # PlotDataCsv(r, 1212, 'xte2ETGM', xname='rho'),
        # PlotDataCsv(r, 1213, 'xte2ETGM', xname='rho'),
        # PlotDataCsv(r, 1214, 'xteETGM', xname='rho'),

        # PlotDataCsv(r, 1237, 'kyrhosETGM', xname='rho', legend=r'min$(k_\mathrm{y}\rho_\mathrm{s}) = 0.05$'),
        # PlotDataCsv(r, 1238, 'kyrhosETGM', xname='rho', legend=r'min$(k_\mathrm{y}\rho_\mathrm{s}) = 0.20$'),
        # PlotDataCsv(r, 1239, 'kyrhosETGM', xname='rho', legend=r'min$(k_\mathrm{y}\rho_\mathrm{s}) = 0.50$'),
        # PlotDataCsv(r, 1240, 'kyrhosETGM', xname='rho', legend=r'min$(k_\mathrm{y}\rho_\mathrm{s}) = 1.00$'),

        # zami(4,4) = 0
        # PlotDataCsv(r, 1249, 'kyrhosETGM', xname='rho', legend=r'min$(k_\mathrm{y}\rho_\mathrm{s}) = 0.05$'),
        # PlotDataCsv(r, 1250, 'kyrhosETGM', xname='rho', legend=r'min$(k_\mathrm{y}\rho_\mathrm{s}) = 0.20$'),
        # PlotDataCsv(r, 1251, 'kyrhosETGM', xname='rho', legend=r'min$(k_\mathrm{y}\rho_\mathrm{s}) = 0.50$'),
        # PlotDataCsv(r, 1252, 'kyrhosETGM', xname='rho', legend=r'min$(k_\mathrm{y}\rho_\mathrm{s}) = 1.00$'),

        # PlotDataCsv(r, 1245, 'xteETGM', xname='rho', legend=r'min$(k_\mathrm{y}\rho_\mathrm{s}) = 0.05$'),
        # PlotDataCsv(r, 1246, 'xteETGM', xname='rho', legend=r'min$(k_\mathrm{y}\rho_\mathrm{s}) = 0.20$'),
        # PlotDataCsv(r, 1247, 'xteETGM', xname='rho', legend=r'min$(k_\mathrm{y}\rho_\mathrm{s}) = 0.50$'),
        # PlotDataCsv(r, 1248, 'xteETGM', xname='rho', legend=r'min$(k_\mathrm{y}\rho_\mathrm{s}) = 1.00$'),

        # PlotDataCsv(r, 1237, 'gmaETGM', xname='rho'),
        # PlotDataCsv(r, 1237, 'omgETGM', xname='rho'),

        # PlotDataCsv(r, 1220, 'xte2ETGM', xname='rho', legend=r'$\epsilon_\mathrm{ne} = 2 L_\mathrm{ne} / R$'),
        # PlotDataCsv(r, 1221, 'xte2ETGM', xname='rho', legend=r'$\epsilon_\mathrm{ne} = 2 L_\mathrm{ne} / L_\mathrm{Bu}$'),
        # PlotDataCsv(r, 1217, 'xte2ETGM', xname='rho'),
        # PlotDataCsv(r, 1218, 'xte2ETGM', xname='rho'),
        # PlotDataCsv(r, 1219, 'xte2ETGM', xname='rho'),

        # 1222+ input Eqs. use FL term
        # PlotDataCsv(r, 1221, 'kyrhosETGM', xname='rho', runname=r'$\mathtt{OLD}$'),
        # PlotDataCsv(r, 1222, 'kyrhosETGM', xname='rho', runname=r'$\mathtt{NEW}$'),
        # PlotDataCsv(r, 1223, 'xteETG', xname='rho', runname=' '),

        # Checking contribution of each individual equation 1-4
        # PlotDataCsv(r, 1227, 'gmaETGM', xname='rho', legend='Eq. 1'),
        # PlotDataCsv(r, 1228, 'gmaETGM', xname='rho', legend='Eq. 2'),
        # PlotDataCsv(r, 1229, 'gmaETGM', xname='rho', legend='Eq. 3'),
        # PlotDataCsv(r, 1230, 'gmaETGM', xname='rho', legend='Eq. 4'),
        # PlotDataCsv(r, 1231, 'xteETGM', xname='rho', legend='Eq. 1'),
        # PlotDataCsv(r, 1232, 'xteETGM', xname='rho', legend='Eq. 2'),
        # PlotDataCsv(r, 1233, 'xteETGM', xname='rho', legend='Eq. 3'),
        # PlotDataCsv(r, 1234, 'xteETGM', xname='rho', legend='Eq. 4'),

        # PlotDataCsv(r, 1165, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 1173, 'xte2ETGM', xname='rho'),
        # PlotDataCsv(r, 1174, 'xte2ETGM', xname='rho'),

        # xteETGM is normal xte2ETGM, and xte2ETGM is generalized formula
        # PlotDataCsv(r, 1201, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 1202, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 1202, 'xte2ETGM', xname='rho'),

        # PlotDataCsv(r, 1105, 'gte', xname='rho', legend='base'),  # gamma_s
        # PlotDataCsv(r, 1105, 'gmaETGM', xname='rho', legend='base'),  # gamma_s
        # PlotDataCsv(r, 1105, 'xteETGM', xname='rho', legend=''),  # gamma_s
        # PlotDataCsv(r, 1105, 'xte2ETGM', xname='rho', legend=''),  # gamma_s
        # PlotDataCsv(r, 1090, 'xteETGM', xname='rho', legend='xte |phi|**2'),  # gamma_s
        # PlotDataCsv(r, 1100, 'gmaETGM', xname='rho', legend='xte |phi|**0'),  # gamma_s
        # PlotDataCsv(r, 1095, 'xte2ETGM', xname='rho', legend='xte |phi|**0'),  # gamma_s
        # PlotDataCsv(r, 1095, 'xte2ETGM', xname='rho', legend='xte |phi|**0'),  # gamma_s

        # PlotDataCsv(r, 1068, 'xteETGM', xname='rho', legend=r'max$(\gamma)$ scan'),
        # PlotDataCsv(r, 1077, 'xteETGM', xname='rho', legend='sat scan'),
        # PlotDataCsv(r, 1080, 'xteETGM', xname='rho', legend='sat scan2'),
        # PlotDataCsv(r, 1084, 'xteETGM', xname='rho', legend='sat scan3'),
        # PlotDataCsv(r, 1079, 'satETGM', xname='rho', legend=r'max$(\gamma_s)$ scan'),
        # PlotDataCsv(r, 1077, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 1077, 'xte2ETGM', xname='rho'),

        # PlotDataCsv(r, 993, 'gmaETGM', xname='rho', legend=r'$k_x/k_y = 0.10$'),
        # PlotDataCsv(r, 998, 'gmaETGM', xname='rho', legend=r'$k_x/k_y = 0.10$'),
        # PlotDataCsv(r, 995, 'satETGM', xname='rho', legend=r'$k_x/k_y = 0.10$'),
        # PlotDataCsv(r, 991, 'gmaETGM', xname='rho', legend=r'$k_x/k_y = 0.50$'),
        # PlotDataCsv(r, 990, 'gmaETGM', xname='rho', legend=r'$k_x/k_y = 0.75$'),
        # PlotDataCsv(r, 992, 'gmaETGM', xname='rho', legend=r'$k_x/k_y = 1.0$'),



        # PlotDataCsv(r, 959, 'xteETGM', xname='rho', legend=r'$\gamma^0$'),
        # PlotDataCsv(r, 949, 'fte', xname='rho', legend=r'$\gamma$'),
        # PlotDataCsv(r, 952, 'xteETGM', xname='rho', legend=r'$\gamma^2$'),
        # PlotDataCsv(r, 953, 'xteETGM', xname='rho', legend=r'$\gamma^2$'),

        # PlotDataCsv(r, 933, 'kyrhosETGM', xname='rho', legend=r'Eq. 1'),
        # PlotDataCsv(r, 931, 'kyrhosETGM', xname='rho', legend=r'Eq. 2'),
        # PlotDataCsv(r, 934, 'kyrhosETGM', xname='rho', legend=r'Eq. 3'),
        # PlotDataCsv(r, 932, 'kyrhosETGM', xname='rho', legend=r'Eq. 4'),

        # PlotDataCsv(r, 879, 'gmaETGM', xname='rho', legend=r''),  # kyrhos 0.5, 100 (100) CL0
        # PlotDataCsv(r, 880, 'gmaETGM', xname='rho', legend=r''),  # kyrhos 0.5, 100 (100) CL0
        # PlotDataCsv(r, 881, 'gmaETGM', xname='rho', legend=r''),  # kyrhos 0.5, 100 (100) CL0

        # PlotDataCsv(r, 791, 'gmaETGM', xname='rho', legend=r'using each $\gamma_j$'),
        # PlotDataCsv(r, 774, 'gmaETGM', xname='rho', legend=r'$k_y\rho_s$ max = 100'),
        # PlotDataCsv(r, 758, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 733, 'gmaETGM', xname='rho', legend='negative sign'),
        # PlotDataCsv(r, 737, 'gmaETGM', xname='rho'),
        # PlotDataCsv(r, 725, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 727, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 694, 'xte2ETGM', xname='rho', legend=r'$k_y\rho_s$ Scan [0.05, 0.5]'),
        # PlotDataCsv(r, 693, 'xte2ETGM', xname='rho', legend=r'$k_y\rho_s$ Scan [0.1, 1]'),
        # PlotDataCsv(r, 692, 'xte2ETGM', xname='rho', legend=r'$k_y\rho_s$ Scan [0.1, 4]'),
        # PlotDataCsv(r, 691, 'xte2ETGM', xname='rho', legend=r'$k_y\rho_s$ Scan [1, 100]'),

        # PlotDataCsv(r, 672, 'xteETGM', xname='rho', legend=r'$0.00\, \omega_\mathrm{E \zvals B}$'),
        # PlotDataCsv(r, 668, 'xteETGM', xname='rho', legend=r'$0.50\, \omega_\mathrm{E \zvals B}$'),
        # PlotDataCsv(r, 670, 'xteETGM', xname='rho', legend=r'$0.25\, \omega_\mathrm{E \zvals B}$'),
        # PlotDataCsv(r, 667, 'xteETGM', xname='rho', legend=r'$0.9\, \omega_\mathrm{E \zvals B}$'),
        # PlotDataCsv(r, 660, 'gmaETGM', xname='rho', legend=r'formula change'),
        # PlotDataCsv(r, 655, 'gmaETGM', xname='rho', legend=r'$\omega_\mathrm{E \zvals B}$ on'),
        # PlotDataCsv(r, 657, 'xteETGM', xname='rho', legend=r'$\omega_\mathrm{E \zvals B}$ (Weiland)'),
        # PlotDataCsv(r, 659, 'xteETGM', xname='rho', legend=r'$\omega_\mathrm{E \zvals B}$ reduced'),

        # PlotDataCsv(r, 653, 'wexb', xname='rho'),

        # PlotDataCsv(r, 625, 'gmaETGM', xname='rho', legend=r'$k_x/k_y = 0.01$'),

        # PlotDataCsv(r, 623, 'xteETGM', xname='rho', legend=r'$\omega_{E\zvals B}$ on'),
        # PlotDataCsv(r, 630, 'xteETGM', xname='rho', legend=r'$\omega_{E\zvals B}$ off'),
        # PlotDataCsv(r, 621, 'gmaETGM', xname='rho', legend=r'$g_{Ti}, g_{ni} = 0$'),

        # PlotDataCsv(r, n1, 'xteETGM', xname='rho', legend='collisional'),
        # PlotDataCsv(r, n2, 'xteETGM', xname='rho', legend='collisionless'),
        # PlotDataCsv(r, n1, 'xte2ETGM', xname='rho', legend='collisional'),
        # PlotDataCsv(r, n2, 'xte2ETGM', xname='rho', legend='collisionless'),
        # PlotDataCsv(r, n1, 'omgETGM', xname='rho', legend='collisional'),
        # PlotDataCsv(r, n2, 'omgETGM', xname='rho', legend='collisionless'),
        # PlotDataCsv(r, n1, 'omgETGM', xname='rho', legend='collisional'),
        # PlotDataCsv(r, n2, 'omgETGM', xname='rho', legend='collisionless'),

        # PlotDataCsv(r, 597, 'gmaETGM', xname='rho', legend='collisional'),
        # PlotDataCsv(r, 598, 'gmaETGM', xname='rho', legend='collisionless'),
        # PlotDataCsv(r, 597, 'omgETGM', xname='rho', legend='collisional'),
        # PlotDataCsv(r, 598, 'omgETGM', xname='rho', legend='collisionless'),
        # PlotDataCsv(r, 597, 'xteETGM', xname='rho', legend='collisional'),
        # PlotDataCsv(r, 598, 'xteETGM', xname='rho', legend='collisionless'),



        # # PlotDataCsv(r, 571, 'shear', xname='rho'),
        # PlotDataCsv(r, 557, 'omgETGM', xname='rho', legend='old collisional'),
        # PlotDataCsv(r, 518, 'gave_gxi', xname='rho'),
        # PlotDataCdf(r, 0.629, 'gave', xname='rho'),
        # PlotDataCsv(r, 495, 'xteETGM'),
        # PlotDataCsv(r, 476, 'shearETGM'),
        # PlotDataCsv(r, 431, 'kyrhoeETGM', runname='default'),
        # PlotDataCsv(r, 431, 'kyrhosETGM', runname='alternate'),Jenko ETG model

        # PlotDataCdf(r, 0.629, 'fle_etgm', xname='rho'),  #, input_points=201, apply_smoothing=True, uniform_rho=True),
        # PlotDataCdf(r, 0.629, 'fle_etgm_shat', xname='rho'),  #, input_points=201, apply_smoothing=True, uniform_rho=True),
        # PlotDataCdf(r, 0.629, 'fle_etgm_shat_gxi', xname='rho'),  #, input_points=201, apply_smoothing=True, uniform_rho=True),
        # PlotDataCdf(r, 0.629, 'fle_etgm_shat_gxi', xname='rho'),

        # PlotDataCdf(r, 0.629, 'gmaxunit_gte', xname='rho'),
        # PlotDataCdf(r, 0.629, 'gmaxunit_gne', xname='rho'),

        # PlotDataCsv(r, 463, 'gave', xname='rho', legend=r'with $\beta^\prime$'),
        # PlotDataCsv(r, 463, 'gave_noss', xname='rho', legend=r'without $\beta^\prime$'),

        # PlotDataCsv(r, 200, 'gmaETGM', xname='rho'),
        # PlotDataCsv(r, 200, 'gmaETGM', xname='rmina'),

        # Old vs New (is the same)
        # PlotDataCsv(r, 97, 'gmaETGM', xname='rho'),
        # PlotDataCsv(r, 213, 'gmaETGM', xname='rmina'),

        # rho vs r/a
        # PlotDataCsv(r, 213, 'gmaETGM', xname='rmina'),
        # PlotDataCsv(r, 213, 'gmaETGM', xname='rho'),

        # Updating Btor definition
        # PlotDataCsv(r, 213, 'gmaETGM', xname='rho', legend=r'Old $B_0$'),
        # PlotDataCsv(r, 214, 'gmaETGM', xname='rho', legend=r'New $B_0$'),

        # Updating Gave definition
        # PlotDataCsv(r, 214, 'gmaETGM', xname='rho', legend=r'$\overline{G}$'),
        # PlotDataCsv(r, 216, 'gmaETGM', xname='rho', legend=r'$\overline{G}^{\star}_{\nabla \rho}$'),

        # Updating Gave definition
        # PlotDataCsv(r, 216, 'gmaETGM', xname='rho', legend=r'$B_\phi$'),
        # PlotDataCsv(r, 217, 'gmaETGM', xname='rho', legend=r'$B_\mathrm{unit}$'),

        # Total change from old to new
        # PlotDataCsv(r, 213, 'gmaETGM', xname='rmina', runname=r'Old'),
        # PlotDataCsv(r, 218, 'gte', xname='rho', runname=r'New'),

        # PlotDataCsv(r, 196, 'bu_btor', xname='rho', legend=r'$B_0 \kappa$'),
        # PlotDataCsv(r, 196, 'bu', xname='rho', legend=r'$B_0 \kappa$'),

        # PlotDataCsv(r, 141, 'xdiETGM', xname='rho' runname='51 Points'),
        # PlotDataCsv(r, 142, 'xdiETGM', xname='rho' runname='201 Points'),
        # PlotDataCsv(r, 143, 'xdiETGM', xname='rho' runname='1001 Points'),
        # PlotDataCsv(r, 144, 'xdiETGM', xname='rho' runname='10001 Points'),

        # PlotDataCsv(r, 218, 'gave', xname='rho', runname='Base'),
        # PlotDataCsv(r, 231, 'gave', xname='rho', runname=r'$g_\mathrm{ne} = 0$'),

        # PlotDataCsv(r, 157, 'shear', xname='rho'),
        # PlotDataCsv(r, 158, 'shat', xname='rho'),
        # PlotDataCsv(r, 159, 'shat_gxi', xname='rho'),

        # PlotDataCsv(r, 169, 'gave', xname='rho', legend=r'$\overline{G}^\star$'),
        # PlotDataCsv(r, 170, 'gave_shat', xname='rho', legend=r'$\overline{G}^{\star}_\kappa$'),
        # PlotDataCsv(r, 171, 'gave_shat_gxi', xname='rho', legend=r'$\overline{G}^{\star}_{\nabla \rho}$'),

        # PlotDataCsv(r, 175, 'gmaETGM', xname='rho', runname=r'$\overline{G}$'),
        # PlotDataCsv(r, 169, 'gmaETGM', xname='rho', runname=r'$\overline{G}^\star$'),

        # PlotDataCsv(r, 176, 'gmaETGM', xname='rho', runname=r'$\overline{G}_\kappa$'),
        # PlotDataCsv(r, 170, 'gmaETGM', xname='rho', runname=r'$\overline{G}^{\star}_\kappa$'),

        # PlotDataCsv(r, 177, 'gmaETGM', xname='rho', runname=r'$\overline{G}_{\nabla \rho}$'),
        # PlotDataCsv(r, 171, 'gmaETGM', xname='rho', runname=r'$\overline{G}^{\star}_{\nabla \rho}$'),

        # PlotDataCsv(r, 175, 'gave', xname='rho', legend=r'$\overline{G}$'),
        # PlotDataCsv(r, 169, 'gave', xname='rho', legend=r'$\overline{G}^\star$'),

        # PlotDataCsv(r, 176, 'gave_shat', xname='rho', legend=r'$\overline{G}_\kappa$'),
        # PlotDataCsv(r, 170, 'gave_shat', xname='rho', legend=r'$\overline{G}^{\star}_\kappa$'),

        # PlotDataCsv(r, 177, 'gave_shat_gxi', xname='rho', legend=r'$\overline{G}_{\nabla \rho}$'),
        # PlotDataCsv(r, 171, 'gave_shat_gxi', xname='rho', legend=r'$\overline{G}^{\star}_{\nabla \rho}$'),

        # PlotDataCsv(r, 175, 'gave', xname='rho'),
        # PlotDataCsv(r, 176, 'gave_shat', xname='rho'),
        # PlotDataCsv(r, 177, 'gave_shat_gxi', xname='rho'),

        # PlotDataCsv(r, 169, 'gave', xname='rho', legend=r'$\overline{G}^\star$'),
        # PlotDataCsv(r, 170, 'gave_shat', xname='rho', legend=r'$\overline{G}^{\star}_\kappa$'),
        # PlotDataCsv(r, 171, 'gave_shat_gxi', xname='rho', legend=r'$\overline{G}^{\star}_{\nabla \rho}$'),

        # PlotDataCsv(r, 1270, 'xte2ETGM', xname='rho', legend='Old'),
        # PlotDataCsv(r, 1282, 'xte2ETGM', xname='rho', legend='New'),

        # PlotDataCsv(r, 1304, 'xte2ETGM', xname='rho', legend=r'$\chi_\mathrm{e, max}$'),


        # MAX vs SUMMED MODES
        # PlotDataCsv(r, 1304, 'xteETGM', xname='rho', legend=r'$\chi_\mathrm{e, max}$'),
        # PlotDataCsv(r, 1431, 'xteETGM', xname='rho', ymult=0.05, legend=r'$\frac{1}{20}\chi_\mathrm{e, sum}$'),
        # PlotDataCsv(r, 1324, 'xte2ETGM', xname='rho', legend=r'$\chi^*_\mathrm{e, max}$'),
        # PlotDataCsv(r, 1431, 'xte2ETGM', xname='rho', ymult=0.1666, legend=r'$\frac{1}{6}\chi^*_\mathrm{e, sum}$'),

        # PlotDataCsv(r, 1652, 'Apara2ETGM', xname='rho',),
        # PlotDataCsv(r, 1660, 'Apara2ETGM', xname='rho',),

        # G_ave change in xte2, min kyrhos = 0.2
        # PlotDataCsv(r, 1648, 'xte2ETGM', xname='rho', legend=r'using $\omega_\mathrm{De} / \overline{G}$'),  # MAX
        # PlotDataCsv(r, 1652, 'xte2ETGM', xname='rho', legend=r'using $\omega_\mathrm{De}$'),
        # PlotDataCsv(r, 1655, 'xte2ETGM', xname='rho', legend=r'using $\omega_\mathrm{De} / \overline{G}$'),  # SUMMED
        # PlotDataCsv(r, 1656, 'xte2ETGM', xname='rho', legend=r'using $\omega_\mathrm{De}$'),
        # PlotDataCsv(r, 1655, 'xteETGM', xname='rho', legend=r'$\chi_\mathrm{e}$'),  # SUMMED
        # PlotDataCsv(r, 1656, 'xte2ETGM', xname='rho', ymult=20, legend=r'$\chi^*_\mathrm{e} \cdot 20$'),
        # PlotDataCsv(r, 1657, 'xte2ETGM', xname='rho', ymult=20, legend=r'$\chi^*_\mathrm{e} \cdot 20\,(g_\mathrm{Bu} = 1)$'),


        # KYRHOS VS SATURATION EXPO
        # PlotDataCsv(r, 1312, 'xteETGM', xname='rho', legend='n=50'),

        # CONVERGENCE SPEED (sat expo = 2)
        # PlotDataCsv(r, 1431, 'xte2ETGM', xname='rho', legend=r'Sum (n=$10^5$, Exponential)'),
        # PlotDataCsv(r, 1401, 'xte2ETGM', xname='rho', legend=r'n=$10^5$, Linear'),
        # PlotDataCsv(r, 1430, 'xte2ETGM', xname='rho', legend=r'n=$10^3$, Linear'),
        # PlotDataCsv(r, 1432, 'xte2ETGM', xname='rho', legend=r'n=$50,\,\,$ Exponential'),

        # CONVERGENCE SPEED (sat expo = 1)
        # PlotDataCsv(r, 1435, 'xte2ETGM', xname='rho', legend=r'n=$10^5$, Linear'),
        # PlotDataCsv(r, 1436, 'xte2ETGM', xname='rho', legend=r'n=$10^3$, Linear'),
        # PlotDataCsv(r, 1437, 'xte2ETGM', xname='rho', legend=r'n=$50,\,\,$ Exponential'),

        # SINGLE KYRHOS VALUES vs SUM (sat expo = 2)
        # PlotDataCsv(r, 1431, 'xteETGM', xname='rho', legend=r'Summed Modes'),
        # PlotDataCsv(r, 1373, 'xteETGM', xname='rho', ymult=2, legend=r'$\mathrm{cf} = 2,\,\,\, k_\mathrm{y}\rho_\mathrm{s}=0.2$'),
        # PlotDataCsv(r, 1427, 'xteETGM', xname='rho', ymult=10, legend=r'$\mathrm{cf} = 10, k_\mathrm{y}\rho_\mathrm{s}=0.5$'),
        # PlotDataCsv(r, 1374, 'xteETGM', xname='rho', ymult=20, legend=r'$\mathrm{cf} = 20, k_\mathrm{y}\rho_\mathrm{s}=1$'),
        # PlotDataCsv(r, 1373, 'xte2ETGM', xname='rho', ymult=0.1, legend=r'$\mathrm{cf} = 0.1, k_\mathrm{y}\rho_\mathrm{s}=0.2$'),
        # PlotDataCsv(r, 1427, 'xte2ETGM', xname='rho', ymult=1.4, legend=r'$\mathrm{cf} = 1.4, k_\mathrm{y}\rho_\mathrm{s}=0.5$'),
        # PlotDataCsv(r, 1374, 'xte2ETGM', xname='rho', ymult=10, legend=r'$\mathrm{cf} = 10,\,\, k_\mathrm{y}\rho_\mathrm{s}=1$'),

        # SINGLE KYRHOS VALUES vs SUM (sat expo = 1)
        # PlotDataCsv(r, 1435, 'xteETGM', xname='rho', legend=r'Summed Modes'),  # n=$10^5$ (Linear)
        # PlotDataCsv(r, 1438, 'xteETGM', xname='rho', ymult=4, legend=r'$\mathrm{cf} = 4, k_\mathrm{y}\rho_\mathrm{s}=0.2$'),
        # PlotDataCsv(r, 1439, 'xteETGM', xname='rho', ymult=5, legend=r'$\mathrm{cf} = 5, k_\mathrm{y}\rho_\mathrm{s}=0.5$'),
        # PlotDataCsv(r, 1440, 'xteETGM', xname='rho', ymult=8, legend=r'$\mathrm{cf} = 8, k_\mathrm{y}\rho_\mathrm{s}=1$'),
        # PlotDataCsv(r, 1438, 'xte2ETGM', xname='rho', ymult=0.25, legend=r'$\mathrm{cf} = 0.25, k_\mathrm{y}\rho_\mathrm{s}=0.2$'),
        # PlotDataCsv(r, 1439, 'xte2ETGM', xname='rho', ymult=1.0, legend=r'$\mathrm{cf} = 1.00, k_\mathrm{y}\rho_\mathrm{s}=0.5$'),
        # PlotDataCsv(r, 1440, 'xte2ETGM', xname='rho', ymult=5.0, legend=r'$\mathrm{cf} = 5.00, k_\mathrm{y}\rho_\mathrm{s}=1$'),

        # SINGLE KYRHOS VALUES (sat expo = 2)
        # PlotDataCsv(r, 1373, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s}=0.2$'),
        # PlotDataCsv(r, 1374, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s}=1$'),
        # PlotDataCsv(r, 1375, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s}=5$'),

        # SINGLE KYRHOS VALUES (sat expo = 1)
        # PlotDataCsv(r, 1370, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s}=0.2$'),
        # PlotDataCsv(r, 1371, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s}=1$'),
        # PlotDataCsv(r, 1372, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s}=5$'),

        # SUMMED MODE CONVERGENCE (Range)
        # PlotDataCsv(r, 1334, 'xte2ETGM', xname='rho', legend='n=10000 (Base)'),
        # PlotDataCsv(r, 1382, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 40]$'),
        # PlotDataCsv(r, 1383, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 30]$'),
        # PlotDataCsv(r, 1384, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 20]$'),
        # PlotDataCsv(r, 1392, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.5, 50]$'),
        # PlotDataCsv(r, 1393, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [1.0, 50]$'),
        # PlotDataCsv(r, 1394, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [2.0, 50]$'),

        # SUMMED MODE CONVERGENCE (Exponential vs Linear)
        # PlotDataCsv(r, 1401, 'xte2ETGM', xname='rho', legend=r'n=$10^5$ (Linear)'),
        # PlotDataCsv(r, 1408, 'xte2ETGM', xname='rho', legend=r'n=$50\,\,\,$ (Linear)'),
        # PlotDataCsv(r, 1402, 'xte2ETGM', xname='rho', legend=r'n=$100$ (Linear)'),
        # PlotDataCsv(r, 1403, 'xte2ETGM', xname='rho', legend=r'n=$200$ (Linear)'),
        # PlotDataCsv(r, 1409, 'xteETGM', xname='rho', legend=r'n=$50\,\,\,$ (Exponential)'),
        # PlotDataCsv(r, 1405, 'xteETGM', xname='rho', legend=r'n=$100$ (Exponential)'),
        # PlotDataCsv(r, 1406, 'xteETGM', xname='rho', legend=r'n=$200$ (Exponential)'),

        # CALIBRATION (m=2)
        # PlotDataCsv(r, 1478, 'xteETGM', xname='rho', legend=r'$\chi_\mathrm{e,max}$, cf=2e-3'),
        # PlotDataCsv(r, 1478, 'xte2ETGM', xname='rho', legend=r'$\chi^*_\mathrm{e,max}$, cf=4e-3'),
        # PlotDataCsv(r, 1480, 'xteETGM', xname='rho', legend=r'$\chi_\mathrm{e,sum}$, cf=1e-4'),
        # PlotDataCsv(r, 1480, 'xte2ETGM', xname='rho', legend=r'$\chi^*_\mathrm{e,sum}$, cf=1e-3'),

        # MAX MODE CONVERGENCE (Exponential vs Linear)
        # PlotDataCsv(r, 1410, 'xteETGM', xname='rho', legend=r'n=$10^5$ (Linear)'),
        # PlotDataCsv(r, 1411, 'xte2ETGM', xname='rho', legend=r'n=$50\,\,\,$ (Linear)'),
        # PlotDataCsv(r, 1412, 'xte2ETGM', xname='rho', legend=r'n=$100$ (Linear)'),
        # PlotDataCsv(r, 1413, 'xte2ETGM', xname='rho', legend=r'n=$200$ (Linear)'),
        # PlotDataCsv(r, 1414, 'xteETGM', xname='rho', legend=r'n=$50\,\,\,$ (Exponential)'),
        # PlotDataCsv(r, 1415, 'xteETGM', xname='rho', legend=r'n=$100$ (Exponential)'),
        # PlotDataCsv(r, 1416, 'xteETGM', xname='rho', legend=r'n=$200$ (Exponential)'),

        # LOWER KYRHOS RANGE (SUM m=2)
        # PlotDataCsv(r, 1441, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 50]$'),
        # PlotDataCsv(r, 1442, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.5, 50]$'),
        # PlotDataCsv(r, 1443, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [1.0, 50]$'),
        # PlotDataCsv(r, 1444, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [2.0, 50]$'),

        # LOWER KYRHOS RANGE (SUM m=1)
        # PlotDataCsv(r, 1445, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 50]$'),
        # PlotDataCsv(r, 1446, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.5, 50]$'),
        # PlotDataCsv(r, 1449, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [1.0, 50]$'),
        # PlotDataCsv(r, 1450, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [2.0, 50]$'),

        # UPPER KYRHOS RANGE (SUM m=2)
        # PlotDataCsv(r, 1441, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 50]$'),
        # PlotDataCsv(r, 1451, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 40]$'),
        # PlotDataCsv(r, 1452, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 30]$'),
        # PlotDataCsv(r, 1453, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 20]$'),

        # UPPER KYRHOS RANGE (SUM m=1)
        # PlotDataCsv(r, 1445, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 50]$'),
        # PlotDataCsv(r, 1454, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 40]$'),
        # PlotDataCsv(r, 1455, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 30]$'),
        # PlotDataCsv(r, 1456, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 20]$'),

        # LOWER KYRHOS RANGE (MAX m=2)
        # PlotDataCsv(r, 1457, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 50]$'),
        # PlotDataCsv(r, 1458, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.5, 50]$'),
        # PlotDataCsv(r, 1459, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [1.0, 50]$'),
        # PlotDataCsv(r, 1460, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [2.0, 50]$'),

        # LOWER KYRHOS RANGE (MAX m=1)
        # PlotDataCsv(r, 1461, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 50]$'),
        # PlotDataCsv(r, 1462, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.5, 50]$'),
        # PlotDataCsv(r, 1463, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [1.0, 50]$'),
        # PlotDataCsv(r, 1464, 'xteETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [2.0, 50]$'),

        # UPPER KYRHOS RANGE (MAX m=2)
        # PlotDataCsv(r, 1457, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 50]$'),
        # PlotDataCsv(r, 1465, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 40]$'),
        # PlotDataCsv(r, 1466, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 30]$'),
        # PlotDataCsv(r, 1467, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 20]$'),

        # UPPER KYRHOS RANGE (MAX m=1)
        # PlotDataCsv(r, 1461, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 50]$'),
        # PlotDataCsv(r, 1468, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 40]$'),
        # PlotDataCsv(r, 1469, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 30]$'),
        # PlotDataCsv(r, 1470, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{y}\rho_\mathrm{s} \in [0.2, 20]$'),

        # wexb EFFECTS (SUM)
        # PlotDataCsv(r, 1441, 'xteETGM', xname='rho', legend=r'$\omega_{\mathrm{E}\zvals\mathrm{B}}$ Off'),  # m=2
        # PlotDataCsv(r, 1471, 'xteETGM', xname='rho', legend=r'$\omega_{\mathrm{E}\zvals\mathrm{B}}$ On'),   # m=2
        # PlotDataCsv(r, 1445, 'xteETGM', xname='rho', legend=r'$\omega_{\mathrm{E}\zvals\mathrm{B}}$ Off'),  # m=1
        # PlotDataCsv(r, 1472, 'xteETGM', xname='rho', legend=r'$\omega_{\mathrm{E}\zvals\mathrm{B}}$ On'),   # m=1

        # wexb EFFECTS (MAX)
        # PlotDataCsv(r, 1457, 'xte2ETGM', xname='rho', legend=r'$\omega_{\mathrm{E}\zvals\mathrm{B}}$ Off'),  # m=2
        # PlotDataCsv(r, 1473, 'xte2ETGM', xname='rho', legend=r'$\omega_{\mathrm{E}\zvals\mathrm{B}}$ On'),   # m=2
        # PlotDataCsv(r, 1461, 'xte2ETGM', xname='rho', legend=r'$\omega_{\mathrm{E}\zvals\mathrm{B}}$ Off'),  # m=1
        # PlotDataCsv(r, 1474, 'xte2ETGM', xname='rho', legend=r'$\omega_{\mathrm{E}\zvals\mathrm{B}}$ On'),   # m=1

        # MAX MODE CONVERGENCE (xte or xte2)
        # PlotDataCsv(r, 1336, 'xte2ETGM', xname='rho', legend='n=50'),
        # PlotDataCsv(r, 1337, 'xte2ETGM', xname='rho', legend='n=100'),
        # PlotDataCsv(r, 1338, 'xte2ETGM', xname='rho', legend='n=200'),
        # PlotDataCsv(r, 1339, 'xte2ETGM', xname='rho', legend='n=400'),

        # PlotDataCsv(r, 1401, 'xte2ETGM', xname='rho', legend='n=100k (Linear)'),
        # PlotDataCsv(r, 1400, 'xte2ETGM', xname='rho', legend='n=200'),

        # SUM METHOD VERIFICATION
        # PlotDataCsv(r, 1417, 'xteETGM', xname='rho', legend='n=50'),
        # PlotDataCsv(r, 1418, 'xteETGM', xname='rho', legend='n=100'),
        # PlotDataCsv(r, 1417, 'xte2ETGM', xname='rho', legend=r'$(k_\mathrm{y}\rho_\mathrm{s})_\mathrm{max}=50$'),
        # PlotDataCsv(r, 1419, 'xte2ETGM', xname='rho', legend=r'$(k_\mathrm{y}\rho_\mathrm{s})_\mathrm{max}=100$'),
        # PlotDataCsv(r, 1417, 'xte2ETGM', xname='rho', legend=r'n=50, $\,\,\,(k_\mathrm{y}\rho_\mathrm{s})_\mathrm{max}=50$'),
        # PlotDataCsv(r, 1420, 'xte2ETGM', xname='rho', legend=r'n=100, $(k_\mathrm{y}\rho_\mathrm{s})_\mathrm{max}=100$'),

        # SATURATION EXPONENT
        # PlotDataCsv(r, 1431, 'xte2ETGM', xname='rho', legend='m = 2'),  # summed modes
        # PlotDataCsv(r, 1435, 'xte2ETGM', xname='rho', legend='m = 1'),  # summed modes
        # PlotDataCsv(r, 1367, 'xte2ETGM', xname='rho', legend='m = 2'),  # max mode
        # PlotDataCsv(r, 1369, 'xte2ETGM', xname='rho', legend='m = 1'),  # max mode

        # OUTPUT COMPARE
        # PlotDataCsv(r, 1270, 'gmaETGM', xname='rho'),
        # PlotDataCsv(r, 1270, 'omgETGM', xname='rho'),
        # PlotDataCsv(r, 1270, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 1282, 'xte2ETGM', xname='rho'),
        # PlotDataCsv(r, 1270, 'kyrhosETGM', xname='rho'),
        # PlotDataCsv(r, 1282, 'satETGM', xname='rho'),

        # PlotDataCsv(r, 1363, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 1363, 'xte2ETGM', xname='rho'),  # gbu = 1
        # PlotDataCsv(r, 1364, 'xteETGM', xname='rho'),
        # PlotDataCsv(r, 1364, 'xte2ETGM', xname='rho'),

        # ELECTROSTATIC LIMIT
        # PlotDataCsv(r, 1270, 'gmaETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1281, 'gmaETGM', xname='rho', legend='Electrostatic'),
        # PlotDataCsv(r, 1270, 'omgETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1281, 'omgETGM', xname='rho', legend='Electrostatic'),
        # PlotDataCsv(r, 1270, 'xteETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1281, 'xteETGM', xname='rho', legend='Electrostatic'),
        # PlotDataCsv(r, 1270, 'xte2ETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1281, 'xte2ETGM', xname='rho', legend='Electrostatic'),
        # PlotDataCsv(r, 1270, 'kyrhosETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1281, 'kyrhosETGM', xname='rho', legend='Electrostatic'),

        # PlotDataCsv(r, 1270, 'gmaETGM', xname='rho', legend=r'using $B_\mathrm{unit}$'),
        # PlotDataCsv(r, 1271, 'gmaETGM', xname='rho', legend=r'using $B_\phi$'),
        # PlotDataCsv(r, 1270, 'omgETGM', xname='rho', legend=r'using $B_\mathrm{unit}$'),
        # PlotDataCsv(r, 1271, 'omgETGM', xname='rho', legend=r'using $B_\phi$'),
        # PlotDataCsv(r, 1270, 'xteETGM', xname='rho', legend=r'using $B_\mathrm{unit}$'),
        # PlotDataCsv(r, 1271, 'xteETGM', xname='rho', legend=r'using $B_\phi$'),
        # PlotDataCsv(r, 1270, 'xte2ETGM', xname='rho', legend=r'using $B_\mathrm{unit}$'),
        # PlotDataCsv(r, 1271, 'xte2ETGM', xname='rho', legend=r'using $B_\phi$'),
        # PlotDataCsv(r, 1270, 'kyrhosETGM', xname='rho', legend=r'using $B_\mathrm{unit}$'),
        # PlotDataCsv(r, 1271, 'kyrhosETGM', xname='rho', legend=r'using $B_\phi$'),

        # PlotDataCsv(r, 1270, 'gmaETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.10$'),
        # PlotDataCsv(r, 1275, 'gmaETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.25$'),
        # PlotDataCsv(r, 1276, 'gmaETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.50$'),
        # PlotDataCsv(r, 1270, 'omgETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.10$'),
        # PlotDataCsv(r, 1275, 'omgETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.25$'),
        # PlotDataCsv(r, 1276, 'omgETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.50$'),
        # PlotDataCsv(r, 1270, 'xteETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.10$'),
        # PlotDataCsv(r, 1275, 'xteETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.25$'),
        # PlotDataCsv(r, 1276, 'xteETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.50$'),
        # PlotDataCsv(r, 1270, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.10$'),
        # PlotDataCsv(r, 1275, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.25$'),
        # PlotDataCsv(r, 1276, 'xte2ETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.50$'),
        # PlotDataCsv(r, 1270, 'kyrhosETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.10$'),
        # PlotDataCsv(r, 1275, 'kyrhosETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.25$'),
        # PlotDataCsv(r, 1276, 'kyrhosETGM', xname='rho', legend=r'$k_\mathrm{x}/k_\mathrm{y} = 0.50$'),

        # PlotDataCsv(r, 1270, 'gmaETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s}$ scan'),
        # PlotDataCsv(r, 1272, 'gmaETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s} = 5$'),
        # PlotDataCsv(r, 1273, 'gmaETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s} = 10$'),
        # PlotDataCsv(r, 1270, 'omgETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s}$ scan'),
        # PlotDataCsv(r, 1272, 'omgETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s} = 5$'),
        # PlotDataCsv(r, 1273, 'omgETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s} = 10$'),
        # PlotDataCsv(r, 1270, 'xteETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s}$ scan'),
        # PlotDataCsv(r, 1272, 'xteETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s} = 5$'),
        # PlotDataCsv(r, 1273, 'xteETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s} = 10$'),
        # PlotDataCsv(r, 1270, 'xte2ETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s}$ scan'),
        # PlotDataCsv(r, 1272, 'xte2ETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s} = 5$'),
        # PlotDataCsv(r, 1273, 'xte2ETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s} = 10$'),
        # PlotDataCsv(r, 1270, 'kyrhosETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s}$ scan'),
        # PlotDataCsv(r, 1272, 'kyrhosETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s} = 5$'),
        # PlotDataCsv(r, 1273, 'kyrhosETGM', xname='rho', legend=r'$k_y\rho_\mathrm{s} = 10$'),

        # PlotDataCsv(r, 1270, 'gmaETGM', xname='rho', legend=r'$\omega_{E \zvals B}$ off'),
        # PlotDataCsv(r, 1277, 'gmaETGM', xname='rho', legend=r'$\omega_{E \zvals B}$ on'),
        # PlotDataCsv(r, 1270, 'omgETGM', xname='rho', legend=r'$\omega_{E \zvals B}$ off'),
        # PlotDataCsv(r, 1277, 'omgETGM', xname='rho', legend=r'$\omega_{E \zvals B}$ on'),
        # PlotDataCsv(r, 1270, 'xteETGM', xname='rho', legend=r'$\omega_{E \zvals B}$ off'),
        # PlotDataCsv(r, 1277, 'xteETGM', xname='rho', legend=r'$\omega_{E \zvals B}$ on'),
        # PlotDataCsv(r, 1270, 'xte2ETGM', xname='rho', legend=r'$\omega_{E \zvals B}$ off'),
        # PlotDataCsv(r, 1277, 'xte2ETGM', xname='rho', legend=r'$\omega_{E \zvals B}$ on'),
        # PlotDataCsv(r, 1270, 'kyrhosETGM', xname='rho', legend=r'$\omega_{E \zvals B}$ off'),
        # PlotDataCsv(r, 1277, 'kyrhosETGM', xname='rho', legend=r'$\omega_{E \zvals B}$ on'),

        # PlotDataCsv(r, 1270, 'gmaETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1278, 'gmaETGM', xname='rho', legend=r'$g_{\mathrm{ne}} = 0$'),
        # PlotDataCsv(r, 1270, 'omgETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1278, 'omgETGM', xname='rho', legend=r'$g_{\mathrm{ne}} = 0$'),
        # PlotDataCsv(r, 1270, 'xteETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1278, 'xteETGM', xname='rho', legend=r'$g_{\mathrm{ne}} = 0$'),
        # PlotDataCsv(r, 1270, 'xte2ETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1278, 'xte2ETGM', xname='rho', legend=r'$g_{\mathrm{ne}} = 0$'),
        # PlotDataCsv(r, 1270, 'kyrhosETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1278, 'kyrhosETGM', xname='rho', legend=r'$g_{\mathrm{ne}} = 0$'),

        # PlotDataCsv(r, 1270, 'gmaETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1279, 'gmaETGM', xname='rho', legend=r'$g_{\mathrm{Te}} = 0$'),
        # PlotDataCsv(r, 1270, 'omgETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1279, 'omgETGM', xname='rho', legend=r'$g_{\mathrm{Te}} = 0$'),
        # PlotDataCsv(r, 1270, 'xteETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1279, 'xteETGM', xname='rho', legend=r'$g_{\mathrm{Te}} = 0$'),
        # PlotDataCsv(r, 1270, 'xte2ETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1279, 'xte2ETGM', xname='rho', legend=r'$g_{\mathrm{Te}} = 0$'),
        # PlotDataCsv(r, 1270, 'kyrhosETGM', xname='rho', legend='Default'),
        # PlotDataCsv(r, 1279, 'kyrhosETGM', xname='rho', legend=r'$g_{\mathrm{Te}} = 0$'),

        # PlotDataCsv(r, 1270, 'gmaETGM', xname='rho', legend='Collisional'),
        # PlotDataCsv(r, 1280, 'gmaETGM', xname='rho', legend='Collisionless'),
        # PlotDataCsv(r, 1270, 'omgETGM', xname='rho', legend='Collisional'),
        # PlotDataCsv(r, 1280, 'omgETGM', xname='rho', legend='Collisionless'),
        # PlotDataCsv(r, 1270, 'xteETGM', xname='rho', legend='Collisional'),
        # PlotDataCsv(r, 1280, 'xteETGM', xname='rho', legend='Collisionless'),
        # PlotDataCsv(r, 1270, 'xte2ETGM', xname='rho', legend='Collisional'),
        # PlotDataCsv(r, 1280, 'xte2ETGM', xname='rho', legend='Collisionless'),
        # PlotDataCsv(r, 1270, 'kyrhosETGM', xname='rho', legend='Collisional'),
        # PlotDataCsv(r, 1280, 'kyrhosETGM', xname='rho', legend='Collisionless'),

        # PlotDataCsv(r, 457, 'gmaETGM', xname='rho', legend=r'0$\beta^\prime$'),
        # PlotDataCsv(r, 450, 'gmaETGM', xname='rho', legend=r'1$\beta^\prime$'),
        # PlotDataCsv(r, 469, 'gmaETGM', xname='rho', legend=r'2$\beta^\prime$'),
        # PlotDataCsv(r, 470, 'gmaETGM', xname='rho', legend=r'3$\beta^\prime$'),
        # PlotDataCsv(r, 450, 'omgETGM', xname='rho', legend=r'with $\beta^\prime$'),
        # PlotDataCsv(r, 457, 'omgETGM', xname='rho', legend=r'without $\beta^\prime$'),
        # PlotDataCsv(r, 450, 'xteETGM', xname='rho', legend=r'with $\beta^\prime$'),
        # PlotDataCsv(r, 457, 'xteETGM', xname='rho', legend=r'without $\beta^\prime$'),
        # PlotDataCsv(r, 450, 'xte2ETGM', xname='rho', legend=r'with $\beta^\prime$'),
        # PlotDataCsv(r, 457, 'xte2ETGM', xname='rho', legend=r'without $\beta^\prime$'),
        # PlotDataCsv(r, 450, 'kyrhoeETGM', xname='rho', legend=r'with $\beta^\prime$'),
        # PlotDataCsv(r, 457, 'kyrhoeETGM', xname='rho', legend=r'without $\beta^\prime$'),
        # PlotDataCsv(r, 450, 'kyrhosETGM', xname='rho', legend=r'with $\beta^\prime$'),
        # PlotDataCsv(r, 457, 'kyrhosETGM', xname='rho', legend=r'without $\beta^\prime$'),

        # PlotDataCsv(r, 1270, 'gmaETGM', xname='rho', legend=r'using $\hat{s}_{\nabla\rho}$'),
        # PlotDataCsv(r, 458, 'gmaETGM', xname='rho', legend=r'using $s$'),
        # PlotDataCsv(r, 1270, 'omgETGM', xname='rho', legend=r'using $\hat{s}_{\nabla\rho}$'),
        # PlotDataCsv(r, 458, 'omgETGM', xname='rho', legend=r'using $s$'),
        # PlotDataCsv(r, 1270, 'xteETGM', xname='rho', legend=r'using $\hat{s}_{\nabla\rho}$'),
        # PlotDataCsv(r, 458, 'xteETGM', xname='rho', legend=r'using $s$'),
        # PlotDataCsv(r, 1270, 'xte2ETGM', xname='rho', legend=r'using $\hat{s}_{\nabla\rho}$'),
        # PlotDataCsv(r, 458, 'xte2ETGM', xname='rho', legend=r'using $s$'),
        # PlotDataCsv(r, 1270, 'kyrhoeETGM', xname='rho', legend=r'using $\hat{s}_{\nabla\rho}$'),
        # PlotDataCsv(r, 458, 'kyrhoeETGM', xname='rho', legend=r'using $s$'),
        # PlotDataCsv(r, 1270, 'kyrhosETGM', xname='rho', legend=r'using $\hat{s}_{\nabla\rho}$'),
        # PlotDataCsv(r, 458, 'kyrhosETGM', xname='rho', legend=r'using $s$'),

        # kryhos scan compare
        # PlotDataCsv(r, 460, 'gmaETGM', xname='rho', legend=r'20'),
        # PlotDataCsv(r, 461, 'gmaETGM', xname='rho', legend=r'100'),
        # PlotDataCsv(r, 450, 'gmaETGM', xname='rho', legend=r'500'),
        # PlotDataCsv(r, 460, 'omgETGM', xname='rho', legend=r'20'),
        # PlotDataCsv(r, 461, 'omgETGM', xname='rho', legend=r'100'),
        # PlotDataCsv(r, 450, 'omgETGM', xname='rho', legend=r'500'),
        # PlotDataCsv(r, 460, 'xteETGM', xname='rho', legend=r'20'),
        # PlotDataCsv(r, 461, 'xteETGM', xname='rho', legend=r'100'),
        # PlotDataCsv(r, 450, 'xteETGM', xname='rho', legend=r'500'),
        # PlotDataCsv(r, 460, 'xte2ETGM', xname='rho', legend=r'20'),
        # PlotDataCsv(r, 461, 'xte2ETGM', xname='rho', legend=r'100'),
        # PlotDataCsv(r, 450, 'xte2ETGM', xname='rho', legend=r'500'),
        # PlotDataCsv(r, 460, 'kyrhoeETGM', xname='rho', legend=r'20'),
        # PlotDataCsv(r, 461, 'kyrhoeETGM', xname='rho', legend=r'100'),
        # PlotDataCsv(r, 450, 'kyrhoeETGM', xname='rho', legend=r'500'),
        # PlotDataCsv(r, 460, 'kyrhosETGM', xname='rho', legend=r'20'),
        # PlotDataCsv(r, 461, 'kyrhosETGM', xname='rho', legend=r'100'),
        # PlotDataCsv(r, 450, 'kyrhosETGM', xname='rho', legend=r'500'),

        # kryhoe scan compare
        # PlotDataCsv(r, 450, 'gmaETGM', xname='rho', legend=r'20'),
        # PlotDataCsv(r, 462, 'gmaETGM', xname='rho', legend=r'500'),
        # PlotDataCsv(r, 450, 'omgETGM', xname='rho', legend=r'20'),
        # PlotDataCsv(r, 462, 'omgETGM', xname='rho', legend=r'500'),
        # PlotDataCsv(r, 450, 'xteETGM', xname='rho', legend=r'20'),
        # PlotDataCsv(r, 462, 'xteETGM', xname='rho', legend=r'500'),
        # PlotDataCsv(r, 450, 'xte2ETGM', xname='rho', legend=r'20'),
        # PlotDataCsv(r, 462, 'xte2ETGM', xname='rho', legend=r'500'),
        # PlotDataCsv(r, 450, 'kyrhoeETGM', xname='rho', legend=r'20'),
        # PlotDataCsv(r, 462, 'kyrhoeETGM', xname='rho', legend=r'500'),
        # PlotDataCsv(r, 450, 'kyrhosETGM', xname='rho', legend=r'20'),
        # PlotDataCsv(r, 462, 'kyrhosETGM', xname='rho', legend=r'500'),

        # PlotDataCsv(r, 298, 'xteETGM', xname='rho', legend=r'no scan, $k_y\rho_s = 0.33$'),
        # PlotDataCsv(r, 347, 'xteETGM', xname='rho', legend=r'normalized using $k_y\rho_s$ from scan'),
        # PlotDataCsv(r, 348, 'xteETGM', xname='rho', legend=r'normalized using $k_y\rho_s$ = 0.33'),

        # PlotDataCsv(r, 347, 'gmaETGM', xname='rho', legend=r'normalized using $k_y\rho_s$ from scan'),
        # PlotDataCsv(r, 345, 'gmaETGM', xname='rho', legend=r'normalized using $k_y\rho_s$ = 0.33'),
        # PlotDataCsv(r, 346, 'gmaETGM', xname='rho', legend=r'normalized using $k_y\rho_s$ = 0.33'),

        # PlotDataCsv(r, 340, 'xteETGM', xname='rho', legend=r'$k_y\rho_s = 0.06$'),
        # PlotDataCsv(r, 341, 'xteETGM', xname='rho', legend=r'$k_y\rho_s = 0.07$'),
        # PlotDataCsv(r, 342, 'xteETGM', xname='rho', legend=r'$k_y\rho_s = 0.08$'),
        # PlotDataCsv(r, 343, 'xteETGM', xname='rho', legend=r'$k_y\rho_s = 0.09$'),

        # PlotDataCsv(r, 139, 'gmaETGM', xname='shear', rho_value=0.2),
        # PlotDataCsv(r, 103, 'gmaETGM', xname='nuei', rho_value=0.20, legend='nuei scan'),
        # PlotDataCsv(r, 106, 'gmaETGM', xname='nuei', rho_value=0.20, legend='te scan'),
        # PlotDataCsv(r, 122, 'gmaETGM', xname='nuei', rho_value=0.20, legend='ne scan'),
        # PlotDataCsv(r, 103, 'alphamhd', xname='nuei', rho_value=0.20, legend='nuei scan'),
        # PlotDataCsv(r, 106, 'alphamhd', xname='nuei', rho_value=0.20, legend='te scan'),
        # PlotDataCsv(r, 122, 'alphamhd', xname='nuei', rho_value=0.20, legend='ne scan'),
        # PlotDataCsv(r, 124, 'alphamhd', xname='nuei', rho_value=0.20, legend='ne scan'),
        # PlotDataCsv(r, 39, 'alphamhd', xname='etae', rho_value=0.8, legend='gte  3, gne  1'),
        # PlotDataCsv(r, 36, 'alphamhd', xname='etae', rho_value=0.8, legend='gte  2, gne  0'),
        # PlotDataCsv(r, 38, 'alphamhd', xname='etae', rho_value=0.8, legend='gte  1, gne -1'),
        # PlotDataCsv(r, 37, 'alphamhd', xname='etae', rho_value=0.8, legend='gte  0, gne -2'),
        # PlotDataCsv(r, 40, 'alphamhd', xname='etae', rho_value=0.8, legend='gte -1, gne -3'),
        # PlotDataCsv(r, 35, 'gmaETGM', xname='etae', scan_num=35, rho_value=0.10),
        # PlotDataCsv(r, 47, 'gmaETGM', xname='alphamhd', rho_value=0.1, legend=r'$\mathtt{tau\_scan}$'),
        # PlotDataCsv(r, 48, 'gmaETGM', xname='alphamhd', rho_value=0.1, legend=r'$\mathtt{te\_scan}$'),
        # PlotDataCsv(r, 49, 'gmaETGM', xname='alphamhd', rho_value=0.1, legend=r'$\mathtt{ti\_scan}$'),
        # CSV: Growth rate vs Average Magnetic Surface Curvature
        # PlotDataCsv(r, 69, 'gne', xname='alphamhd', rho_value=0.20, runname=''),
        # PlotDataCsv(r, 69, 'gne', xname='alphamhd', rho_value=0.40, runname=''),
        # PlotDataCsv(r, 69, 'gne', xname='alphamhd', rho_value=0.60, runname=''),
        # PlotDataCsv(r, 15, 'gmaETGM', xname='var_to_scan', rho_value=0.39, runname=r'OLD'),
        # PlotDataCsv(r, 55, 'gmaETGM', xname='var_to_scan', rho_value=0.39, runname=r'NEW'),
        # PlotDataCsv(r, 15, 'gmaETGM', xname='var_to_scan', rho_value=0.61, runname=''),
        # PlotDataCsv(r, 15, 'gmaETGM', xname='var_to_scan', rho_value=0.62, runname=''),
        # PlotDataCsv(r, 15, 'gmaETGM', xname='var_to_scan', rho_value=0.7, runname=''),
        # PlotDataCsv(r, 21, 'gmaETGM', xname='var_to_scan', rho_value=0.94, runname=''),
        # PlotDataCsv(r, 3, 'gmaETGM', xname='var_to_scan', rho_value=0.3),
        # PlotDataCsv(r, 3, 'alphamhd', xname='var_to_scan', rho_value=0.4),
    )

    main(fig_data)
