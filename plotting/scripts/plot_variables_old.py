#!/usr/bin/python3

"""Calls plotting/plot_variables.py

Contains all original plotting calls for
"""

# Standard Packages
import sys; sys.path.insert(0, '../'); sys.path.insert(0, '../../')
import logging

# 3rd Party Packages
import matplotlib.pyplot as plt
from PyQt5.QtGui import QImage
from PyQt5.QtWidgets import QApplication

# Local Packages
import modules.options  # needed, even though linter doesn't like it
from plotting.plot_variables import FigData, PlotDataCsv, PlotDataCdf, main
import modules.utils as utils
from modules.enums import SaveType
from plotting.modules.plotstyles import PlotStyles, StyleType


_log = logging.getLogger(__name__)


if __name__ == '__main__':
    """Run this module directly to plot variable data stored in CDF or CSV files"""

    utils.init_logging()

    # Initialize visual styles for the plot
    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.AIP3,
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
        'legend.fontsize': 9,
        'legend.handlelength': 1.8,  # length of lines
        'figure.dpi': 200,
    })

    # Define settings for the plot
    fig_data = FigData(
        replace_offset_text=0,
        allow_title_runid=1,
        allow_title_time=1,
        allow_title_factor=1,
        allow_title_rho=1,
        invert_y_axis=False,
        invert_x_axis=False,
        nomralize_y_axis=False,
        nomralize_x_axis=False,
        title_override=' ',
        ylabel_override='',
        xlabel_override='',
        savefig=False,
        savedata=False,
    )

    # Define data for the plot (Examples shown below)
    fig_data.set(
        # CDF: Same y-variable, different x-variables
        # PlotDataCdf(runid='138536A01', yname='te', xname='rmina', zval=0.50),
        # PlotDataCdf(runid='138536A01', yname='te', xname='rho', zval=0.50),
        # CDF: Different y-variable units
        # PlotDataCdf(runid='138536A01', yname='te', xname='rho', zval=0.50),
        # PlotDataCdf(runid='138536A01', yname='ti', xname='rho', zval=0.50),
        # PlotDataCdf(runid='138536A01', yname='btor', xname='rho', zval=0.50),
        # CDF: Different runid, y-variables, and zvals
        # PlotDataCdf(runid='120982A09', yname='ne', xname='rho', zval=0.60),
        # PlotDataCdf(runid='120968A02', yname='ni', xname='rho', zval=0.50),
        # PlotDataCdf(runid='129041A10', yname='nd', xname='rho', zval=0.40),
        # CDF: Compare as a function of time
        # PlotDataCdf(runid='138536A01', yname='te', xname='time', zval=0.40, timeplot=True),
        # PlotDataCdf(runid='138536A01', yname='te', xname='time', zval=0.50, timeplot=True),
        # PlotDataCdf(runid='138536A01', yname='te', xname='time', zval=0.60, timeplot=True),
        # CDF: Compare TRANSP and MMM calculations (must be defined in calculations.py)
        # PlotDataCdf(runid='141552A01', yname='loge', xname='rho', zval=0.629, source='cdf'),
        # PlotDataCdf(runid='141552A01', yname='loge', xname='rho', zval=0.629),
        # CDF and CSV: Compare same variable from different data sources
        # PlotDataCdf(runid='138536A01', yname='ne', xname='rho', zval=0.629),
        # PlotDataCsv(runid='138536A01', yname='ne', xname='rho', scan_num=2),
        # CSV: Different scanned variables with same scan factor
        # PlotDataCsv(runid='138536A01', yname='ne', xname='rho', scan_num=1, scan_factor=2.5),
        # PlotDataCsv(runid='138536A01', yname='ne', xname='rho', scan_num=2, scan_factor=2.5),
        # PlotDataCsv(runid='138536A01', yname='ne', xname='rho', scan_num=5, scan_factor=2.5),
        # CSV: Comparing output results
        # PlotDataCsv(runid='138536A01', yname='xteETGM', xname='rho', scan_num=17, legend='exbs = 0'),
        # PlotDataCsv(runid='138536A01', yname='xteETGM', xname='rho', scan_num=54, legend='exbs = 1'),
        # # CSV: Growth rate vs \eta_e
        # PlotDataCsv(runid='138536A01', yname='gmaETGM', xname='etae', scan_num=35, rho_value=0.10),
        # CSV: Growth rate as a function of different scanned variables
        # PlotDataCsv(runid='138536A01', yname='gmaETGM', xname='var_to_scan', scan_num=15, rho_value=0.39, runname=r'OLD'),
        # PlotDataCsv(runid='138536A01', yname='gmaETGM', xname='var_to_scan', scan_num=55, rho_value=0.39, runname=r'NEW'),
    )

    fig_data.set(

        # PlotDataCdf(runid='129016A03', yname='te', xname='rho', zval=0.395, legend='Experiment'),
        # PlotDataCdf(runid='129016Q50', yname='te', xname='rho', zval=0.395, legend='Off'),
        # PlotDataCdf(runid='129016Q93', yname='te', xname='rho', zval=0.395, legend='W20 +'),
        # PlotDataCdf(runid='129016Q94', yname='te', xname='rho', zval=0.395, legend='W20 +-'),
        # PlotDataCdf(runid='129016Q50', yname='ti', xname='rho', zval=0.395, legend='Off'),
        # PlotDataCdf(runid='129016W47', yname='ti', xname='rho', zval=0.395, legend='DBM Failure'),
        # PlotDataCdf(runid='129016W52', yname='ti', xname='rho', zval=0.395, legend='DBM Failure'),
        # PlotDataCdf(runid='129016X31', yname='xkew20', xname='rho', zval=0.4, legend='DBM Failure'),

        # PlotDataCdf(runid='129016X35', yname='bpol', xname='xb', zval=0.391, source='mmm', legend='MMM'),
        # PlotDataCdf(runid='129016X31', yname='bpol', xname='xb', zval=0.391, source='cdf', legend='CDF'),
        # xmax=0.79,

        # PlotDataCdf(runid='129016X31', yname='dbp', xname='xb', zval=0.391, source='mmm', legend='MMM Definition'),
        # PlotDataCdf(runid='129016X31', yname='d2bp', xname='xb', zval=0.391, source='mmm', legend='MMM Definition'),



        # PlotDataCdf(runid='129016Q93', yname='xkemtm', xname='rho', zval=0.4),
        # PlotDataCdf(runid='129016Q93', yname='wexb', xname='rho', zval=0.4),

        # PlotDataCdf(runid='129016A03', yname='gte', xname='rmina', zval=0.46),
        # PlotDataCdf(runid='129016A03', yname='gte', xname='rmina', zval=0.462),
        # PlotDataCdf(runid='129016A03', yname='gte', xname='rmina', zval=0.465),
        # PlotDataCdf(runid='129016A03', yname='gte', xname='rmina', zval=0.49),
        # xmin=0.59,xmax=0.61,
        
        # PlotDataCdf(runid='138536A01', yname='vtorx', xname='xb', zval=0.629, legend=r'RMAJ$\cdot$OMEGA'),
        # PlotDataCdf(runid='138536A01', yname='vtorxnc', xname='xb', zval=0.629, legend=r'VTORX$\_$NC'),

        # PlotDataCdf(runid='120982A09', yname='ti', xname='xb', zval=0.629),
        # PlotDataCdf(runid='120982A09', yname='nf', xname='xb', zval=0.629),
        # PlotDataCdf(runid='120982A09', yname='pf', xname='xb', zval=0.629),

        # PlotDataCdf(runid='138536A01', yname='gti', xname='xb', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='gnf', xname='xb', zval=0.629),
        # PlotDataCdf(runid='120982A09', yname='gpf', xname='xb', zval=0.629),
        # PlotDataCdf(runid='120982A09', yname='gpf2', xname='xb', zval=0.629),
        

        # PlotDataCdf(runid='138536A01', yname='nf', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='nf2', xname='rho', zval=0.629),

        # PlotDataCdf(runid='121123K55', yname='tmhdf', xname='rho', zval=0.629),

        # PlotDataCdf(runid='138536A01', yname='pmhdf', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='p', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='pmhd', xname='rho', zval=0.629),
        # # PlotDataCdf(runid='138536A01', yname='pf', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='pr', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='pmhdr', xname='rho', zval=0.629),

        # PlotDataCdf(runid='18696R06', yname='pmhd', xname='rho', zval=0.629),
        # PlotDataCdf(runid='38265R80', yname='pmhd', xname='rho', zval=0.629),
        # PlotDataCdf(runid='50000A10', yname='pmhd', xname='rho', zval=0.629),
        # PlotDataCdf(runid='80200A13', yname='pmhd', xname='rho', zval=0.629),
        # PlotDataCdf(runid='80300A02', yname='pmhd', xname='rho', zval=0.629),


        # PlotDataCdf(runid='138536A01', yname='ti', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='epapl', xname='rho', zval=0.629, source='cdf'),
        # PlotDataCdf(runid='138536A01', yname='epapp', xname='rho', zval=0.629, source='cdf'),
        # PlotDataCdf(runid='138536A01', yname='epa', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='tmhdf', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='pmhdf', xname='rho', zval=0.629),

        # PlotDataCdf(runid='138536A01', yname='nf', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='nfd', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='nfmp', xname='rho', zval=0.629),

        # PlotDataCdf(runid='138536A01', yname='nfmp', xname='rmajm', zval=0.629),

        # PlotDataCdf(runid='138536A01', yname='ebeamsum', xname='rmajm', zval=0.629),

        # PlotDataCdf(runid='138536A01', yname='dvol', xname='rho', zval=0.629),

        # PlotDataCdf(runid='138536A01', yname='rmaj', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='rmajm', xname='rho', zval=0.629),

        # PlotDataCdf(runid='138536A01', yname='ebeam', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='ebeampp', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='ebeampl', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='ebeamsum', xname='rho', zval=0.629),

        # PlotDataCdf(runid='138536A01', yname='tf', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='tmhdf', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='tbtbe', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='tfast', xname='rho', zval=0.629),

        # PlotDataCdf(runid='138536A01', yname='bpshi', xname='time', zval=0.629),

        # PlotDataCdf(runid='138536A01', yname='gr2i', xname='rho', zval=0.629),

        # PlotDataCdf(runid='138536A01', yname='pmhdf', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='p', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='ebeamr', xname='rho', zval=0.58),

        # PlotDataCdf(runid='138536A01', yname='tf', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='tmhdf', xname='rho', zval=0.629),
        # PlotDataCdf(runid='138536A01', yname='tbtbe', xname='rho', zval=0.629),

        # PlotDataCdf(runid='138536A01', yname='ebeam', xname='rho', zval=0.0),
        # PlotDataCdf(runid='138536A01', yname='tf', xname='rho', zval=0.8),
        # PlotDataCdf(runid='138536A01', yname='gtf', xname='rho', zval=0.8),

        # PlotDataCdf(runid='132017T01', yname='pmhd', xname='rho', zval=0.629),
        # PlotDataCdf(runid='132017T01', yname='pmhdf', xname='rho', zval=0.629),
        # PlotDataCdf(runid='132017T01', yname='pmhdr', xname='rho', zval=0.629),
        # PlotDataCdf(runid='132017T01', yname='p', xname='rho', zval=0.629),
        # PlotDataCdf(runid='132017T01', yname='pf', xname='rho', zval=0.629),

        # PlotDataCdf(runid='141031A01', yname='ebeam', xname='rho', zval=0.629),
        # PlotDataCdf(runid='141031A01', yname='epa', xname='rho', zval=0.629),
        # PlotDataCdf(runid='141031A01', yname='tmhdf', xname='rho', zval=0.629),


        # PlotDataCdf(runid='141552A01', yname='ebeam2', xname='rho', zval=0.629),
        # PlotDataCdf(runid='141552A01', yname='tmhdf', xname='rho', zval=0.629),

        # PlotDataCdf(runid='129017A04', yname='epapl', xname='rho', zval=0.629),
        # PlotDataCdf(runid='129017A04', yname='epapp', xname='rho', zval=0.629),
        # PlotDataCdf(runid='129017A04', yname='epa', xname='rho', zval=0.629),
        # PlotDataCdf(runid='129017A04', yname='ebeam', xname='rho', zval=0.629),
        



        # PlotDataCdf(runid='129016A03', yname='ne', xname='rho', zval=0.46),
        # PlotDataCdf(runid='129016A03', yname='ni', xname='rho', zval=0.46),
        # PlotDataCdf(runid='129016A03', yname='ni2', xname='rho', zval=0.46),

        # PlotDataCdf(runid='138536A01', yname='gq', xname='rho', input_points=1001, apply_smoothing=True, zval=0.63),
        # PlotDataCdf(runid='129041A10', yname='gq', xname='rho', input_points=1001, apply_smoothing=True, zval=0.49),
        # PlotDataCdf(runid='120982A09', yname='gq', xname='rho', input_points=1001, apply_smoothing=True, zval=0.62),
        # PlotDataCdf(runid='120968A02', yname='gq', xname='rho', input_points=1001, apply_smoothing=True, zval=0.56),

        # PlotDataCdf(runid='129041A10', yname='gq', xname='rho', apply_smoothing=True, zval=0.49),
        # PlotDataCdf(runid='129041A10', yname='gq', xname='rho', input_points=101, apply_smoothing=True, zval=0.49),
        # PlotDataCdf(runid='129041A10', yname='gq', xname='rho', input_points=1001, apply_smoothing=True, zval=0.49),
        
        # PlotDataCdf(runid='85126T02', yname='curdoh', zval=4),
        # PlotDataCdf(runid='85610T01', yname='curdoh', zval=4),
        # PlotDataCdf(runid='85126T02', yname='area', zval=4),
        # PlotDataCdf(runid='85610T01', yname='area', zval=4),
       
        # PlotDataCdf(runid='85126T02', yname='curoh', xname='time', timeplot=True, zval=0),
        # PlotDataCdf(runid='85610T01', yname='curoh', xname='time', timeplot=True, zval=0),
        # PlotDataCdf(runid='85126T02', yname='curohrho', xname='rho', zval=4),
        # PlotDataCdf(runid='85610T01', yname='curohrho', xname='rho', zval=4),
        
        # PlotDataCdf(runid='120968A02', yname='icur', xname='time', timeplot=True, zval=0),
        # PlotDataCdf(runid='120968A02', yname='jcur', xname='rho', zval=0.559),
        # PlotDataCdf(runid='120968A02', yname='area', xname='rho', zval=0.559),

        # PlotDataCdf(runid='120968A02', yname='q', xname='rho', zval=0.56),
        # PlotDataCdf(runid='120982A09', yname='wexb', xname='rho', zval=0.62),
        # PlotDataCdf(runid='129041A10', yname='q', xname='rho', zval=0.49),
        # PlotDataCdf(runid='138536A01', yname='gnh', xname='rho', zval=0.629, apply_smoothing=1, input_points=101),
        # PlotDataCsv(runid='138536A01', yname='q', xname='rho', scan_num=253),
        # PlotDataCdf(runid='138536A01', yname='etanc', xname='rho', zval=0.629),


        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=547, legend=r'Old', scan_factor=0.001),
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=551, legend=r'New', scan_factor=0.001),
        # PlotDataCsv(runid='138536A01', yname='kyrhosEPM', xname='rho', scan_num=655, legend=r'Default'),
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=656, legend=r'Kyrhos Correlation Loop'),
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=657, legend=r'Kyrhos Correlation Loop'),
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=443),
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=444),
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=449),

        # PlotDataCsv(runid='129016A03', yname='omgnETGM', xname='kyrhosETGM', scan_num=72, rho_value=0.477, legend=r'BR34 = (1 + FL) / kpc'),
        # PlotDataCsv(runid='129016A03', yname='omgnETGM', xname='kyrhosETGM', scan_num=71, rho_value=0.477, legend=r'BR34 = 1 / kpc'),
        # PlotDataCsv(runid='129016A03', yname='gmanETGM', xname='kyrhosETGM', scan_num=37, rho_value=0.477, legend=r'1.0 fsa'),
        # PlotDataCsv(runid='129016A03', yname='gmanETGM', xname='kyrhosETGM', scan_num=42, rho_value=0.477, legend=r'0.5 fsa'),
        # title_override=r'$r/a = 0.6$',xmin=np.log10(0.1),
        # PlotDataCsv(runid='129016A03', yname='omgnETGM', xname='kyrhosETGM', scan_num=37, rho_value=0.477, legend=r'0.50 kx/ky'),
        # PlotDataCsv(runid='129016A03', yname='gmanETGM', xname='kyrhosETGM', scan_num=43, rho_value=0.477, legend=r'0.25 kx/ky'),
        # title_override=r'$r/a = 0.6$',xmin=np.log10(0.1),
        # PlotDataCsv(runid='129016A03', yname='gmanETGM', xname='kyrhosETGM', scan_num=37, rho_value=0.572, legend=r'1.0 fsa'),
        # PlotDataCsv(runid='129016A03', yname='gmanETGM', xname='kyrhosETGM', scan_num=42, rho_value=0.572, legend=r'0.5 fsa'),
        # title_override=r'$r/a = 0.7$',xmin=np.log10(0.1),

        # PlotDataCsv(runid='129016A03', yname='gmaDBM', xname='kyrhosDBM, scan_num=93, rho_value=0.57),

        # PlotDataCsv(runid='129016A03', yname='gmanETGM', xname='kyrhosETGM', scan_num=42, rho_value=0.57),
        # PlotDataCsv(runid='129016A03', yname='omgnETGM', xname='kyrhosETGM', scan_num=42, rho_value=0.57),
        # title_override=r'$r/a = 0.7$',xmin=np.log10(0.1),

        # PlotDataCsv(runid='129016A03', yname='fte', xname='gte', scan_num=78, rho_value=0.48),
        # PlotDataCsv(runid='129016A03', yname='fte', xname='gte', scan_num=78, rho_value=0.57),

        # PlotDataCsv(runid='129016A03', yname='gmaETGM', xname='gte', scan_num=83, rho_value=0.57),  # kyrhos = 13
        # PlotDataCsv(runid='129016A03', yname='gmaETGM', xname='gte', scan_num=84, rho_value=0.57),  # kyrhos = 18
        # PlotDataCsv(runid='129016A03', yname='gmaETGM', xname='gte', scan_num=85, rho_value=0.57),  # kyrhos = 25
        # PlotDataCsv(runid='129016A03', yname='gmaETGM', xname='gte', scan_num=86, rho_value=0.57),  # kyrhos = 33

        # PlotDataCsv(runid='129016A03', yname='gmaETGM', xname='kyrhosETGM', scan_num=87, rho_value=0.48),  # kyrhos scan
        # PlotDataCsv(runid='129016A03', yname='gmaETGM', xname='kyrhosETGM', scan_num=87, rho_value=0.57),  # kyrhos scan

        # PlotDataCsv(runid='129016A03', yname='gmanETGM', xname='gte', scan_num=88, rho_value=0.48, legend=r'$k_{\rm y}\rho_{\rm s} = 2$'),  # kyrhos = 2
        # PlotDataCsv(runid='129016A03', yname='gmanETGM', xname='gte', scan_num=89, rho_value=0.48, legend=r'$k_{\rm y}\rho_{\rm s} = 4$'),  # kyrhos = 4
        # PlotDataCsv(runid='129016A03', yname='gmanETGM', xname='gte', scan_num=90, rho_value=0.48, legend=r'$k_{\rm y}\rho_{\rm s} = 8$'),  # kyrhos = 8
        # PlotDataCsv(runid='129016A03', yname='gmanETGM', xname='gte', scan_num=91, rho_value=0.48, legend=r'$k_{\rm y}\rho_{\rm s} = 16$'),  # kyrhos = 16

        
        # PlotDataCsv(runid='129016A03', yname='vcz', xname='rho', scan_num=75,),
        # PlotDataCsv(runid='129016A03', yname='vct', xname='rho', scan_num=75,),
        # PlotDataCsv(runid='129016A03', yname='vcp', xname='rho', scan_num=75,),

        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=751, legend='Guess = 0'),
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=752, legend='Default Guess'),
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=750, legend='Increased Guess'),
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=753, legend='Enormous Guess'),


        # PlotDataCsv(runid='138536A01', yname='gaveEPM', xname='rho', scan_num=757, legend='Guess = 0'),
        # PlotDataCsv(runid='138536A01', yname='gaveEPM', xname='rho', scan_num=756, legend='Default Guess'),
        # PlotDataCsv(runid='138536A01', yname='gaveEPM', xname='rho', scan_num=755, legend='Increased Guess'),
       
        # PlotDataCsv(runid='138536A01', yname='gmaW20i', xname='rho', scan_num=1048, legend='Guess = 0'),
        # PlotDataCsv(runid='138536A01', yname='gmaW20i', xname='rho', scan_num=1051, legend='Guess = 0'),
        

        # PlotDataCsv(runid='138536A01', yname='fti', xname='rho', scan_num=1130, scan_factor=0.7, legend='BaseNew'),
        # PlotDataCsv(runid='138536A01', yname='fti', xname='rho', scan_num=1248, scan_factor=0.7, legend='New'),
        # PlotDataCsv(runid='138536A01', yname='fti', xname='rho', scan_num=1244, scan_factor=0.7, legend='BaseOld'),
        # PlotDataCsv(runid='138536A01', yname='fti', xname='rho', scan_num=1241, scan_factor=0.7, legend='Old'),
        # PlotDataCsv(runid='138536A01', yname='fte', xname='rho', scan_num=1203, scan_factor=0.8, legend='++'),
        

        # PlotDataCsv(runid='138536A01', yname='gmaW20i', xname='weiland_gmult', scan_num=1059, rho_value=0.99, legend='1 + i'),
        # PlotDataCsv(runid='138536A01', yname='gmaW20i', xname='weiland_gmult', scan_num=1059, rho_value=0.99, legend='1 + i'),


        # PlotDataCsv(runid='118341T54', yname='omgW20i', xname='weiland_gmult', scan_num=565, rho_value=0.5, legend='Default'),
        # PlotDataCsv(runid='118341T54', yname='omgW20i', xname='weiland_gmult', scan_num=564, rho_value=0.5, legend='1 + i'),
        # xmin=-1, ymax=2.5e5,

        # PlotDataCsv(runid='118341T54', yname='kyrhosW20e', xname='rho', scan_num=564),

        # PlotDataCsv(runid='118341T54', yname='kyrhosW20i', xname='rho', scan_num=564),
        # ymin=0.19,

        # PlotDataCsv(runid='138536A01', yname='kyrhosW20i', xname='rho', scan_num=1058),
        # PlotDataCsv(runid='138536A01', yname='kyrhosW20e', xname='rho', scan_num=1058),

        
        # PlotDataCsv(runid='118341T54', yname='gmaW20i', xname='weiland_gmult', scan_num=564, rho_value=0.6),
        # PlotDataCsv(runid='118341T54', yname='gmaW20e', xname='weiland_gmult', scan_num=564, rho_value=0.5),
        # PlotDataCsv(runid='138536A01', yname='gmaW20e', xname='weiland_gmult', scan_num=1056, rho_value=0.88),
        
        # # 21908, 10, 18.14
        # PlotDataCsv(runid='138536A01', yname='omgEPM', xname='rho', scan_num=786, legend='Previous Matching'),
        # # 8728, 1, 7.07
        # PlotDataCsv(runid='138536A01', yname='omgEPM', xname='rho', scan_num=802, legend='Updated Matching'),
        # title_override='Default Guess',

        # # 25506, 17, 22.08
        # PlotDataCsv(runid='138536A01', yname='omgEPM', xname='rho', scan_num=789, legend='Previous Matching'),
        # # 13841, 6, 11.03
        # PlotDataCsv(runid='138536A01', yname='omgEPM', xname='rho', scan_num=803, legend='Updated Matching'),
        # title_override='Larger Guess',

        # # # 8879, 6, 6.88
        # PlotDataCsv(runid='138536A01', yname='omgEPM', xname='rho', scan_num=791, legend='Previous Matching'),
        # # 3352, 0, 2.31
        # PlotDataCsv(runid='138536A01', yname='omgEPM', xname='rho', scan_num=805, legend='Updated Matching'),
        # title_override='Guess = 0',

        # # 41414, 24, 26.38
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=795, legend='++G, dN'),
        # # 20455, 0, 12.01
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=796, legend='++G, 0N'),
        # # 20364, 0, 11.84
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=804, legend='++G, 0N, R/A'),
          
        # # 19793, 7, 16.65
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=797, legend='-G, dN'),
        # # 12149, 8, 8.26
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=798, legend='-G, 0N'),
        # # 7294, 2, 5.94
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=806, legend='-G, 0N, R/A'),
          
        # ------------------ 20% previous guess and 80% new guess below

         # # 23668, 15
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=768, legend='default Norm'),
        # # 21206, 13
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=769, legend='fixed Norm'),
        # # 15534, 14
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=770, legend='no Norm'),
       
        # # 16668, 10
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=773, legend='default Norm'),
        # # 16640, 9
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=772, legend='fixed Norm'),
        # # 8919, 4
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=771, legend='no Norm'),

        # # 20881, 17 -- BIGGER GUESS
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=774, legend='Lrg, default Norm'),
        # # 25047, 11, 14.7
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=774, legend='ELrg, default Norm'),
        

        # 21150, 17
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=775, legend='Lrg, fixed Norm'),
        # # 26601, 11, 15.44
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=775, legend='ELrg, fixed Norm'),
        # 8078, 0
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=776, legend='Lrg, no Norm'),
        # # 13979, 1, 7.91 -- HUGE GUESS
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=777, legend='ELrg, no Norm'),

        # # 16668, 10, 13.32 -- DEFAULT GUESS
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=780, legend='def, default Norm'),
        # # 3139, 0, 2.14
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=785, legend='0, default Norm'),
        
        # # 16640, 9. 13.13
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=781, legend='def, fixed Norm'),
        # # 8010, 4, 5.97
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=784, legend='0, fixed Norm'),


        # # 8919, 4, 6.66
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=782, legend='def, no Norm'),
        # # 8156, 6, 6.27 -- NO GUESS
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=783, legend='0, no Norm'),

        # PlotDataCsv(runid='138536A01', yname='gaveEPM', xname='rho', scan_num=766, legend='Default Guess'),
        # PlotDataCsv(runid='138536A01', yname='gaveEPM', xname='rho', scan_num=765, legend='Increased Guess'),
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=759, legend='0.1 Increased Guess'),
        # PlotDataCsv(runid='138536A01', yname='gaveEPM', xname='rho', scan_num=754, legend='Enormous Guess'),
        # PlotDataCsv(runid='120982A09', yname='gmaW20i', xname='rho', scan_num=12006, legend='cubic'),
        # PlotDataCsv(runid='120982A09', yname='gmaW20i', xname='rho', scan_num=12011, legend='quadratic'),

        # PlotDataCsv(runid='120968A02', yname='kyrhosDBM', xname='rho', scan_num=48, legend='100 loops'),
        # PlotDataCsv(runid='120968A02', yname='kyrhosDBM', xname='rho', scan_num=51, legend='10,000 loops'),
        # PlotDataCsv(runid='120982A09', yname='xteETG', xname='rho', scan_num=12012, legend='quadratic'),

        # PlotDataCsv(runid='129041A10', yname='xteDBM', xname='rho', scan_num=115, legend='exp 10'),
        # PlotDataCsv(runid='129041A10', yname='xteDBM', xname='rho', scan_num=116, legend='exp 20'),
        # PlotDataCsv(runid='129041A10', yname='xteDBM', xname='rho', scan_num=117, legend='exp 50'),
        # PlotDataCsv(runid='129041A10', yname='xteDBM', xname='rho', scan_num=118, legend='exp 100'),

        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=513, legend='converged'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=500, legend='linear 250'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=501, legend='linear 500'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=502, legend='linear 1000'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=504, legend='linear 5000'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=513, legend='converged'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=503, legend=2000),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=505, legend='exp 250'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=506, legend='exp 500'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=507, legend='exp 1000'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=513, legend='converged'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=510, legend='exp 50'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=511, legend='exp 100'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=512, legend='exp 200'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=513, legend='linear 100k'),

        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=513, legend='converged'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=509, legend='10000'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=502, legend='linear 1000'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=503, legend='linear 2000'),
        # PlotDataCsv(runid='129041A10', yname='xteMTM', xname='rho', scan_num=511, legend='exp 100'),

        # PlotDataCsv(runid='120982A09', yname='fvp', xname='rho', scan_num=38, legend=r'$\ \ \,$VPOL'),
        # PlotDataCsv(runid='120982A09', yname='fvp', xname='rho', scan_num=39, legend=r'$-$VPOL'),

        # PlotDataCsv(runid='129016A03', yname='xte', xname='rmina', scan_num=6),
        # PlotDataCsv(runid='129016A04', yname='xte', xname='rmina', scan_num=18),
        # PlotDataCsv(runid='16325T10', yname='xtiW20', xname='rho', scan_num=13200, scan_factor=9.5),
        # PlotDataCsv(runid='16325T10', yname='xtiW20', xname='rho', scan_num=13201, scan_factor=9.5),
        
        # PlotDataCsv(runid='16325T10', yname='xtiW20', xname='rho', scan_num=49,), #original, ifort, s
        # PlotDataCsv(runid='16325T10', yname='xtiW20', xname='rho', scan_num=50,), #opt, ifort, s
        # PlotDataCsv(runid='16325T10', yname='xtiW20', xname='rho', scan_num=51,), #original, ifort, nos
        # PlotDataCsv(runid='16325T10', yname='xtiW20', xname='rho', scan_num=52,), #opt, ifort, nos
        # PlotDataCsv(runid='16325T10', yname='xtiW20', xname='rho', scan_num=54,), #opt2, ifort, nos

        # PlotDataCsv(runid='120968A02', yname='gmaDBM', xname='rho', scan_num=13216,),
        # PlotDataCsv(runid='120968A02', yname='gmaDBM', xname='rho', scan_num=13214,),
        # PlotDataCsv(runid='120968A02', yname='gmaDBM', xname='rho', scan_num=52,),

 
        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=9,  legend=r'$n_{_\mathrm{L}} = 1000$'),
        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=10, legend=r'$n_{_\mathrm{L}} = 2000$'),
        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=11, legend=r'$n_{_\mathrm{L}} = 4000$'),
        # title_override=r'Linear $k_y\rho_s$ Increments',
        # allow_title_runid=0,
        # allow_title_time=0,
        # ymax=9,

        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=12, legend=r'$n_{_\mathrm{E}} = 50$'),
        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=13, legend=r'$n_{_\mathrm{E}} = 100$'),
        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=14, legend=r'$n_{_\mathrm{E}} = 200$'),
        # title_override=r'Exponential $k_y\rho_s$ Increments',
        # allow_title_runid=0,
        # allow_title_time=0,
        # ymax=9,

        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=8, legend=r'Converged'),
        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=10, legend=r'$n_{_\mathrm{L}} = 2000$'),
        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=13, legend=r'$n_{_\mathrm{E}} = 100\,$'),
        # title_override=r'Linear vs. Exponential',
        # allow_title_runid=0,
        # allow_title_time=0,
        # ymin=6.2, ymax=8.1, xmin=0.8, xmax=0.93,

        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=10, legend=r'$n = 100\ \,$ (Exp.)'),
        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=11, legend=r'$n = 200\ \,$ (Exp.)'),

        # # NEW NSTX COLLISIONALITY COMPARISON (9.0.4)
        # PlotDataCsv(runid='120968A02', yname='gmaMTM', xname='rmina', scan_num=13300, legend=r'High'),
        # PlotDataCsv(runid='120982A09', yname='gmaMTM', xname='rmina', scan_num=13300, legend=r'Low'),
        # xmax=0.8,

   
        # PlotDataCdf(runid='138536W03', yname='ti', xname='rho', zval=0.755, legend=r'MMM'),
        # PlotDataCdf(runid='138536A01', yname='ti', legend=r'Analysis', zval=0.755, source='mmm'),
        # xmax=1,

        # PlotDataCdf(runid='138536W03', yname='xki', xname='rho', zval=0.755, legend=r'$\chi_{\rm i, tot}$'),
        # PlotDataCdf(runid='138536W03', yname='xkimmm', xname='rho', zval=0.755, legend=r'$\chi_{\rm i, MMM}$'),
        # PlotDataCdf(runid='138536W03', yname='xkemmm', xname='rho', zval=0.755, legend=r'$\chi_{\rm e, MMM}$'),
        # xmax=1, ymin=0, ymax=35,


        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=10, legend=r'$\mathtt{Exp.}$ 200'),
        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=11, legend=r'$\mathtt{Lin.}$ 4000'),
        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=12, legend=r'$\mathtt{Converged}$'),
        # PlotDataCsv(runid='121123K55', yname='xteMTM', xname='rho', scan_num=13, legend=r'$\mathtt{Converged}$'),
        
        # PlotDataCsv(runid='138536A01', yname='xteETGM', xname='rho',  scan_num=49),
        # PlotDataCsv(runid='138536A01', yname='xte2ETGM', xname='rho', scan_num=49),
        # PlotDataCsv(runid='138536A01', yname='xteETGM', xname='rho',  scan_num=81),
        # PlotDataCsv(runid='138536A01', yname='xte2ETGM', xname='rho', scan_num=81),
        # PlotDataCsv(runid='138536A01', yname='xteETGM', xname='rho',  scan_num=84),
        # PlotDataCsv(runid='138536A01', yname='xte2ETGM', xname='rho', scan_num=84),

        # PlotDataCsv(runid='129016A03', yname='gmaETGM', xname='kyrhosETGM', scan_num=42, rho_value=0.6),
        # # PlotDataCsv(runid='138536A01', yname='gmaETGM', xname='rho', scan_num=99),

        # PlotDataCsv(runid='120968A02', yname='xteETGM', xname='rho',  scan_num=40),
        # PlotDataCsv(runid='120968A02', yname='xteETGM', xname='rho', scan_num=43),


        # PlotDataCsv(runid='138536A01', yname='gmaETGM', xname='ai', scan_num=304, rho_value=0.3, legend=r'kyrhos scan'),
        # PlotDataCsv(runid='138536A01', yname='gmaETGM', xname='ai', scan_num=305, rho_value=0.3, legend=r'kyrhos = 8'),
        # PlotDataCsv(runid='138536A01', yname='gmaETGM', xname='ai', scan_num=306, rho_value=0.3, legend=r'kyrhos = 16'),
        # PlotDataCsv(runid='138536A01', yname='gmaETGM', xname='ai', scan_num=307, rho_value=0.3),



        # PlotDataCsv(runid='138536A01', yname='xtiDBM', xname='rho',  scan_num=5004),
        # PlotDataCsv(runid='138536A01', yname='xti2DBM', xname='rho', scan_num=5004),
        # PlotDataCsv(runid='138536A01', yname='xteDBM', xname='rho',  scan_num=5004),
        # PlotDataCsv(runid='138536A01', yname='xte2DBM', xname='rho', scan_num=5004),
        # PlotDataCsv(runid='138536A01', yname='xdeDBM', xname='rho',  scan_num=5009),
        # PlotDataCsv(runid='138536A01', yname='xde2DBM', xname='rho', scan_num=5009),
        # PlotDataCsv(runid='138536A01', yname='satDBM', xname='rho', scan_num=5006),
        # PlotDataCsv(runid='138536A01', yname='satDBM', xname='rho', scan_num=5007),

        # PlotDataCsv(runid='129041A10', yname='xteDBM', xname='rho',  scan_num=5003),
        # PlotDataCsv(runid='129041A10', yname='xte2DBM', xname='rho', scan_num=5003),
        # PlotDataCsv(runid='129041A10', yname='xtiDBM', xname='rho',  scan_num=5002),
        # PlotDataCsv(runid='129041A10', yname='xti2DBM', xname='rho', scan_num=5002),
        # PlotDataCsv(runid='129041A10', yname='xdeDBM', xname='rho',  scan_num=5008),
        # PlotDataCsv(runid='129041A10', yname='xde2DBM', xname='rho', scan_num=5008),
        # PlotDataCsv(runid='129041A10', yname='satDBM', xname='rho', scan_num=5006),
        # PlotDataCsv(runid='129041A10', yname='satDBM', xname='rho', scan_num=5007),


        # PlotDataCsv(runid='120982A09', yname='xtiDBM', xname='rho',  scan_num=5003),
        # PlotDataCsv(runid='120982A09', yname='xti2DBM', xname='rho', scan_num=5003),
        # PlotDataCsv(runid='120982A09', yname='xteDBM', xname='rho',  scan_num=5003),
        # PlotDataCsv(runid='120982A09', yname='xte2DBM', xname='rho', scan_num=5003),
        # PlotDataCsv(runid='120982A09', yname='xdeDBM', xname='rho',  scan_num=5009),
        # PlotDataCsv(runid='120982A09', yname='xde2DBM', xname='rho', scan_num=5009),

        # PlotDataCsv(runid='120968A02', yname='xteDBM', xname='rho',  scan_num=5002),
        # PlotDataCsv(runid='120968A02', yname='xte2DBM', xname='rho', scan_num=5002),
        # PlotDataCsv(runid='120968A02', yname='xtiDBM', xname='rho',  scan_num=5001),
        # PlotDataCsv(runid='120968A02', yname='xti2DBM', xname='rho', scan_num=5001),
        # PlotDataCsv(runid='120968A02', yname='xdeDBM', xname='rho',  scan_num=5009),
        # PlotDataCsv(runid='120968A02', yname='xde2DBM', xname='rho', scan_num=5009),





        # PlotDataCsv(runid='129041A10', yname='gmaDBM', xname='kyrhosDBM', rho_value=0.6, scan_num=26),
        # PlotDataCsv(runid='129041A10', yname='gmaDBM', xname='kyrhosDBM', rho_value=0.7, scan_num=26),
        # PlotDataCsv(runid='129041A10', yname='gmaDBM', xname='kyrhosDBM', rho_value=0.8, scan_num=26),
        # PlotDataCsv(runid='129041A10', yname='gmaDBM', xname='kyrhosDBM', rho_value=0.9, scan_num=26),

        # PlotDataCsv(runid='120968A02', yname='xti2DBM', xname='rho', scan_num=24, legend='wexb off'),
        # PlotDataCsv(runid='120968A02', yname='xti2DBM', xname='rho', scan_num=25, legend='wexb on'),
        # PlotDataCsv(runid='120968A02', yname='xteDBM', xname='rho', scan_num=19),
        # PlotDataCsv(runid='120968A02', yname='xteDBM', xname='rho', scan_num=17),
        # PlotDataCsv(runid='120968A02', yname='xteDBM', xname='rho', scan_num=17),

        # PlotDataCsv(runid='120968A02', yname='xdeDBM',  xname='rho', scan_num=31),
        # PlotDataCsv(runid='120968A02', yname='xde2DBM', xname='rho', scan_num=31),
        # PlotDataCsv(runid='120982A09', yname='xdeDBM',  xname='rho', scan_num=36),
        # PlotDataCsv(runid='120982A09', yname='xde2DBM', xname='rho', scan_num=36),
        # PlotDataCsv(runid='132411T02', yname='xdeDBM',  xname='rho', scan_num=50),
        # PlotDataCsv(runid='132411T02', yname='xde2DBM', xname='rho', scan_num=50),
        # PlotDataCsv(runid='138536A01', yname='xdeDBM',  xname='rho', scan_num=67),
        # PlotDataCsv(runid='138536A01', yname='xde2DBM', xname='rho', scan_num=67),
        # PlotDataCsv(runid='118341T54', yname='xdeDBM',  xname='rho', scan_num=519),
        # PlotDataCsv(runid='118341T54', yname='xde2DBM', xname='rho', scan_num=519),

        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=338),
        # PlotDataCsv(runid='138536A01', yname='gmaEPM', xname='rho', scan_num=341),

        # PlotDataCsv(runid='132411T02', yname='ti', xname='rho', scan_num=16),
        # PlotDataCsv(runid='132411T02', yname='ti', xname='rho', scan_num=17),

        # PlotDataCsv(runid='132411T02', yname='xtiW20',  xname='rho', scan_num=23),
        # PlotDataCsv(runid='132411T02', yname='xtiW20',  xname='rho', scan_num=24),
        # PlotDataCsv(runid='132411T02', yname='xtiDBM',  xname='rho', scan_num=24),
        # PlotDataCsv(runid='132411T02', yname='xti2DBM', xname='rho', scan_num=24),
        # ymax=4,

        # PlotDataCsv(runid='132411T02', yname='xtiDBM', xname='rho', scan_num=37),
        # PlotDataCsv(runid='132411T02', yname='xtiDBM', xname='rho', scan_num=38),
        # PlotDataCsv(runid='132411T02', yname='xtiDBM', xname='rho', scan_num=39),
        # PlotDataCsv(runid='132411T02', yname='xtiDBM', xname='rho', scan_num=40),
        # PlotDataCsv(runid='132411T02', yname='xtiDBM', xname='rho', scan_num=41),
        # PlotDataCsv(runid='132411T02', yname='xtiDBM', xname='rho', scan_num=42),
        # PlotDataCsv(runid='132411T02', yname='gmaDBM',  xname='rho', scan_num=43),

        # PlotDataCsv(runid='132411T02', yname='xteW20',  xname='rho', scan_num=24),
        # PlotDataCsv(runid='132411T02', yname='xteDBM',  xname='rho', scan_num=24),
        # PlotDataCsv(runid='132411T02', yname='xte2DBM', xname='rho', scan_num=24),
        # PlotDataCsv(runid='132411T02', yname='xteETG',  xname='rho', scan_num=24),
        # PlotDataCsv(runid='132411T02', yname='xteETGM', xname='rho', scan_num=24),
        # PlotDataCdf(runid='132411T02', yname='condewnc', xname='rho',zval=0.56),
        # ymax=10,

        # PlotDataCsv(runid='132411T02', yname='xteDBM', xname='rho', scan_num=9),
        # PlotDataCsv(runid='132411T02', yname='xteW20', xname='rho', scan_num=9),
        # PlotDataCsv(runid='132411T02', yname='xteETGM', xname='rho',scan_num=9),
        # PlotDataCsv(runid='132411T02', yname='xtiDBM', xname='rho',  scan_num=6),
        # PlotDataCsv(runid='132411T02', yname='xti2DBM', xname='rho', scan_num=6),

        # PlotDataCsv(runid='101381T31', yname='gmaDBM', xname='kyrhosDBM', scan_num=9),
        # PlotDataCsv(runid='101381T31', yname='gmaDBM', xname='kyrhosDBM', scan_num=9),
        # PlotDataCsv(runid='101381T31', yname='gmaDBM', xname='kyrhosDBM', scan_num=9),
        # PlotDataCsv(runid='120982A09', yname='gmaDBM', xname='kyrhosDBM', rho_value=0.8, scan_num=31, legend='Ti/Te = 1'),


        # PlotDataCsv(runid='129016A04', yname='xteMTM', scan_num=23),

        # PlotDataCdf(runid='129016Q34', yname='walltime', xname='time', zval=0.5, timeplot=True, legend='08 Processors'),
        # PlotDataCdf(runid='129016Q36', yname='walltime', xname='time', zval=0.5, timeplot=True, legend='16 Processors'),
        # PlotDataCdf(runid='129016Q35', yname='walltime', xname='time', zval=0.5, timeplot=True, legend='32 Processors'),

        # PlotDataCdf(runid='129016Q32', yname='walltime', xname='time', zval=0.5, timeplot=True, legend='Models Disabled'),
        # PlotDataCdf(runid='129016Q33', yname='walltime', xname='time', zval=0.5, timeplot=True, legend='Models Enabled'),

        
        # PlotDataCdf(runid='129016Q50', yname='walltime', xname='time', zval=0.5, timeplot=True, legend='v9.0.1 Enabled'),
        # PlotDataCdf(runid='129016Q51', yname='walltime', xname='time', zval=0.5, timeplot=True, legend='v9.0.1 Disabled'),

        # PlotDataCdf(runid='129016W50', yname='walltime', xname='time', zval=0.5, timeplot=True, legend='v9.0.6 Enabled'),
        # PlotDataCdf(runid='129016W51', yname='walltime', xname='time', zval=0.5, timeplot=True, legend='v9.0.6 Disabled'),
        # xmax=0.35

        # PlotDataCdf(runid='129016Q69', yname='walltime', xname='time', zval=0.5, timeplot=True, legend='v8.2.3 Enabled'),
        # PlotDataCdf(runid='129016Q70', yname='walltime', xname='time', zval=0.5, timeplot=True, legend='v8.2.3 Disabled'),
        
        # PlotDataCdf(runid='129016Q77', yname='walltime', xname='time', zval=0.5, timeplot=True, legend='v8.2.3 Enabled'),
        # PlotDataCdf(runid='129016Q70', yname='walltime', xname='time', zval=0.5, timeplot=True, legend='v8.2.3 Disabled'),

        # PlotDataCdf(runid='129016Q77', yname='vcz', xname='time', zval=0.5, timeplot=True, legend='v8.2.3'),

        # Sat Expo 2 vs 1
        # PlotDataCdf(runid='129016A03', yname='te', xname='rho', zval=0.415, legend='Experiment'),
        # PlotDataCdf(runid='129016Q58', yname='te', xname='rho', zval=0.415, legend=r'Sat Expo = 2'),
        # PlotDataCdf(runid='129016Q67', yname='te', xname='rho', zval=0.415, legend=r'Sat Expo = 1'),

        # PlotDataCdf(runid='129016A03', yname='te', xname='rho', zval=0.3, legend='Experiment'),
        # PlotDataCdf(runid='129016Q51', yname='te', xname='rho', zval=0.3, legend='MMM Off'),
        # PlotDataCdf(runid='129016Q58', yname='te', xname='rho', zval=0.3, legend=r'$\chi_{\rm e, etgm}$'),
        # PlotDataCdf(runid='129016Q59', yname='te', xname='rho', zval=0.3, legend=r'$_{{^\sum}}\chi_{\rm e, etgm}$'),

        # PlotDataCdf(runid='129016A03', yname='te', xname='rho', zval=0.4, legend='Experiment'),
        # PlotDataCdf(runid='129016Q51', yname='te', xname='rho', zval=0.4, legend='MMM Off'),
        # PlotDataCdf(runid='129016Q58', yname='te', xname='rho', zval=0.4, legend=r'$\chi_{\rm e, etgm}$'),
        # PlotDataCdf(runid='129016Q63', yname='te', xname='rho', zval=0.4, legend=r'2$_{{^\sum}}\chi_{\rm e, etgm}$'),

        # PlotDataCdf(runid='129016A03', yname='xkemmm07', xname='rho', zval=0.4, legend='Experiment'),
        # PlotDataCdf(runid='129016Q51', yname='xkemmm07', xname='rho', zval=0.4, legend='MMM Off'),
        # PlotDataCdf(runid='129016Q58', yname='xkemmm07', xname='rho', zval=0.4, legend=r'$\chi_{\rm e, etgm}$'),
        # PlotDataCdf(runid='129016Q63', yname='xkemmm07', xname='rho', zval=0.4, legend=r'2$_{{^\sum}}\chi_{\rm e, etgm}$'),

        # PlotDataCdf(runid='129016A03', yname='te', xname='rho', zval=0.415, legend='Experiment'),
        # PlotDataCdf(runid='129016Q51', yname='te', xname='rho', zval=0.415, legend='MMM Off'),
        # PlotDataCdf(runid='129016Q60', yname='te', xname='rho', zval=0.415, legend=r'$\chi_{\rm e, etgm}^*$'),
        # PlotDataCdf(runid='129016Q61', yname='te', xname='rho', zval=0.415, legend=r'$_{{^\sum}}\chi_{\rm e, etgm}^*$'),

        # PlotDataCdf(runid='129016Q58', yname='xkemmm07', xname='rho', zval=0.3, legend=r'$\chi_{\rm e, etgm}$'),
        # PlotDataCdf(runid='129016Q59', yname='xkemmm07', xname='rho', zval=0.3, legend=r'$_{{^\sum}}\chi_{\rm e, etgm}$'),
        # PlotDataCdf(runid='129016Q60', yname='xkemmm07', xname='rho', zval=0.3, legend=r'$\chi_{\rm e, etgm}^*$'),
        # PlotDataCdf(runid='129016Q61', yname='xkemmm07', xname='rho', zval=0.3, legend=r'$_{{^\sum}}\chi_{\rm e, etgm}^*$'),

        # ETGM VS MTM
        # PlotDataCdf(runid='129016A03', yname='te', xname='rho', zval=0.3, legend='Experiment'),
        # PlotDataCdf(runid='129016Q51', yname='te', xname='rho', zval=0.3, legend='MMM Off'),
        # PlotDataCdf(runid='129016Q58', yname='te', xname='rho', zval=0.3, legend=r'$\chi_{\rm e, etgm}$'),
        # PlotDataCdf(runid='129016Q64', yname='te', xname='rho', zval=0.3, legend=r'$\chi_{\rm e, mtm}$'),

        # MTM OLD vs NEW
        # PlotDataCdf(runid='129016Q69', yname='xkemtm', xname='rho', zval=0.4, legend=r'OLD'),
        # PlotDataCdf(runid='129016Q64', yname='xkemtm', xname='rho', zval=0.4, legend=r'NEW'),

        # MTM OLD vs NEW
        # PlotDataCdf(runid='129016Q75', yname='xkidrbm', xname='xb', zval=0.31, source=r'raw'),
        # PlotDataCdf(runid='129016Q79', yname='xkemmm07', xname='rho', zval=0.4, legend=r'OLD'),
        # PlotDataCdf(runid='129016Q64', yname='xkemtm', xname='rho', zval=0.4, legend=r'NEW'),

        # PlotDataCdf(runid='129016Q64', yname='xkemtm', xname='rho', zval=0.3, legend=r'$\chi_{\rm e, mtm}$'),


        # xteETGM: TRANSP vs Stand Alone with wexb = 0
        # PlotDataCdf(runid='129016Q65', yname='xkemmm07', xname='rho', zval=0.4, legend=r'TRANSP'),
        # PlotDataCsv(runid='129016Q65', yname='xteETGM', xname='rho', scan_num=1, legend='Stand Alone'),

        # # xte2ETGM: TRANSP vs Stand Alone with wexb = 0
        # PlotDataCdf(runid='129016Q66', yname='xkemmm07', xname='rho', zval=0.4, legend=r'TRANSP'),
        # PlotDataCsv(runid='129016Q66', yname='xte2ETGM', xname='rho', scan_num=2, legend='Stand Alone (negative gBu)'),


        # PlotDataCsv(runid='129016Q50', yname='xtiW20', xname='rho', scan_num=12, legend='default'),
        # PlotDataCsv(runid='129016Q50', yname='xtiW20', xname='rho', scan_num=13, legend='btor = 0'),
        # PlotDataCsv(runid='129016Q50', yname='xtiW20', xname='rho', scan_num=14, legend='ne = 1e16'),
        # PlotDataCsv(runid='129016Q50', yname='xtiW20', xname='rho', scan_num=15, legend='both'),


        # PlotDataCsv(runid='129016A04', yname='xtiW20', xname='rho', scan_num=68,),
        # PlotDataCsv(runid='129016A04', yname='xdeW20', xname='rho', scan_num=68,),
        # PlotDataCsv(runid='129016A04', yname='xteW20', xname='rho', scan_num=68,),
        # PlotDataCsv(runid='129016A04', yname='xdz', xname='rho', scan_num=68,),
        # PlotDataCsv(runid='129016A04', yname='xvt', xname='rho', scan_num=68,),
        # PlotDataCsv(runid='129016A04', yname='xvp', xname='rho', scan_num=68,),

        # PlotDataCsv(runid='129016A04', yname='xtiDBM', xname='rho', scan_num=67,),
        # PlotDataCsv(runid='129016A04', yname='xdeDBM', xname='rho', scan_num=67,),
        # PlotDataCsv(runid='129016A04', yname='xteDBM', xname='rho', scan_num=67,),
        # PlotDataCsv(runid='129016A04', yname='xti2DBM', xname='rho', scan_num=68,),
        # PlotDataCsv(runid='129016A04', yname='xde2DBM', xname='rho', scan_num=68,),
        # PlotDataCsv(runid='129016A04', yname='xte2DBM', xname='rho', scan_num=68,),

        # PlotDataCsv(runid='129016A04', yname='xteETGM', xname='rho',  scan_num=68,),
        # PlotDataCsv(runid='129016A04', yname='xte2ETGM', xname='rho', scan_num=68,),
        # PlotDataCsv(runid='129016A04', yname='xteETGM', xname='rho',  scan_num=67,),
        # PlotDataCsv(runid='129016A04', yname='xte2ETGM', xname='rho', scan_num=67,),

        # PlotDataCsv(runid='129016A04', yname='xteMTM', xname='rho',  scan_num=67,),

        # PlotDataCdf(runid='129016A03', yname='ti', xname='rho', zval=0.389),
        # PlotDataCdf(runid='129016Z11', yname='ti', xname='rho', zval=0.389),
        # PlotDataCdf(runid='129016A03', yname='ti', xname='rho', zval=0.39, legend='Experiment'),
        # PlotDataCdf(runid='129016Z11', yname='ti', xname='rho', zval=0.4, legend='MMM On'),
        
        # PlotDataCdf(runid='129016W47', yname='te', xname='rho', zval=0.4, legend='DBM Failure'),
        # PlotDataCdf(runid='129016Q50', yname='ti', xname='rho', zval=0.4, legend='MMM Off'),

        # PlotDataCdf(runid='129016Z18', yname='wexb', xname='rho', zval=0.32, legend='MMM 8in9 '),
        # PlotDataCdf(runid='129016Z19', yname='ti', xname='rho', zval=0.32, legend='MMM 9'),
        # PlotDataCdf(runid='129016Z15', yname='xkemmm', xname='rho', zval=0.34, legend='MMM 8'),
        # PlotDataCdf(runid='129016Z21', yname='xkemmm', xname='rho', zval=0.34, legend='MMM 9 Disabled'),
        # PlotDataCdf(runid='129016Z20', yname='ti', xname='rho', zval=0.31, legend='MMM 8 Disabled'),
        # PlotDataCdf(runid='129016Z24', yname='ti', xname='rho', zval=0.32, legend=r'MMM v8 (Pshare), $+\chi$'),
        # PlotDataCdf(runid='129016Z25', yname='wexb', xname='rho', zval=0.3, legend=r'MMM v8 (Pshare), $\pm\chi$'),
        # PlotDataCdf(runid='129016A03', yname='ti', xname='rho', zval=0.3, legend='Experiment'),
        # title_override='129016, 0.3s',

        # PlotDataCdf(runid='129016Q69', yname='wexb', xname='rho', zval=0.32, legend=r'MMM v8 (Pshare), $+\chi$'), # no pphi
        # PlotDataCdf(runid='129016Z27', yname='te', xname='rho', zval=0.32, legend=r'MMM v8 (Pshare), $+\chi$'), # original TR.DAT
        # PlotDataCdf(runid='129016Z28', yname='te', xname='rho', zval=0.32, legend=r'W20 only (Pshare), $+\chi$'), # original TR.DAT
        # PlotDataCdf(runid='129016A03', yname='ti', xname='rho', zval=0.32, legend='Experiment'),

        # Showing ETG at edge affects wexb
        # PlotDataCdf(runid='129016Z30', yname='xkimmm', xname='rho', zval=0.32, legend=r'Axial + Edge Active + ETG'),
        # PlotDataCdf(runid='129016Z31', yname='xkimmm', xname='rho', zval=0.32, legend=r'Axial + Edge Active'),
        # PlotDataCdf(runid='129016Z33', yname='xkimmm', xname='rho', zval=0.32, legend=r'Confinement Only'), 
        
        # PlotDataCdf(runid='120968A02', yname='bpol', zval=0.5, legend='CDF', source='cdf'),
        # PlotDataCsv(runid='120968A02', yname='bpol', scan_num=74, legend=r'modmmm'),
        # PlotDataCsv(runid='120968A02', yname='bpol', scan_num=75, legend=r'pt$\_$mmm$\_$mod'),

        # PlotDataCdf(runid='120968A02', yname='shear', xname='rho', zval=0.5, legend='CDF', source='cdf'),
        # PlotDataCsv(runid='120968A02', yname='gti', xname='rho', scan_num=78, legend=r'interp, trad'),
        # PlotDataCsv(runid='120968A02', yname='gti', xname='rho', scan_num=79, legend=r'interp, interp'),
        # PlotDataCsv(runid='120968A02', yname='gti', xname='rho', scan_num=80, legend=r'akima, trad'),
        # PlotDataCsv(runid='120968A02', yname='gti', xname='rho', scan_num=81, legend=r'akima, interp'),

        # PlotDataCsv(runid='120968A02', yname='gti', xname='rho', scan_num=78, legend=r'interp1d'),
        # PlotDataCsv(runid='120968A02', yname='gti', xname='rho', scan_num=80, legend=r'Akima'),
        # title_override='120968 (0.559s)',

        # PlotDataCdf(runid='120968W34', yname='xke', legend=r'$\chi_{\mathrm{e, total}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120968W34', yname='condepr', legend=r'$\chi_{\mathrm{e, condepr}}$', zval=0.559, source='cdf'),
        # PlotDataCdf(runid='120968W34', yname='xkemmm',legend=r'$\chi_{\mathrm{e, mmm}}$', zval=0.559, source='cdf'),
        # PlotDataCdf(runid='120968W34', yname='xkeetgm',legend=r'$\chi_{\mathrm{e, etgm}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120968W34', yname='xkepaleo', legend=r'$\chi_{\mathrm{e, paleo}}$', zval=0.559, source='cdf'),
        # PlotDataCdf(runid='120968W34', yname='xkemtm', legend=r'$\chi_{\mathrm{e, mtm}}$', zval=0.559, source='cdf'),
        # # PlotDataCdf(runid='120968W34', yname='xkew20', legend=r'$\chi_{\mathrm{e, w20}}$', zval=0.559, source='cdf'),
        # # PlotDataCdf(runid='120968W34', yname='condewnc', xname='rho', zval=0.559, source='cdf'),
        # title_override='120968', xmin=0.02,

        # PlotDataCdf(runid='120968W34', yname='xki', legend=r'$\chi_{\mathrm{i, total}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120968W34', yname='condipr', legend=r'$\chi_{\mathrm{i, condipr}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120968W34', yname='xkimmm', legend=r'$\chi_{\mathrm{i, mmm}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120968W34', yname='condiwnc', legend=r'$\chi_{\mathrm{i, neo}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120968W34', yname='xkiw20', legend=r'$\chi_{\mathrm{i, w20}}$', zval=0.559, source='cdf'),
        # # PlotDataCdf(runid='120968W34', yname='xkiw20', legend=r'$\chi_{\mathrm{i, w20}}$', zval=0.559, source='cdf'),
        # # PlotDataCdf(runid='120968W34', yname='condewnc', xname='rho', zval=0.559, source='cdf'),
        # title_override='120968', xmin=0.02,

        # PlotDataCdf(runid='120982W31', yname='xki', legend=r'$\chi_{\mathrm{i, total}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120982W31', yname='condipr', legend=r'$\chi_{\mathrm{i, condipr}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120982W31', yname='xkimmm', legend=r'$\chi_{\mathrm{i, mmm}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120982W31', yname='condiwnc', legend=r'$\chi_{\mathrm{i, neo}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120982W31', yname='xkiw20', legend=r'$\chi_{\mathrm{i, w20}}$', zval=0.559, source='cdf'),
        # # PlotDataCdf(runid='120968W34', yname='xkiw20', legend=r'$\chi_{\mathrm{i, w20}}$', zval=0.559, source='cdf'),
        # # PlotDataCdf(runid='120968W34', yname='condewnc', xname='rho', zval=0.559, source='cdf'),
        # title_override='120982', xmin=0.02,

        # PlotDataCdf(runid='120982W31', yname='xki', legend=r'$\chi_{\mathrm{i, total}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120982W31', yname='condipr', legend=r'$\chi_{\mathrm{i, condipr}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120982W31', yname='xkimmm', legend=r'$\chi_{\mathrm{i, mmm}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120982W31', yname='condiwnc', legend=r'$\chi_{\mathrm{i, neo}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120982W31', yname='fkch2', legend=r'$\chi_{\mathrm{i, fkch2}}$', zval=0.559, source='cdf'),
        # PlotDataCdf(runid='120982W31', yname='fkchh', legend=r'$\chi_{\mathrm{i, fkchh}}$', zval=0.559, source='cdf'),
        # PlotDataCdf(runid='120982W31', yname='fkchz', legend=r'$\chi_{\mathrm{i, fkchz}}$', zval=0.559, source='cdf'),
        # # PlotDataCdf(runid='120982W31', yname='xkiw20', legend=r'$\chi_{\mathrm{i, w20}}$', zval=0.559, source='cdf'),
        # # PlotDataCdf(runid='120968W34', yname='xkiw20', legend=r'$\chi_{\mathrm{i, w20}}$', zval=0.559, source='cdf'),
        # # PlotDataCdf(runid='120968W34', yname='condewnc', xname='rho', zval=0.559, source='cdf'),
        # title_override='120982', xmin=0.02,

        # PlotDataCdf(runid='18476T02', yname='xki', legend=r'$\chi_{\mathrm{i, total}}$', zval=7, source='mmm'),
        # PlotDataCdf(runid='18476T02', yname='condipr', legend=r'$\chi_{\mathrm{i, condipr}}$', zval=7, source='cdf'),
        # PlotDataCdf(runid='18476T02', yname='xkimmm', legend=r'$\chi_{\mathrm{i, mmm}}$', zval=7, source='cdf'),
        # PlotDataCdf(runid='18476T02', yname='condiwnc', legend=r'$\chi_{\mathrm{i, neo}}$', zval=7, source='cdf'),
        # # PlotDataCdf(runid='18476T02', yname='fkch2', legend=r'$\chi_{\mathrm{i, fkch2}}$', zval=0.559, source='cdf'),
        # # PlotDataCdf(runid='18476T02', yname='fkchh', legend=r'$\chi_{\mathrm{i, fkchh}}$', zval=0.559, source='cdf'),
        # # PlotDataCdf(runid='18476T02', yname='fkchz', legend=r'$\chi_{\mathrm{i, fkchz}}$', zval=0.559, source='cdf'),
        # # PlotDataCdf(runid='120982W31', yname='xkiw20', legend=r'$\chi_{\mathrm{i, w20}}$', zval=0.559, source='cdf'),
        # # PlotDataCdf(runid='120968W34', yname='xkiw20', legend=r'$\chi_{\mathrm{i, w20}}$', zval=0.559, source='cdf'),
        # # PlotDataCdf(runid='120968W34', yname='condewnc', xname='rho', zval=0.559, source='cdf'),
        # title_override='18476', xmin=0.02, ymax=6, xmax=0.8,

        # PlotDataCdf(runid='120968W34', yname='xki', legend=r'$\chi_{\mathrm{i, total}}$', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120968W34', yname='condipr', legend=r'$\chi_{\mathrm{i, condipr}}$', zval=0.559, source='cdf'),
        # PlotDataCdf(runid='120968W34', yname='xkimmm', legend=r'$\chi_{\mathrm{i, mmm}}$', zval=0.559, source='cdf'),
        # PlotDataCdf(runid='120968W34', yname='condiwnc', legend=r'$\chi_{\mathrm{i, neo}}$', zval=0.559, source='cdf'),
        # title_override='120968', xmin=0.02, xmax=0.8,

        # PlotDataCdf(runid='129016W16', yname='xkemmm',  zval=0.533, source='cdf'),
        # PlotDataCdf(runid='129016W16', yname='xkew20',  zval=0.533, source='cdf'),
        # PlotDataCdf(runid='129016W16', yname='xkeetgm', zval=0.533, source='mmm'),
        # PlotDataCdf(runid='129016W16', yname='xkemtm',  zval=0.533, source='cdf'),
        # PlotDataCdf(runid='129016W16', yname='xkedrbm',  zval=0.533, source='cdf'),
        # title_override='129016', xmin=0.02, xmax=0.8,

        # PlotDataCdf(runid='129016W12', yname='mmmtime', xname='time', legend='Default', source='mmm'),
        # PlotDataCdf(runid='129016W21', yname='mmmtime', xname='time', legend=r'$g \geq 1$', source='mmm'),
        # PlotDataCdf(runid='129016W22', yname='mmmtime', xname='time', legend=r'$g \geq 1$, 0 PT', source='mmm'),
        # title_override='129016',

        # PlotDataCdf(runid='129016W12', yname='mmmtime', xname='time', legend=r'9.0.7', source='mmm'),
        # PlotDataCdf(runid='129016W24', yname='mmmtime', xname='time', legend=r'9.0.9', source='mmm'),
        # PlotDataCdf(runid='129016W23', yname='mmmtime', xname='time', legend=r'9.0.9 $g \geq 1$', source='mmm'),
        # title_override='129016', 

        # PlotDataCdf(runid='120982A09', yname='ne', xname='rho', source='cdf'),
        # PlotDataCdf(runid='120982A09', yname='ne', xname='rho', source='mmm'),
        # PlotDataCdf(runid='120982A09', yname='ni', xname='rho', source='mmm'),
        # title_override='129016', 

        # PlotDataCdf(runid='120982A09', yname='nf', xname='rho', source='cdf'),
        # # PlotDataCdf(runid='120982A09', yname='nh', xname='rho', source='mmm'),
        # PlotDataCdf(runid='120982A09', yname='nz', xname='rho', source='mmm'),
        # title_override='120982A09', 

        # PlotDataCdf(runid='129016W12', yname='condepr', xname='rho', zval=0.5, source='cdf'),
        # PlotDataCdf(runid='129016W12', yname='nh0', xname='rho',   zval=0.5,  source='cdf'),
        # PlotDataCdf(runid='129016W12', yname='nd', xname='rho',   zval=0.5,  source='cdf'),
        # PlotDataCdf(runid='129016W12', yname='nt', xname='rho',   zval=0.5,  source='cdf'),
        # PlotDataCdf(runid='129016W12', yname='nhe3', xname='rho',   zval=0.5,  source='cdf'),
        # PlotDataCdf(runid='129016W12', yname='nhe4', xname='rho',   zval=0.5,  source='cdf'),
        # title_override='129016', 

        # PlotDataCdf(runid='120968A02', yname='bpol', zval=0.5, legend='cdf', source='cdf'),
        # PlotDataCdf(runid='120968A02', yname='bpol', zval=0.5, legend='mmm', source='mmm'),
        # title_override='120968A02', 

        # PlotDataCdf(runid='129016W25', yname='omega', zval=0.51, legend='9.0.7 Only Confinement', source='cdf'),
        # PlotDataCdf(runid='129016W26', yname='omega', zval=0.51, legend='9.0.7 All three regions', source='cdf'),
        # PlotDataCdf(runid='129016W13', yname='omega', zval=0.51, legend='Analysis', source='cdf'),
        # title_override='129016', 

        # PlotDataCdf(runid='147634T61', yname='gne', zval=2.53, legend='9.0.7 Only Confinement', source='mmm'),
        # PlotDataCdf(runid='129016W26', yname='omega', zval=0.51, legend='9.0.7 All three regions', source='cdf'),
        # PlotDataCdf(runid='129016W13', yname='omega', zval=0.51, legend='Analysis', source='cdf'),
        # title_override='129016', 

        # PlotDataCdf(runid='120968W34', yname='ti', legend=r'MMM', zval=0.51),
        # PlotDataCdf(runid='120968A02', yname='ti', legend=r'Analysis', zval=0.51),



        # ------------------------------------------------------------------------------- #


        # PlotDataCdf(runid='120968W34', yname='ti', xname='x', legend=r'MMM', zval=0.56, source='raw', ymult=1e-3),
        # PlotDataCdf(runid='120968A02', yname='ti', xname='x', legend=r'Analysis', zval=0.56, source='raw', ymult=1e-3),
        # title_override='120968 (W34, A02), t = 0.56$\,$s',
        # xmin=0, xmax=1, ymin=0, ymax=1, ylabel_override=r'$T_{\rm i}$ (keV)'

        # PlotDataCdf(runid='120968W34', yname='te', xname='x', legend=r'MMM',      zval=0.56, source='raw', ymult=1e-3),
        # PlotDataCdf(runid='120968A02', yname='te', xname='x', legend=r'Analysis', zval=0.56, source='raw', ymult=1e-3),
        # title_override='120968 (W34, A02), t = 0.56$\,$s',
        # xmin=0, xmax=1, ymin=0, ymax=1, ylabel_override=r'$T_{\rm e}$ (keV)'

        # PlotDataCdf(runid='120982W31', yname='ti', xname='x', legend=r'MMM',      zval=0.62, source='raw', ymult=1e-3),
        # PlotDataCdf(runid='120982A09', yname='ti', xname='x', legend=r'Analysis', zval=0.62, source='raw', ymult=1e-3),
        # title_override='120968 (W34, A02), t = 0.62$\,$s',
        # xmin=0, xmax=1, ymin=0, ymax=1.6, ylabel_override=r'$T_{\rm i}$ (keV)',
        # yticks=np.arange(0, 1.8, step=0.4)

        # PlotDataCdf(runid='120982W31', yname='te', xname='x', legend=r'MMM',      zval=0.62, source='raw', ymult=1e-3),
        # PlotDataCdf(runid='120982A09', yname='te', xname='x', legend=r'Analysis', zval=0.62, source='raw', ymult=1e-3),
        # title_override='120968 (W34, A02), t = 0.62$\,$s',
        # xmin=0, xmax=1, ymin=0, ymax=1, ylabel_override=r'$T_{\rm e}$ (keV)',
        # yticks=np.arange(0, 1.8, step=0.4)

        # PlotDataCdf(runid='138536W03', yname='ti', xname='x', legend=r'MMM', zval=0.755, source='raw', ymult=1e-3),
        # PlotDataCdf(runid='138536A01', yname='ti', xname='x', legend=r'Analysis', zval=0.755, source='raw', ymult=1e-3),
        # title_override='138536 (W03, A01), t = 0.755$\,$s',
        # xmin=0, xmax=1, ymin=0, ymax=1.2, ylabel_override=r'$T_{\rm i}$ (keV)'

        # PlotDataCdf(runid='138536W03', yname='te', xname='x', legend=r'MMM', zval=0.755, source='raw', ymult=1e-3),
        # PlotDataCdf(runid='138536A01', yname='te', xname='x', legend=r'Analysis', zval=0.755, source='raw', ymult=1e-3),
        # title_override='138536 (W03, A01), t = 0.755$\,$s',
        # xmin=0, xmax=1, ymin=0, ymax=1.2, ylabel_override=r'$T_{\rm e}$ (keV)'
  
        # PlotDataCdf(runid='138536W03', yname='condipr', xname='x', legend=r'$\chi_{\rm i, tot}$', zval=0.755, source='raw', ymult=1e-4),
        # PlotDataCdf(runid='138536W03', yname='xkimmm', xname='x',  legend=r'$\chi_{\rm i, MMM}$', zval=0.755, source='raw', ymult=1e-4),
        # PlotDataCdf(runid='138536W03', yname='xkemmm', xname='x',  legend=r'$\chi_{\rm e, MMM}$', zval=0.755, source='raw', ymult=1e-4),
        # title_override='138536 (W03), t = 0.755$\,$s',
        # xmin=0, xmax=1, ymin=0, ylabel_override=r'(m$^2$/s)'

        # PlotDataCdf(runid='120968W34', yname='condipr', xname='x', legend=r'$\chi_{\rm i, tot}$', zval=0.51, source='raw', ymult=1e-4),
        # PlotDataCdf(runid='120968W34', yname='xkimmm', xname='x', legend=r'$\chi_{\rm i, MMM}$', zval=0.51, source='raw', ymult=1e-4),
        # PlotDataCdf(runid='120968W34', yname='xkemmm', xname='x', legend=r'$\chi_{\rm e, MMM}$', zval=0.51, source='raw', ymult=1e-4),
        # title_override=r'120968 (W34), t = 0.51$\,$s',
        # xmin=0, xmax=1, ymin=0, ylabel_override=r'(m$^2$/s)'


        ## FIG 10

        # PlotDataCdf(runid='120982W31', yname='ti', xname='x', legend=r'MMM',      zval=0.62, source='raw', ymult=1e-3),
        # PlotDataCdf(runid='120982A09', yname='ti', xname='x', legend=r'Analysis', zval=0.62, source='raw', ymult=1e-3),
        # title_override=r'120982 (W31, A09), t = 0.62$\,$s',
        # xmin=0, xmax=1, ymin=0, ymax=1.6, ylabel_override=r'$T_{\rm i}$ (keV)',
        # yticks=np.arange(0, 1.8, step=0.4)

        # PlotDataCdf(runid='120982W31', yname='te', xname='x', legend=r'MMM',      zval=0.755, source='raw', ymult=1e-3),
        # PlotDataCdf(runid='120982A09', yname='te', xname='x', legend=r'Analysis', zval=0.755, source='raw', ymult=1e-3),
        # title_override=r'120982 (W31, A09), t = 0.62$\,$s',
        # xmin=0, xmax=1, ymin=0, ymax=1.6, ylabel_override=r'$T_{\rm e}$ (keV)',
        # yticks=np.arange(0, 1.8, step=0.4)
  
        # PlotDataCdf(runid='120982W31', yname='condipr', xname='x', legend=r'$\chi_{\rm i, tot}$', zval=0.62, source='raw', ymult=1e-4),
        # PlotDataCdf(runid='120982W31', yname='xkimmm',  xname='x', legend=r'$\chi_{\rm i, MMM}$', zval=0.62, source='raw', ymult=1e-4),
        # PlotDataCdf(runid='120982W31', yname='xkemmm',  xname='x', legend=r'$\chi_{\rm e, MMM}$', zval=0.62, source='raw', ymult=1e-4),
        # title_override=r'120982 (W31), t = 0.62$\,$s',
        # xmin=0, xmax=1, ymin=0, ylabel_override=r'(m$^2$/s)'
  
        ##### ETGM OPTIMIZATION PLOTS

        ## 1: Default plots, no new optimization - show convergence inconsistencies wrt PT Solver convergence
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20003, legend=r'10000',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20002, legend=r'1000',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20000, legend=r'50',),

        # # 2: 10x5 comparable to 10000
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20003, legend=r'10000',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20002, legend=r'1000',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20071, legend=r'10x5',),

        # # 3: 2% tolerance + gamma cutoff = no change in results
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20071, legend=r'(5000) 10x5',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20089, legend=r'(4050) 10x5 2%',),

        ## 4: 5% tolerance even faster, but with some error 
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20071, legend=r'(5000) 10x5',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20089, legend=r'(4050) 10x5 2%',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20086, legend=r'(3300) 10x5 5%',),

        ## 5: 5% tolerance on par with 1000 default scans
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20089, legend=r'(4050) 10x5 2%',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20086, legend=r'(3300) 10x5 5%',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20002, legend=r'(1E5)  1000 Default',),

        ## 6: New discharge, showing 10x5 = 10000 default
        # PlotDataCsv(runid='129017A04', yname='xteETGM', xname='rho', scan_num=20022, legend=r'(1E7)  10000 Default',),
        # PlotDataCsv(runid='129017A04', yname='xteETGM', xname='rho', scan_num=20000, legend=r'(5000) 10x5',),

        ## 7: Showing different tolerances
        # PlotDataCsv(runid='129017A04', yname='xteETGM', xname='rho', scan_num=20000, legend=r'(5000) 10x5',),
        # PlotDataCsv(runid='129017A04', yname='xteETGM', xname='rho', scan_num=20019, legend=r'(2550) 10x5 2%',),
        # PlotDataCsv(runid='129017A04', yname='xteETGM', xname='rho', scan_num=20020, legend=r'(2120) 10x5 5%',),

        # 8: Different tolerance compared to 50 default
        # PlotDataCsv(runid='129017A04', yname='xteETGM', xname='rho', scan_num=20019, legend=r'(2550) 10x5 2%',),
        # PlotDataCsv(runid='129017A04', yname='xteETGM', xname='rho', scan_num=20020, legend=r'(2120) 10x5 5%',),
        # PlotDataCsv(runid='129017A04', yname='xteETGM', xname='rho', scan_num=20029, legend=r'(1967) 7x10 2%',),
        # PlotDataCsv(runid='129017A04', yname='xteETGM', xname='rho', scan_num=25004, legend=r'(1932) 7x10 2%',),
        # PlotDataCsv(runid='129017A04', yname='xteETGM', xname='rho', scan_num=20021, legend=r'(5000) 50 Default',),

        # PlotDataCsv(runid='129017A04', yname='xteETGM', xname='rho', scan_num=20028, legend=r'(2560) delta kyrhos',),

        ## 9: kyrhos flat-line with gamma cutoff
        # PlotDataCsv(runid='121123K55', yname='gmaETGM', xname='rho', scan_num=20071, legend=r'(5000) 10x5',),
        # PlotDataCsv(runid='121123K55', yname='gmaETGM', xname='rho', scan_num=20089, legend=r'(4050) 10x5 2%',),
        # PlotDataCsv(runid='121123K55', yname='gmaETGM', xname='rho', scan_num=20086, legend=r'(3300) 10x5 5%',),

        # 10: Multi-tiered tolerance check (using kyrhos value)
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20089, legend=r'(4050) 10x5 2%',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20086, legend=r'(3300) 10x5 5%',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20099, legend=r'(3610) 10x5 2%, 5%',),

        # # 11: Different kyrhos ranges
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20109, legend=r'(4050) kyrhos 1-50',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20108, legend=r'(4010) kyrhos 1-100',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20110, legend=r'(3960) kyrhos 1-25',),

        ## 12: Exp vs linear (no gc)
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20003, legend=r'10000',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20112, legend=r'(4860) 10x5 2% no gc Linear',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20113, legend=r'(4640) 10x5 2% no gc Exp',),

        ## 13: Exp vs linear (w/ gc)
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20003, legend=r'10000',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20114, legend=r'(4100) 10x5 2% Linear',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20089, legend=r'(4050) 10x5 2% Exp',),

        ## 14: Exp vs linear (w/ gc, 5% tol)
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20003, legend=r'10000',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20115, legend=r'(3620) 10x5 5% Linear',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20086, legend=r'(3300) 10x5 5% Exp',),

        ## Output inconsistency
        # PlotDataCsv(runid='129016A04', yname='gmaETGM', xname='rho', scan_num=25009, legend=r'10000x1',),
        # PlotDataCsv(runid='129016A04', yname='gmaETGM', xname='rho', scan_num=25005, legend=r'7x7',),
        # PlotDataCsv(runid='129016A04', yname='gmaETGM', xname='rho', scan_num=25012, legend=r'8x10',),

        ## Inconsistency gma vs kyrhos
        # PlotDataCsv(runid='129016A04', yname='gmaETGM', xname='kyrhosETGM', scan_num=25027, rho_value=0.13),
        # PlotDataCsv(runid='129016A04', yname='gmaETGM', xname='kyrhosETGM', scan_num=25027, rho_value=0.14),

        ## Testing
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20089, legend=r'(4050) 10x5 2%',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20086, legend=r'(3300) 10x5 5%',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20111, legend=r'(4100) 10x5 2% Linear',),

        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20092, legend=r'(3810) $<$1E-2',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20096, legend=r'(3970) Two tolerances 5e5',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20099, legend=r'(3610) kyrhos $>$ 5',),

        # PlotDataCsv(runid='121123K55', yname='kyrhosETGM', xname='rho', scan_num=20098, legend=r'(3760) Two tolerances 1e6',),

        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20089, legend=r'(4050) 10x5 2%',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20086, legend=r'(3300) 10x5 5%',),
        # # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20123, legend=r'(2802) 6x10 5%',),
        # # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20124, legend=r'(3162) 6x10 4%',),
        # # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20127, legend=r'(3345) 5x10 2%',),
        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20132, legend=r'(3325) 7x10 2%',),
        # # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20130, legend=r'(2955) 5x10 test',),

        # PlotDataCsv(runid='121123K55', yname='xteETGM', xname='rho', scan_num=20086, legend=r'(3300) 10x5 5%',),

        # ### DBM Testing
        # PlotDataCsv(runid='120968A02', yname='xtiDBM', xname='rho', scan_num=25031, legend=r'5000x1',),
        # PlotDataCsv(runid='120968A02', yname='xtiDBM', xname='rho', scan_num=25029, legend=r'1000x1',),
        # PlotDataCsv(runid='120968A02', yname='xtiDBM', xname='rho', scan_num=25028, legend=r'(2000) 20x1',),
        # PlotDataCsv(runid='120968A02', yname='xtiDBM', xname='rho', scan_num=25032, legend=r'(2107) 7x10',),
        # PlotDataCsv(runid='120968A02', yname='xtiDBM', xname='rho', scan_num=25033, legend=r'(2107) 7x10 gma-wexb',),

    #     PlotDataCsv(runid='120968A02', yname='xtiDBM', xname='rho', scan_num=25038, legend=r'(1526) 7x10 1e3-wexb',),
    #     PlotDataCsv(runid='120968A02', yname='xtiDBM', xname='rho', scan_num=25039, legend=r'(1435) 7x10 1e4-wexb',),
    #     PlotDataCsv(runid='120968A02', yname='xtiDBM', xname='rho', scan_num=25037, legend=r'(2000) 20x1 gma-wexb',),

        # PlotDataCsv(runid='138536A01', yname='gmaDBM', xname='kyrhosDBM', scan_num=25059, rho_value=0.85,),
        # PlotDataCsv(runid='120968A02', yname='gmaDBM', xname='rho', scan_num=26026),
        # PlotDataCsv(runid='120968A02', yname='gmaDBM', xname='rho', scan_num=26027),
        # PlotDataCsv(runid='120968A02', yname='gmaDBM', xname='rmina', scan_num=26032, legend='wexb = 0'),
        # xmax=0.8,

        # PlotDataCsv(runid='120968A02', yname='nR8TOMSQZ', xname='rho', scan_num=26044),
        # PlotDataCsv(runid='15334T03', yname='xteETGM', xname='rho', scan_num=26210),
        # PlotDataCsv(runid='15334T03', yname='xteETGM', xname='rho', scan_num=26217),

        # PlotDataCsv(runid='80200A13', yname='xtiW20', xname='rho', scan_num=25),
        # PlotDataCsv(runid='80200A13', yname='xtiW20', xname='rho', scan_num=17),
        
        # PlotDataCsv(runid='120968A02', yname='ne', xname='rho',  scan_num=41),
        # PlotDataCsv(runid='120968A02', yname='zni', xname='rho', scan_num=41),
        # PlotDataCsv(runid='120968A02', yname='ni', xname='rho',  scan_num=41),
        
        # PlotDataCsv(runid='120968A02', yname='zave', xname='rho',  scan_num=42),
        # PlotDataCsv(runid='120968A02', yname='zeff', xname='rho',  scan_num=42),

        # PlotDataCsv(runid='120968A02', yname='gmaDBM', legend='akima', scan_num=21),
        # PlotDataCsv(runid='120968A02', yname='gmaDBM', legend='interp (cubic)', scan_num=24),
        # PlotDataCsv(runid='120968A02', yname='gmaDBM', legend='ptsolver', scan_num=23),
        # ymax=2e5,

        # PlotDataCsv(runid='120968A02', yname='gmaDBM', legend='akima', scan_num=28),
        # PlotDataCsv(runid='120968A02', yname='gmaDBM', legend='interp (cubic)', scan_num=29),
        # PlotDataCsv(runid='120968A02', yname='gmaDBM', legend='ptsolver', scan_num=30),
        # ymax=2e5,

        # PlotDataCdf(runid='129016A04', yname='tf', zval=0.629),
        # PlotDataCdf(runid='129016A04', yname='tfast', zval=0.629),
        # PlotDataCdf(runid='129016A04', yname='tmhdf', zval=0.629),

        # PlotDataCdf(runid='120968A02', yname='tf', zval=0.629),
        # PlotDataCdf(runid='120968A02', yname='tfpa', zval=0.629),
        # PlotDataCdf(runid='120968A02', yname='tfpp', zval=0.629),

        # PlotDataCdf(runid='129016A04', yname='ufastpa', zval=0.629),
        # PlotDataCdf(runid='129016A04', yname='ufastpp', zval=0.629),
        
        # PlotDataCdf(runid='129016A04', yname='gxi', zval=0.629, source='cdf'),
        # PlotDataCdf(runid='129016A04', yname='gxi', zval=0.629),
        # # ymax=2e5,

        # MTM Scan Counts
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E 10000', scan_num=129),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E 100000', scan_num=136),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E 1000000', scan_num=139),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='L 100000', scan_num=135),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='L 500000', scan_num=137),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='L 1000000', scan_num=138),

        # PlotDataCsv(runid='138536A01', yname='xte', legend='Converged', scan_num=138),  # 1M Linear
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E 200', scan_num=126),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='L 4000',  scan_num=140),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E+ 200',  scan_num=141),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E+ 100',  scan_num=144),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E 2000', scan_num=127),
        # # PlotDataCsv(runid='138536A01', yname='xte', legend='E 5000', scan_num=128),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E 10000', scan_num=129),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='L 200',   scan_num=130),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='L 2000',  scan_num=131),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='L 5000',  scan_num=132),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='L 10000', scan_num=133),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='L 100000', scan_num=135),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E 100000', scan_num=136),

        # PlotDataCsv(runid='138536A01', yname='xte', legend='Converged', scan_num=138),  # 1M Linear
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E+ 100',  scan_num=144),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E 100',  scan_num=145),
        
        # PlotDataCsv(runid='138536A01', yname='xte', legend='Converged', scan_num=138),  # 1M Linear
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E 200 ++',  scan_num=147),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E 200', scan_num=126),
        
        # PlotDataCsv(runid='138536A01', yname='xte', legend='Converged', scan_num=138),  # 1M Linear
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E 200 ++',  scan_num=147),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='L 4000',  scan_num=140),

        # PlotDataCsv(runid='138536A01', yname='xte', legend='Converged', scan_factor=0.252, scan_num=10000),  # 1M Linear
        # PlotDataCsv(runid='138536A01', yname='xte', legend='E 2E2',     scan_factor=0.252, scan_num=10001),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='L 1E4',     scan_factor=0.252, scan_num=10006),

        # PlotDataCsv(runid='138536A01', yname='xte', legend='Converged',            scan_num=10000),  # 1M Linear
        # PlotDataCsv(runid='138536A01', yname='xte', legend=r'$n_{\rm E} = 200$',   scan_num=10001),
        # PlotDataCsv(runid='138536A01', yname='xte', legend=r'$n_{\rm L} = 4000$', scan_num=10004),
        # PlotDataCsv(runid='138536A01', yname='xte', legend=r'$n_{\rm L} = 10000$', scan_num=10006),

        # PlotDataCsv(runid='141716A80', yname='fte', legend='Converged',          scan_num=10000),  # 1E6 E
        # PlotDataCsv(runid='141716A80', yname='fte', legend=r'$n_{\rm E} = 50$',  scan_num=10001),
        # PlotDataCsv(runid='141716A80', yname='fte', legend=r'$n_{\rm E} = 100$', scan_num=10002),
        # PlotDataCsv(runid='141716A80', yname='fte', legend=r'$n_{\rm E} = 200$', scan_num=10003),

        # PlotDataCsv(runid='141716A80', yname='fte', legend='Converged',           scan_num=10000),  # 1E6 E
        # PlotDataCsv(runid='141716A80', yname='fte', legend=r'$n_{\rm L} = 2000$', scan_num=10004),
        # PlotDataCsv(runid='141716A80', yname='fte', legend=r'$n_{\rm L} = 4000$', scan_num=10005),
        # PlotDataCsv(runid='141716A80', yname='fte', legend=r'$n_{\rm L} = 8000$', scan_num=10006),

        # PlotDataCsv(runid='141716A80', yname='fte', legend='Converged',           scan_factor=0.4, scan_num=10000),  # 1E6 E
        # PlotDataCsv(runid='141716A80', yname='fte', legend=r'$n_{\rm L} = 8000$', scan_factor=0.4, scan_num=10006),
        # PlotDataCsv(runid='141716A80', yname='fte', legend=r'$n_{\rm E} = 200$',  scan_factor=0.4, scan_num=10003),
        # xmin=0.36, xmax=0.82, ymin=0.6, ymax=1.6,

        # PlotDataCsv(runid='176523L01', yname='fte', legend='Converged',            scan_factor=0.6, scan_num=10000),  # 1E6 E
        # PlotDataCsv(runid='176523L01', yname='fte', legend=r'$n_{\rm L} = 16000$', scan_factor=0.6, scan_num=10010),
        # PlotDataCsv(runid='176523L01', yname='fte', legend=r'$n_{\rm E} = 200$',   scan_factor=0.6, scan_num=10003),
        # xmin=0.3, xmax=0.85, ymin=0.15, ymax=0.4,
        
        # PlotDataCsv(runid='138536A01', yname='xte', legend='Converged', scan_num=138),  # 1M Linear
        # PlotDataCsv(runid='138536A01', yname='xte', legend='L 500',   scan_num=149),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='L 2000',  scan_num=131),
        # PlotDataCsv(runid='138536A01', yname='xte', legend='L 4000',  scan_num=140),   
        # # ymax=2e5,


        # PlotDataCsv(runid='176523L01', yname='fte', legend='Converged',           scan_factor=0.6, scan_num=10000),  # 1E6 E
        # # PlotDataCsv(runid='176523L01', yname='fte', legend=r'$n_{\rm L} = 16000$', scan_factor=0.6, scan_num=10010),
        # PlotDataCsv(runid='176523L01', yname='fte', legend=r'$n_{\rm E} = 200$',  scan_factor=0.6, scan_num=10003),
        # # PlotDataCsv(runid='176523L01', yname='fte', legend=r'$n_{\rm \mathcal{E}} = 149$',  scan_factor=0.6, scan_num=10023),
        # PlotDataCsv(runid='176523L01', yname='fte', legend=r'$\overline{n}_{\rm \mathcal{E}} = 245$',  scan_factor=0.6, scan_num=10026),
        # xmin=0.58, xmax=0.78, ymin=0.35, ymax=0.392,


        # PlotDataCsv(runid='141716A80', yname='fte', legend='Converged',           scan_factor=0.5, scan_num=10000),  # 1E6 E
        # PlotDataCsv(runid='141716A80', yname='fte', legend=r'$n_{\rm L} = 8000$', scan_factor=0.5, scan_num=10006),
        # PlotDataCsv(runid='141716A80', yname='fte', legend=r'$n_{\rm E} = 200$',  scan_factor=0.5, scan_num=10003),
        # PlotDataCsv(runid='141716A80', yname='fte', legend=r'$\overline{n}_{\rm \mathcal{E}} = 246$',  scan_factor=0.5, scan_num=10027),
        # xmin=0.45, xmax=0.575, ymin=2.4, ymax=3.2,


        # PlotDataCsv(runid='141716A80', yname='fte', legend='Converged',           scan_factor=0.4, scan_num=10000),  # 1E6 E
        # PlotDataCsv(runid='141716A80', yname='fte', legend=r'$n_{\rm L} = 8000$', scan_factor=0.4, scan_num=10006),
        # PlotDataCsv(runid='141716A80', yname='fte', legend=r'$n_{\rm E} = 200$',  scan_factor=0.4, scan_num=10003),
        # PlotDataCsv(runid='141716A80', yname='fte', legend=r'$\overline{n}_{\rm \mathcal{E}} = 246$',  scan_factor=0.4, scan_num=10027),
        # xmin=0.58, xmax=0.8, ymin=1.2, ymax=1.52,


        # PlotDataCsv(runid='120968A02', yname='gxi', legend='Akima',     scan_num=41),  # 1E6 E
        # PlotDataCsv(runid='120968A02', yname='gxi', legend='PT Solver', scan_num=42),
        
        # PlotDataCsv(runid='141716A80', yname='xteMTM', scan_num=20, scan_factor=1.5),
        # PlotDataCsv(runid='141716A80', yname='xteMTM', scan_num=19, scan_factor=1.5),
        # PlotDataCsv(runid='141716A80', yname='xteMTM', scan_num=21, scan_factor=1.5),

        # PlotDataCsv(runid='120968A02', xname='rmina', yname='xne', legend='High', scan_num=45),
        # PlotDataCsv(runid='120982A09', xname='rmina', yname='xne', legend='Low', scan_num=9),
        # xmin=0.0, xmax=0.8,

        # wexb = 0
        PlotDataCsv(runid='120968A02', xname='rmina', yname='nuste', legend='High', scan_num=50),
        PlotDataCsv(runid='120982A09', xname='rmina', yname='nuste', legend='Low', scan_num=12),
        xmin=0.0, xmax=0.8, ymin=0, logy=True, logx=True,# ymax=1.6e6,

        # PlotDataCdf(runid='120968A02', xname='rmina', yname='nuei', legend='High', zval=0.56),
        # PlotDataCdf(runid='120982A09', xname='rmina', yname='nuei', legend='Low',zval=0.62),
        # xmin=0.0, xmax=0.8, ymin=0, ymax=1.6e6,
    )

    main(fig_data)
