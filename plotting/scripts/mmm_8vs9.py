# Standard Packages
import sys; sys.path.insert(0, '../'); sys.path.insert(0, '../../')
import logging

# Third Party Packages
from matplotlib.pyplot import rcParams

# Local Packages
import modules.utils as utils
from plotting.modules.plotstyles import PlotStyles, StyleType
from plotting.plot_variables import FigData, PlotDataCdf, PlotDataCsv, main


_log = logging.getLogger(__name__)


if __name__ == '__main__':
    """Run this module directly to plot variable data stored in CDF or CSV files"""

    utils.init_logging()

    # Initialize visual styles for the plot
    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.AIP2,
    )

    # Define settings for the plot
    fig_data = FigData(
        replace_offset_text=False,
        allow_title_runid=1,
        allow_title_time=1,
        allow_title_factor=True,
        allow_title_rho=True,
        invert_y_axis=False,
        invert_x_axis=False,
        nomralize_y_axis=False,
        nomralize_x_axis=False,
        savename_append='',
        title_override=' ',
        ylabel_override='',
        xlabel_override='',
        # xmax=1,
    )

    rcParams.update({
        'legend.fontsize': 8,
    })

    savefig = 1
    savedata = 0

    fig_data.set(

        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------
        ## 129016: 0.3 to 0.5 (some to 0.4)
        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------

        ## 32000 vs 8000 PTCLS
        # PlotDataCdf(runid='129016Z33', yname='ti', xname='rho', zval=0.500, legend='NPTCLS=32000'),
        # PlotDataCdf(runid='129016Z35', yname='ti', xname='rho', zval=0.500, legend='NPTCLS=8000'),
        ## ---------------------------------------------------------------------------


        ## TESTING PPHI
        # PlotDataCdf(runid='129016Z29', yname='xkemtm', xname='rho', zval=0.32, legend=r'$\mathtt{lpredict\_pphi=0}$'),
        # PlotDataCdf(runid='129016Z24', yname='xkemtm', xname='rho', zval=0.32, legend=r'$\mathtt{lpredict\_pphi=1}$'),

        # PlotDataCdf(runid='129016Z42', yname='xkemtm', xname='rho', zval=0.32, legend=r'$\mathtt{lpredict\_pphi=0}$'),
        # PlotDataCdf(runid='129016Z41', yname='xkemtm', xname='rho', zval=0.32, legend=r'$\mathtt{lpredict\_pphi=1}$'),
        
        
        # PlotDataCdf(runid='129016Z31', yname='xkemtm', xname='rho', zval=0.32, legend=r''),
        # PlotDataCdf(runid='129016Z30', yname='xkemtm', xname='rho', zval=0.32, legend=r'All Regions + pphi enabled'),
        # PlotDataCdf(runid='129016Z33', yname='te', xname='rho', zval=0.32, legend=r'v9 Confinement Only'), 

        ## MMM 8 vs 9
        # PlotDataCdf(runid='129016Z29', yname='xkiw20', xname='rho', zval=0.301, legend='MMM 8.2.1'),
        # PlotDataCdf(runid='129016Z33', yname='xkiw20', xname='rho', zval=0.301, legend='MMM 9.0.7'),

        # PlotDataCdf(runid='129016Z29', yname='xkew20', xname='rho', zval=0.301, legend='MMM 8.2.1'),
        # PlotDataCdf(runid='129016Z33', yname='xkew20', xname='rho', zval=0.301, legend='MMM 9.0.7'),

        # PlotDataCdf(runid='129016Z29', yname='xkemtm', xname='rho', zval=0.301, legend='MMM 8.2.1'),
        # PlotDataCdf(runid='129016Z33', yname='xkemtm', xname='rho', zval=0.301, legend='MMM 9.0.7'),

        # PlotDataCdf(runid='129016Z29', yname='xkemmm', xname='rho', zval=0.33, legend='MMM 8.2.1'),
        # PlotDataCdf(runid='129016Z33', yname='xkemmm', xname='rho', zval=0.33, legend='MMM 9.0.7'),

        # PlotDataCdf(runid='129016Z29', yname='ti', xname='rho', zval=0.301, legend='MMM 8.2.1'),
        # PlotDataCdf(runid='129016Z33', yname='ti', xname='rho', zval=0.301, legend='MMM 9.0.7'),

        # PlotDataCdf(runid='129016Z29', yname='te', xname='rho', zval=0.33, legend='MMM 8.2.1'),
        # PlotDataCdf(runid='129016Z33', yname='te', xname='rho', zval=0.33, legend='MMM 9.0.7'),

        # PlotDataCdf(runid='129016Z29', yname='te', xname='time', zval=0.02, timeplot=True, legend='MMM 8.2.1'),
        # PlotDataCdf(runid='129016Z33', yname='te', xname='time', zval=0.02, timeplot=True, legend='MMM 9.0.7'),
        # xmax=0.37, 

        # PlotDataCdf(runid='129016Z29', yname='xkemmm', xname='time', zval=0.02, timeplot=True, legend='MMM 8.2.1'),
        # PlotDataCdf(runid='129016Z33', yname='xkemmm', xname='time', zval=0.02, timeplot=True, legend='MMM 9.0.7'),
        # xmax=0.37, 

        # title_override='129016',
        # PlotDataCdf(runid='129016A03', yname='ti', xname='rho', zval=0.35, legend='Analysis'),
        # PlotDataCdf(runid='129016Z21', yname='ti', xname='rho', zval=0.35, legend='MMM Disabled'),
        # title_override='129016',
        ## ---------------------------------------------------------------------------

        ## MMM 8 vs 9 + ETGM
        # PlotDataCdf(runid='129016Z36', yname='ti', xname='rho', zval=0.35, legend='9.0.7 +ETGM'),
        # # PlotDataCdf(runid='129016Z43', yname='ti', xname='rho', zval=0.35, legend='9.0.7 +Horton'),
        # # PlotDataCdf(runid='129016Z33', yname='ti', xname='rho', zval=0.35, legend='9.0.7'),
        # PlotDataCdf(runid='129016Z29', yname='ti', xname='rho', zval=0.40, legend='8.2.1'),
        # PlotDataCdf(runid='129016Z21', yname='ti', xname='rho', zval=0.46, legend='Only NC'),
        # PlotDataCdf(runid='129016A03', yname='ti', xname='rho', zval=0.35, legend='Analysis'),
        # xmax=0.8, title_override='129016 (0.46s)'


        # PlotDataCdf(runid='129016Z29', yname='te', xname='time', zval=0.3, legend='MMM 8.2.1'),
        # PlotDataCdf(runid='129016Z36', yname='te', xname='time', zval=0.3, legend='MMM 9.0.7'),
        # PlotDataCdf(runid='129016Z21', yname='te', xname='time', zval=0.3, legend='Only NC'),
        # PlotDataCdf(runid='129016A03', yname='te', xname='time', zval=0.3, legend='Analysis'),
        # title_override=r'129016, $\hat{\rho} = 0.3$',
        # ## -----------------------------------------------

        # PlotDataCdf(runid='129016W03', yname='te', xname='time', zval=0.3, legend='MMM 9.0.7'),
        # PlotDataCdf(runid='129016W01', yname='te', xname='time', zval=0.3, legend='MMM 8.2.1'),
        # PlotDataCdf(runid='129016Z21', yname='te', xname='time', zval=0.3, legend='Only NC'),
        # PlotDataCdf(runid='129016A03', yname='te', xname='time', zval=0.3, legend='Analysis'),
        # title_override=r'129016, $\hat{\rho} = 0.3$',
        ## ---------------------------------------------------------------------------

        ## WALLTIME
        # PlotDataCdf(runid='129016Z29', yname='walltime', xname='time', timeplot=True, zval=0.37, legend='WALLTIME'),
        # PlotDataCdf(runid='129016Z29', yname='cpmcfi', xname='time', timeplot=True, zval=0.37, legend='CPMCFI'),
        # title_override='MMM v8.2.1', xmax=0.37, ymax=1.75,

        # PlotDataCdf(runid='129016Z33', yname='walltime', xname='time', timeplot=True, zval=0.37, legend='WALLTIME'),
        # PlotDataCdf(runid='129016Z33', yname='cpmcfi', xname='time', timeplot=True, zval=0.37, legend='CPMCFI'),
        # title_override='MMM v9.0.7', xmax=0.37, ymax=1.75,

        # PlotDataCdf(runid='129016Z35', yname='walltime', xname='time', timeplot=True, zval=0.37, legend='WALLTIME'),
        # PlotDataCdf(runid='129016Z35', yname='cpmcfi', xname='time', timeplot=True, zval=0.37, legend='CPMCFI'),
        # title_override='MMM v9.0.7 (8000 PTCLS)',
        ## ---------------------------------------------------------------------------


        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------
        ## 120968: 0.1 to 0.6
        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------
        ## MMM 8 vs 9 + ETGM
        # PlotDataCdf(runid='120968A02', yname='te', xname='rho', zval=0.511, legend='Analysis'),
        # PlotDataCdf(runid='120968W34', yname='te', xname='rho', zval=0.511, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='120968W35', yname='te', xname='rho', zval=0.511, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='120968W33', yname='te', xname='rho', zval=0.511, legend='9.0.7'),
        # xmax=1, xmin=0.02, title_override='120968',

        # PlotDataCdf(runid='120968W34', yname='xkew20',  xname='xb', zval=0.559, legend='W20', source='mmm'),
        # PlotDataCdf(runid='120968W34', yname='xkemtm',  xname='xb', zval=0.559, legend='MTM', source='mmm'),
        # PlotDataCdf(runid='120968W34', yname='xkeetgm', xname='xb', zval=0.559, legend='ETGM', source='mmm'),
        # xmax=0.8, xmin=0.02, title_override='120968', allow_title_runid=0,

        # PlotDataCdf(runid='120968W34', yname='xkemmm',  xname='xb', zval=0.559, source='mmm'),
        # PlotDataCdf(runid='120968W34', yname='xkimmm',  xname='xb', zval=0.559, source='mmm'),
        # xmax=0.8, xmin=0.02, title_override='120968', allow_title_runid=0,

        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------
        ## 120982: 0.15 to 0.65
        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------
        ## MMM 8 vs 9 + ETGM
        # PlotDataCdf(runid='120982W31', yname='te', xname='rho', zval=0.6, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='120982W33', yname='te', xname='rho', zval=0.6, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='120982W32', yname='te', xname='rho', zval=0.6, legend='9.0.7'),
        # PlotDataCdf(runid='120982W30', yname='te', xname='rho', zval=0.6, legend='8.2.1'),
        # PlotDataCdf(runid='120982A09', yname='te', xname='rho', zval=0.6, legend='Analysis'),
        # xmax=0.8, title_override='120982 (0.17s)',

        # PlotDataCdf(runid='120982W31', yname='te', xname='time', zval=0.5, legend='9.0.7'), # +ETGM
        # # PlotDataCdf(runid='120982W33', yname='te', xname='time', legend='9.0.7 +Horton'),
        # # PlotDataCdf(runid='120982W32', yname='te', xname='time', legend='9.0.7'),
        # PlotDataCdf(runid='120982W01', yname='te', xname='time', zval=0.5, legend='Only NC'),
        # PlotDataCdf(runid='120982A09', yname='te', xname='time', zval=0.5, legend='Analysis'),
        # title_override=r'120982, $\hat{\rho} = 0.5$'


        # PlotDataCdf(runid='120982W31', yname='te', xname='time', zval=0.5, legend='9.0.7'), # +ETGM
        # # PlotDataCdf(runid='120982W33', yname='te', xname='time', legend='9.0.7 +Horton'),
        # # PlotDataCdf(runid='120982W32', yname='te', xname='time', legend='9.0.7'),
        # PlotDataCdf(runid='120982W01', yname='walltime', xname='time', zval=0.5, legend='Only NC'),
        # PlotDataCdf(runid='120982A09', yname='te', xname='time', zval=0.5, legend='Analysis'),
        # title_override=r'120982'

        # PlotDataCdf(runid='120982W31', yname='te', xname='time', legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='120982W33', yname='te', xname='time', legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='120982W32', yname='te', xname='time', legend='9.0.7'),
        # PlotDataCdf(runid='120982W01', yname='te', xname='time', legend='9.0.7 Disabled'),
        # PlotDataCdf(runid='120982W02', yname='mmmtime', xname='time', legend='9.0.7 Enabled'),
        # PlotDataCdf(runid='120982W01', yname='walltime', xname='time', legend='9.0.7 Disabled'),
        # PlotDataCdf(runid='120982W02', yname='walltime', xname='time', legend='9.0.7 Enabled'),
        # PlotDataCdf(runid='120982W02', yname='cptim', xname='time', legend='9.0.7 Enabled'),
        # title_override='120982' 

        # PlotDataCdf(runid='120982W01', yname='xkemmm', xname='rho', zval=0.65, legend='9.0.7 Disabled'),
        # PlotDataCdf(runid='120982W02', yname='xkemmm', xname='rho', zval=0.65, legend='9.0.7 Enabled'),
        # title_override='120982'

        # PlotDataCdf(runid='120982W31', yname='ti', xname='rho', zval=0.62, legend='MMM 9.0.7 + ETGM'),
        # PlotDataCdf(runid='120982W32', yname='ti', xname='rho', zval=0.62, legend='MMM 9.0.7'),
        # PlotDataCdf(runid='120982A09', yname='ti', xname='rho', zval=0.62, legend='Analysis'),
        # xmax=0.8, title_override='120982',

        # PlotDataCdf(runid='120982W31', yname='xkemmm', xname='rho', zval=0.17, legend='MMM 9.0.7 + ETGM'),
        # PlotDataCdf(runid='120982W32', yname='xkemmm', xname='rho', zval=0.17, legend='MMM 9.0.7'),
        # PlotDataCdf(runid='120982W30', yname='xkemmm', xname='rho', zval=0.17, legend='MMM 8.2.1'),
        # xmax=0.8, title_override='120982 (0.17s)',

        # PlotDataCdf(runid='120982W31', yname='xkemmm', xname='rho', zval=0.62, legend='MMM 9.0.7 + ETGM'),
        # PlotDataCdf(runid='120982W32', yname='xkemmm', xname='rho', zval=0.62, legend='MMM 9.0.7'),
        # xmax=0.8, title_override='120982',

        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------
        ## MTM kyrhos counts
        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------
        # PlotDataCdf(runid='129016Z37', yname='xkemtm', xname='rho', zval=0.35, legend='100'),
        # PlotDataCdf(runid='129016Z38', yname='xkemtm', xname='rho', zval=0.35, legend='200'),
        # PlotDataCdf(runid='129016Z39', yname='xkemtm', xname='rho', zval=0.35, legend='400'),


        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------
        ## 129017: 0.01 to 1
        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------
        # PlotDataCdf(runid='129017W04', yname='te', xname='rho', zval=0.75, legend='MMM 9.0.7 + ETGM'),
        # PlotDataCdf(runid='129017W03', yname='te', xname='rho', zval=0.75, legend='MMM 9.0.7'),
        # PlotDataCdf(runid='129017W02', yname='te', xname='rho', zval=0.75, legend='MMM 8.2.1'),
        # PlotDataCdf(runid='129017W01', yname='te', xname='rho', zval=0.75, legend='Analysis'),
        # xmax=0.8, title_override='129017 (0.25s)',

        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------
        ## 141716: 0.11 to 0.562
        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------
        # PlotDataCdf(runid='141716W01', yname='ti', xname='rho', zval=0.17, legend='Analysis'),
        # PlotDataCdf(runid='141716W04', yname='ti', xname='rho', zval=0.17, legend='MMM 9.0.7 ETGM On'),
        # PlotDataCdf(runid='120982W32', yname='ti', xname='rho', zval=0.17, legend='MMM 9.0.7'),
        # PlotDataCdf(runid='141716W02', yname='ti', xname='rho', zval=0.17, legend='MMM 8.2.1'),
        # xmax=0.8, title_override='120982 (0.17s)',

        # ---------------------------------------------------------------------------
        # ---------------------------------------------------------------------------
        # 138536: 0.04 to 0.79
        # ---------------------------------------------------------------------------
        # ---------------------------------------------------------------------------
        # PlotDataCdf(runid='138536A01', yname='te', xname='rho', zval=0.755, legend='Analysis'),
        # PlotDataCdf(runid='138536W03', yname='te', xname='rho', zval=0.755, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='138536W02', yname='te', xname='rho', zval=0.755, legend='9.0.7'),
        # # PlotDataCdf(runid='138536W01', yname='te', xname='rho', zval=0.17, legend='8.2.1'),
        # xmax=1, title_override='138536',


        # ---------------------------------------------------------------------------
        # ---------------------------------------------------------------------------
        # BEST MATCHES
        # ---------------------------------------------------------------------------
        # ---------------------------------------------------------------------------
        # PlotDataCdf(runid='129016Z36', yname='te', xname='rho', zval=0.4108, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='129016Z43', yname='te', xname='rho', zval=0.4108, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='129016Z33', yname='te', xname='rho', zval=0.4108, legend='9.0.7'),
        # PlotDataCdf(runid='129016A03', yname='te', xname='rho', zval=0.4108, legend='Analysis'),
        # title_override='129016 (0.411s)', xmax=0.8, xmin=0.02, allow_title_time=0,

        # PlotDataCdf(runid='129016W11', yname='te', xname='rho', zval=0.4108, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='129016W10', yname='te', xname='rho', zval=0.4108, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='129016A03', yname='te', xname='rho', zval=0.4108, legend='Analysis'),
        # title_override='129016 (0.411)', xmax=0.8, xmin=0.02, allow_title_time=0,

        # PlotDataCdf(runid='129016W03', yname='walltime', xname='time', zval=0.55, legend=r'9.0.7 $+\chi$'),
        # PlotDataCdf(runid='129016W11', yname='walltime', xname='time', zval=0.55, legend=r'9.0.7 $\pm\chi$'),
        # PlotDataCdf(runid='129016A03', yname='walltime', xname='time', zval=0.55, legend='Analysis'),
        # title_override='129016 (0.55)', xmax=0.8, xmin=0.02, allow_title_time=0,

        # PlotDataCdf(runid='129016W03', yname='ti', xname='rho', zval=0.533, legend='9.0.7 +ETGM'),
        # #PlotDataCdf(runid='129016W04', yname='te', xname='rho', zval=0.533, legend='9.0.7 +Horton'),
        # #PlotDataCdf(runid='129016W02', yname='ti', xname='rho', zval=0.533, legend='9.0.7'),
        # PlotDataCdf(runid='129016A03', yname='ti', xname='rho', zval=0.533, legend='Analysis'),
        # title_override='129016 (0.533s)', xmax=0.8, xmin=0.02, allow_title_time=0,

        # PlotDataCdf(runid='129017W04', yname='te', xname='rho', zval=0.5, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='129017W05', yname='te', xname='rho', zval=0.5, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='129017W03', yname='te', xname='rho', zval=0.5, legend='9.0.7'),
        # PlotDataCdf(runid='129017W02', yname='te', xname='rho', zval=0.5, legend='8.2.1'),
        # PlotDataCdf(runid='129017W01', yname='te', xname='rho', zval=0.5, legend='Analysis'),
        # title_override='129017 (0.500s)', xmax=0.8, xmin=0.02, allow_title_time=0,

        # PlotDataCdf(runid='120968W03', yname='ti', xname='rho', zval=0.511, legend='9.0.7 +ETGM'),
        # # PlotDataCdf(runid='120968W35', yname='te', xname='rho', zval=0.511, legend='9.0.7 +Horton'),
        # # PlotDataCdf(runid='120968W33', yname='te', xname='rho', zval=0.511, legend='9.0.7'),
        # PlotDataCdf(runid='120968A02', yname='ti', xname='rho', zval=0.511, legend='Analysis'),
        # title_override='120968 (0.511s)', xmax=0.8, xmin=0.02, allow_title_time=0,

        # PlotDataCdf(runid='120968W03', yname='te', xname='rho', zval=0.511, legend='Prediction', source='cdf'),
        # PlotDataCdf(runid='120968W02', yname='te', xname='rho', zval=0.511, legend='Analysis', source='cdf'),
        # title_override='120968', xmin=0.02,

        # PlotDataCdf(runid='120982W31', yname='te', xname='rho', zval=0.431, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='120982W33', yname='te', xname='rho', zval=0.431, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='120982W32', yname='te', xname='rho', zval=0.431, legend='9.0.7'),
        # PlotDataCdf(runid='120982A09', yname='te', xname='rho', zval=0.431, legend='Analysis'),
        # title_override='120982 (0.431s)', xmax=0.8, xmin=0.02, allow_title_time=0,

        # PlotDataCdf(runid='120982W31', yname='te', xname='rho', zval=0.431, legend='ETGM Default'),
        # PlotDataCdf(runid='120982W33', yname='te', xname='rho', zval=0.431, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='120982W32', yname='te', xname='rho', zval=0.431, legend='9.0.7'),
        # PlotDataCdf(runid='120982A09', yname='te', xname='rho', zval=0.431, legend='Analysis'),
        # title_override='120982 (0.431s)', xmax=0.8, xmin=0.02, allow_title_time=0,

        # PlotDataCdf(runid='138536W03', yname='ti', xname='rho', zval=0.755, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='138536W04', yname='ti', xname='rho', zval=0.755, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='138536W02', yname='ti', xname='rho', zval=0.755, legend='9.0.7'),
        # PlotDataCdf(runid='138536A01', yname='ti', xname='rho', zval=0.755, legend='Analysis'),
        # title_override='138536 (0.755s)', xmax=0.8, xmin=0.02, allow_title_time=0,

        # PlotDataCdf(runid='138536W03', yname='te', xname='rho', zval=0.755, legend='Prediction'),
        # PlotDataCdf(runid='138536A01', yname='te', xname='rho', zval=0.755, legend='Analysis'),
        # title_override='138536 (0.755s)', xmax=1, xmin=0.02, allow_title_time=0,


        # PlotDataCdf(runid='141716W04', yname='te', zval=0.29, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='141716W05', yname='te', zval=0.29, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='141716W03', yname='te', zval=0.29, legend='9.0.7'),
        # #PlotDataCdf(runid='141716W02', yname='ti', zval=0.29, legend='8.2.1'),
        # PlotDataCdf(runid='141716W01', yname='te', zval=0.29, legend='Analysis'),
        # title_override='141716 (0.29s)', allow_title_time=0, xmin=0.02,

        # PlotDataCdf(runid='129016W12', yname='mmmtime', xname='time', legend=r'Default'),
        # PlotDataCdf(runid='129016W15', yname='mmmtime', xname='time', legend=r'fact, min = 0'),
        # PlotDataCdf(runid='129016W19', yname='te', zval=0.533, legend=r'9.0.9 pphi, No Smoothing'),
        # PlotDataCdf(runid='129016W13', yname='ti', xname='time', zval=0.5, legend=r'Analysis'),
        # title_override='129016', allow_title_time=1,

        PlotDataCdf(runid='129016W13', yname='te', xname='time', zval=0.5, legend=r'Analysis'),
        PlotDataCdf(runid='129016W12', yname='te', xname='time', zval=0.5, legend=r'With Artificial Diffusivity'),
        PlotDataCdf(runid='129016W15', yname='te', xname='time', zval=0.5, legend=r'No Artificial Diffusivity'),
        title_override='129016, rho = 0.5'
        # ---------------------------------------------------------------------------
        # ---------------------------------------------------------------------------
        # TOTAL DIFFUSIVITIES AT BEST MATCHEES
        # ---------------------------------------------------------------------------
        # ---------------------------------------------------------------------------
        # PlotDataCdf(runid='129016Z36', yname='xkimmm', xname='rho', zval=0.4108, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='129016Z43', yname='xkimmm', xname='rho', zval=0.4108, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='129016Z33', yname='xkimmm', xname='rho', zval=0.4108, legend='9.0.7'),
        # title_override='129016 (0.411s)', xmax=0.8, xmin=0.02, allow_title_time=0,

        # PlotDataCdf(runid='129016W03', yname='xkimmm', xname='rho', zval=0.533, legend='9.0.7 +ETGM'),
        # #PlotDataCdf(runid='129016W04', yname='xkimmm', xname='rho', zval=0.533, legend='9.0.7 +Horton'),
        # #PlotDataCdf(runid='129016W02', yname='xkimmm', xname='rho', zval=0.533, legend='9.0.7'),
        # title_override='129016 (0.533)', xmax=0.8, xmin=0.02, allow_title_time=0,

        # PlotDataCdf(runid='129016W03', yname='xkemmm', zval=0.533, legend=r'9.0.7 $+\chi$'),
        # PlotDataCdf(runid='129016W12', yname='xkemmm', zval=0.533, legend=r'9.0.7 $\pm\chi$'),
        # PlotDataCdf(runid='129016W15', yname='xkemmm', zval=0.533, legend=r'9.0.7 $\pm\chi$ 0 PT'),
        # # PlotDataCdf(runid='129016W13', yname='xkemmm', xname='rho', zval=0.533, legend=r'Analysis'),
        # title_override='129016', allow_title_time=0,

        # PlotDataCdf(runid='129017W04', yname='xkimmm', xname='rho', zval=0.5, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='129017W05', yname='xkimmm', xname='rho', zval=0.5, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='129017W03', yname='xkimmm', xname='rho', zval=0.5, legend='9.0.7'),
        # PlotDataCdf(runid='129017W02', yname='xkimmm', xname='rho', zval=0.5, legend='8.2.1'),
        # title_override='129017 (0.500s)', xmax=0.8, xmin=0.02, ymax_cutoff=15, allow_title_time=0,

        # PlotDataCdf(runid='120968W34', yname='xkimmm', xname='rho', zval=0.511, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='120968W35', yname='xkimmm', xname='rho', zval=0.511, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='120968W33', yname='xkimmm', xname='rho', zval=0.511, legend='9.0.7'),
        # title_override='120968 (0.511s)', xmax=0.8, xmin=0.02, allow_title_time=0,

        # PlotDataCdf(runid='120982W31', yname='xkimmm', xname='rho', zval=0.431, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='120982W33', yname='xkimmm', xname='rho', zval=0.431, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='120982W32', yname='xkimmm', xname='rho', zval=0.431, legend='9.0.7'),
        # title_override='120982 (0.431s)', xmax=0.8, xmin=0.02, allow_title_time=0,

        # PlotDataCdf(runid='138536W03', yname='xkemmm', xname='rho', zval=0.755, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='138536W04', yname='xkemmm', xname='rho', zval=0.755, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='138536W02', yname='xkemmm', xname='rho', zval=0.755, legend='9.0.7'),
        # title_override='138536 (0.755s)', xmax=0.8, xmin=0.02, allow_title_time=0,

        # PlotDataCdf(runid='141716W04', yname='xkimmm', zval=0.29, legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='141716W05', yname='xkimmm', zval=0.29, legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='141716W03', yname='xkimmm', zval=0.29, legend='9.0.7'),
        # #PlotDataCdf(runid='141716W02', yname='xkimmm', zval=0.29, legend='8.2.1'),
        # title_override='141716 (0.29s)', allow_title_time=0, xmin=0.02,

        # ---------------------------------------------------------------------------
        # ---------------------------------------------------------------------------
        # WALLTIME
        # ---------------------------------------------------------------------------
        # ---------------------------------------------------------------------------

        # PlotDataCdf(runid='129016W12', yname='mmmtime', xname='time', legend=r'9.0.7 Default'),
        # PlotDataCdf(runid='129016W15', yname='mmmtime', xname='time', legend=r'9.0.7 Default, No smoothing'),
        # PlotDataCdf(runid='129016W20', yname='mmmtime', xname='time', legend=r'9.0.9 pphi'),
        # PlotDataCdf(runid='129016W19', yname='mmmtime', xname='time', legend=r'9.0.9 pphi, No Smoothing'),
        # # PlotDataCdf(runid='129016W13', yname='xkemmm', xname='rho', zval=0.533, legend=r'Analysis'),
        # title_override='129016', allow_title_time=0,

        # PlotDataCdf(runid='129016Z36', yname='mmmtime', xname='time', legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='129016Z44', yname='mmmtime', xname='time', legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='129016Z33', yname='mmmtime', xname='time', legend='9.0.7'),
        # PlotDataCdf(runid='129016Z29', yname='mmmtime', xname='time', legend='8.2.1'),
        # title_override='129016', allow_title_time=0,

        # PlotDataCdf(runid='129016W03', yname='mmmtime', xname='time', legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='129016W04', yname='mmmtime', xname='time', legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='129016W02', yname='mmmtime', xname='time', legend='9.0.7'),
        # PlotDataCdf(runid='129016W01', yname='mmmtime', xname='time', legend='8.2.1'),
        # title_override='129016', allow_title_time=0,

        # PlotDataCdf(runid='129017W04', yname='mmmtime', xname='time', legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='129017W05', yname='mmmtime', xname='time', legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='129017W03', yname='mmmtime', xname='time', legend='9.0.7'),
        # PlotDataCdf(runid='129017W02', yname='mmmtime', xname='time', legend='8.2.1'),
        # title_override='129017', allow_title_time=0,

        # PlotDataCdf(runid='120968W34', yname='mmmtime', xname='time', legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='120968W35', yname='mmmtime', xname='time', legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='120968W33', yname='mmmtime', xname='time', legend='9.0.7'),
        # PlotDataCdf(runid='120968W32', yname='mmmtime', xname='time', legend='8.2.1'),
        # title_override='120968', allow_title_time=0,

        # PlotDataCdf(runid='120982W31', yname='mmmtime', xname='time', legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='120982W33', yname='mmmtime', xname='time', legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='120982W32', yname='mmmtime', xname='time', legend='9.0.7'),
        # PlotDataCdf(runid='120982W30', yname='mmmtime', xname='time', legend='8.2.1'),
        # title_override='120982', allow_title_time=0,

        # PlotDataCdf(runid='138536W03', yname='mmmtime', xname='time', legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='138536W04', yname='mmmtime', xname='time', legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='138536W02', yname='mmmtime', xname='time', legend='9.0.7'),
        # PlotDataCdf(runid='138536W01', yname='mmmtime', xname='time', legend='8.2.1'),
        # title_override='138536', allow_title_time=0,

        # PlotDataCdf(runid='141716W04', yname='mmmtime', xname='time', legend='9.0.7 +ETGM'),
        # PlotDataCdf(runid='141716W05', yname='mmmtime', xname='time', legend='9.0.7 +Horton'),
        # PlotDataCdf(runid='141716W03', yname='mmmtime', xname='time', legend='9.0.7'),
        # PlotDataCdf(runid='141716W02', yname='mmmtime', xname='time', legend='8.2.1'),
        # title_override='141716', allow_title_time=0,


        # ---------------------------------------------------------------------------
        # ---------------------------------------------------------------------------
        # MTM KYRHOS SCANS
        # ---------------------------------------------------------------------------
        # ---------------------------------------------------------------------------
        # PlotDataCdf(runid='129016Z37', yname='walltime', xname='time', legend='100'),
        # PlotDataCdf(runid='129016Z38', yname='walltime', xname='time', legend='200'),
        # PlotDataCdf(runid='129016Z39', yname='walltime', xname='time', legend='400'),
        # title_override='129016', allow_title_time=0,

    )


    main(fig_data, savefig=0, savedata=savedata)


