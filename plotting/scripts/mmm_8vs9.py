# Standard Packages
import sys; sys.path.insert(0, '../'); sys.path.insert(0, '../../')
import logging

# Third Party Packages
from matplotlib.pyplot import rcParams

# Local Packages
import modules.utils as utils
from plotting.modules.plotstyles import PlotStyles, StyleType
from plotting.plot_variables import AllPlotData, PlotDataCdf, PlotDataCsv, main


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
    all_data = AllPlotData(
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

    all_data.set(

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
        # PlotDataCdf(runid='129016A03', yname='ti', xname='rho', zval=0.35, legend='Experiment'),
        # PlotDataCdf(runid='129016Z21', yname='ti', xname='rho', zval=0.35, legend='MMM Disabled'),
        ## ---------------------------------------------------------------------------

        ## MMM 8 vs 9 + ETGM
        # PlotDataCdf(runid='129016Z29', yname='te', xname='rho', zval=0.37, legend='MMM 8.2.1'),
        # PlotDataCdf(runid='129016Z36', yname='te', xname='rho', zval=0.5, legend='MMM 9.0.7'),
        # PlotDataCdf(runid='129016A03', yname='te', xname='rho', zval=0.5, legend='Experiment'),

        # PlotDataCdf(runid='129016Z29', yname='wexb', xname='rho', zval=0.34, legend='MMM 8.2.1'),
        # PlotDataCdf(runid='129016Z36', yname='wexb', xname='rho', zval=0.34, legend='MMM 9.0.7'),
        # PlotDataCdf(runid='129016A03', yname='wexb', xname='rho', zval=0.34, legend='Experiment'),

        # PlotDataCdf(runid='129016Z21', yname='ti', xname='rho', zval=0.35, legend='MMM Disabled'),
        # title_override='129016',
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
        ## 120968
        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------
        ## MMM 8 vs 9 + ETGM
        # PlotDataCdf(runid='120968A02', yname='ti', xname='rho', zval=0.56, legend='Experiment'),
        # PlotDataCdf(runid='120968W34', yname='ti', xname='rho', zval=0.56, legend='MMM 9.0.7 + ETGM'),
        # PlotDataCdf(runid='120968W33', yname='ti', xname='rho', zval=0.56, legend='MMM 9.0.7'),
        # PlotDataCdf(runid='120968W32', yname='te', xname='rho', zval=0.37, legend='MMM 8.2.1'),
        # xmax=0.8, title_override='120968, 0.56s',

        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------
        ## 120982
        ## ---------------------------------------------------------------------------
        ## ---------------------------------------------------------------------------
        ## MMM 8 vs 9 + ETGM
        # PlotDataCdf(runid='120982W31', yname='ti', xname='rho', zval=0.17, legend='MMM 9.0.7 + ETGM'),
        # PlotDataCdf(runid='120982W32', yname='ti', xname='rho', zval=0.17, legend='MMM 9.0.7'),
        # PlotDataCdf(runid='120982W30', yname='ti', xname='rho', zval=0.17, legend='MMM 8.2.1'),
        # PlotDataCdf(runid='120982A09', yname='ti', xname='rho', zval=0.17, legend='Experiment'),
        # xmax=0.8, title_override='120982 (0.17s)',

        # PlotDataCdf(runid='120982W31', yname='ti', xname='rho', zval=0.62, legend='MMM 9.0.7 + ETGM'),
        # PlotDataCdf(runid='120982W32', yname='ti', xname='rho', zval=0.62, legend='MMM 9.0.7'),
        # PlotDataCdf(runid='120982A09', yname='ti', xname='rho', zval=0.62, legend='Experiment'),
        # xmax=0.8, title_override='120982',

        # PlotDataCdf(runid='120982W31', yname='xkemmm', xname='rho', zval=0.17, legend='MMM 9.0.7 + ETGM'),
        # PlotDataCdf(runid='120982W32', yname='xkemmm', xname='rho', zval=0.17, legend='MMM 9.0.7'),
        # PlotDataCdf(runid='120982W30', yname='xkemmm', xname='rho', zval=0.17, legend='MMM 8.2.1'),
        # xmax=0.8, title_override='120982 (0.17s)',

        # PlotDataCdf(runid='120982W31', yname='xkemmm', xname='rho', zval=0.62, legend='MMM 9.0.7 + ETGM'),
        # PlotDataCdf(runid='120982W32', yname='xkemmm', xname='rho', zval=0.62, legend='MMM 9.0.7'),
        # xmax=0.8, title_override='120982',
    )


    main(all_data, savefig=False, savedata=False)


