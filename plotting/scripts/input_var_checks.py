# Standard Packages
import sys; sys.path.insert(0, '../'); sys.path.insert(0, '../../')
import logging

# Third Party Packages
import matplotlib.pyplot as plt

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
        allow_title_runid=0,
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
        xmax=0.8,
    )

    fig_data.set(


    )

    fig_data.set(

        # ## 1:
        # PlotDataCdf(runid='129016X32', yname='xkiw20', xname='xb', zval=0.3, legend=r'$B_{\rm \phi}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016X32', yname='btor', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='T', title_override='btor'
        # PlotDataCdf(runid='129016X32', yname='xdew20', xname='xb', zval=0.3, legend=r'$B_{\rm u}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016X32', yname='bu', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='T',
        # PlotDataCdf(runid='129016X32', yname='xkew20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$g_{\rm Bu}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X32', yname='gbu', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',

        ## 2:
        # PlotDataCdf(runid='129016X33', yname='xkiw20', xname='xb', zval=0.3,  ymult=1e-4, legend=r'$g_{\rm ne}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X33', yname='gne', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ', title_override='gne',
        # PlotDataCdf(runid='129016X33', yname='xdew20', xname='xb', zval=0.3,  ymult=1e-4, legend=r'$g_{\rm nh}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X33', yname='gnh', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016X33', yname='xkew20', xname='xb', zval=0.3,  ymult=1e-4, legend=r'$g_{\rm ni}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X33', yname='gni', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',

        ## 3:
        # PlotDataCdf(runid='129016X34', yname='xkiw20', xname='xb',  ymult=1e-4, zval=0.3, legend=r'$g_{\rm nz}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X34', yname='gnz', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016X34', yname='xdew20', xname='xb',  ymult=1e-4, zval=0.3, legend=r'$g_{\rm Te}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X34', yname='gte', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016X34', yname='xkew20', xname='xb',  ymult=1e-4, zval=0.3, legend=r'$g_{\rm Ti}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X34', yname='gti', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',

        ## 4:
        # PlotDataCdf(runid='129016X35', yname='xkiw20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$g_{\rm q}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X35', yname='gq', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016X35', yname='xdew20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$B_{\rm \theta}^\prime$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X35', yname='dbp', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='T', title_override='dbp',
        # PlotDataCdf(runid='129016X35', yname='xkew20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$B_{\rm \theta}^{\prime\prime}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X35', yname='d2bp', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='T', title_override='d2bp',

        ## 5:
        # PlotDataCdf(runid='129016X36', yname='xkiw20', xname='xb', zval=0.3, legend=r'$v_{\rm \phi}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016X36', yname='vtor', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='m/s',
        # PlotDataCdf(runid='129016X36', yname='xdew20', xname='xb', zval=0.3, legend=r'$v_{\rm \theta}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016X36', yname='vpol', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='m/s', title_override=r'$\ \ $vpol'
        # PlotDataCdf(runid='129016X36', yname='xkew20', xname='xb', zval=0.3, legend=r'$v_{\rm \parallel}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016X36', yname='vpar', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='m/s', title_override=r'$\ \ $vpar'

        ## 6:
        # PlotDataCdf(runid='129016X37', yname='xkiw20', xname='xb', zval=0.3,  ymult=1e-4,legend=r'$g_{\rm v\phi}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X37', yname='gvtor', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016X37', yname='xdew20', xname='xb', zval=0.3,  ymult=1e-4,legend=r'$g_{\rm v\theta}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X37', yname='gvpol', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016X37', yname='xkew20', xname='xb', zval=0.3,  ymult=1e-4,legend=r'$g_{\rm v\parallel}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X37', yname='gvpar', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',

        ## 7:
        # PlotDataCdf(runid='129016X38', yname='xkiw20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$\omega_{\rm E \times B}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X38', yname='wexbsmod', xname='xb', zval=0.3, legend=r'$\omega_{\rm E \times B}$ Expected', source=r'raw'),
        # PlotDataCdf(runid='129016X38', yname='wexb', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='1/s',
        # PlotDataCdf(runid='129016X38', yname='xdew20', xname='xb', ymult=1e-4, zval=0.3, legend=r'$Z_{\rm eff}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X38', yname='zeff', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016X38', yname='xkew20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$\nabla\hat{\rho}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X38', yname='gxi', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='1/m', title_override='gxi',

        # PlotDataCdf(runid='129016X44', yname='xkiw20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$\omega_{\rm E \times B}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X44', yname='wexb', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='1/s',

        ## 8:
        # PlotDataCdf(runid='129016X39', yname='xkiw20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$g_\kappa$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X39', yname='gelong', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016X39', yname='xdew20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$\overline{M}_\mathrm{h}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X39', yname='ah', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='u',
        # PlotDataCdf(runid='129016X39', yname='xkew20', xname='xb', zval=0.4, ymult=1e-4, legend=r'$\overline{M}_\mathrm{z}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X39', yname='az', xname='xb', zval=0.4, source=r'mmm'),
        # ylabel_override='u',

        ## 9:
        # PlotDataCdf(runid='129016X40', yname='xkiw20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$\eta_\mathrm{NC}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X40', yname='etanc', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'$\Omega {\rm m}$',
        # PlotDataCdf(runid='129016X40', yname='xdew20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$n_\mathrm{e}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X40', yname='ne', xname='xb', zval=0.3, ymult=1e-19, legend=r'$n_\mathrm{e}$ Calculated', source=r'mmm'),
        # PlotDataCdf(runid='129016X40', yname='ne', xname='xb', zval=0.3, ymult=1e-19, legend=r'$n_\mathrm{e}$ CDF', source=r'cdf'),
        # ylabel_override=r'1/m$^3$', title_override=r'$\ \ $ne',
        # PlotDataCdf(runid='129016X40', yname='xkew20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$n_\mathrm{h}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X40', yname='nh', xname='xb', zval=0.3, ymult=1e-19, source=r'mmm'),
        # ylabel_override=r'1/m$^3$',

        ## 10:
        # PlotDataCdf(runid='129016X41', yname='xkiw20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$n_\mathrm{z}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X41', yname='nz', xname='xb', zval=0.3, ymult=1e-19, source=r'mmm'),
        # ylabel_override=r'1/m$^3$',
        # PlotDataCdf(runid='129016X41', yname='xdew20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$n_\mathrm{f}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X41', yname='nf', xname='xb', zval=0.3, ymult=1e-19, source=r'mmm'),
        # ylabel_override=r'1/m$^3$', title_override=r'$\ \ $nf',
        # PlotDataCdf(runid='129016X41', yname='xkew20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$\kappa$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X41', yname='elong', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'1/m$^3$',

        ## 11:
        # PlotDataCdf(runid='129016X42', yname='xkiw20', xname='xb', zval=0.3, legend=r'$r$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016X42', yname='rmin', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'm',
        # PlotDataCdf(runid='129016X42', yname='xdew20', xname='xb', zval=0.3, legend=r'$R$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016X42', yname='rmaj', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'm',
        # PlotDataCdf(runid='129016X42', yname='xkew20', xname='xb', zval=0.3, legend=r'$R_0$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016X42', yname='rmaj', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'm',

        ## 12:
        # PlotDataCdf(runid='129016X43', yname='xkiw20', xname='xb', zval=0.3, legend=r'$r$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016X43', yname='rmin', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'm',
        # PlotDataCdf(runid='129016X43', yname='xdew20', xname='xb', zval=0.3, legend=r'$T_{\rm e}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016X43', yname='te', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'keV',
        # PlotDataCdf(runid='129016X43', yname='xkew20', xname='xb', zval=0.3, legend=r'$T_{\rm i}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016X43', yname='ti', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'keV',


        # """
        # WEXB CHECKS
        # """
        ## 7:
        # PlotDataCdf(runid='129016X38', yname='xkiw20', xname='xb', zval=0.3, ymult=1e-4, legend=r'MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016X38', yname='wexbsmod', xname='xb', zval=0.3, legend=r'SREXBSMOD', source=r'raw'),
        # PlotDataCdf(runid='129016X38', yname='wexbsmod', xname='xb', zval=0.3, legend=r'SREXBSMOD', source=r'cdf'),
        # PlotDataCdf(runid='129016X38', yname='wexbsv2', xname='xb', zval=0.3, legend=r'SREXBSV2', source=r'raw'),
        # ylabel_override='1/s', title_override=r'$\omega_{\rm E \times B}$'

        # """
        # OLD CHECKS (v9.0.1)
        # """

        ## 1:
        # PlotDataCdf(runid='129016Q87', yname='xkiw20', xname='xb', zval=0.3, legend=r'$B_{\rm \phi}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016Q87', yname='btor', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='T', title_override='btor'
        # PlotDataCdf(runid='129016Q50', yname='xdew20', xname='xb', zval=0.3, legend=r'$B_{\rm u}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016Q50', yname='bu', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='T',
        # PlotDataCdf(runid='129016Q50', yname='xkew20', xname='xb', zval=0.3, legend=r'$g_{\rm Bu}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016Q50', yname='gbu', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',

        ## 2:
        # PlotDataCdf(runid='129016Q84', yname='xkiw20', xname='xb', zval=0.3, legend=r'$g_{\rm ne}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016Q84', yname='gne', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ', title_override='gne',
        # PlotDataCdf(runid='129016Q84', yname='xdew20', xname='xb', zval=0.3, legend=r'$g_{\rm nh}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016Q84', yname='gnh', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016Q84', yname='xkew20', xname='xb', zval=0.3, legend=r'$g_{\rm ni}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016Q84', yname='gni', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',

        ## 3:
        # PlotDataCdf(runid='129016W33', yname='xkiw20', xname='xb', zval=0.3, legend=r'$g_{\rm nz}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W33', yname='gnz', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016W33', yname='xdew20', xname='xb', zval=0.3, legend=r'$g_{\rm Te}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016W33', yname='gte', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016W33', yname='xkew20', xname='xb', zval=0.3, legend=r'$g_{\rm Ti}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W33', yname='gti', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',

        ## 4:
        # PlotDataCdf(runid='129016W34', yname='xkiw20', xname='xb', zval=0.3, legend=r'$g_{\rm q}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W34', yname='gq', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016W34', yname='xdew20', xname='xb', zval=0.3, legend=r'$B_{\rm \theta}^\prime$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W34', yname='dbp', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='T', title_override='dbp',
        # PlotDataCdf(runid='129016W34', yname='xkew20', xname='xb', zval=0.3, legend=r'$B_{\rm \theta}^{\prime\prime}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W34', yname='d2bp', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='T', title_override='d2bp',

        ## 5:
        # PlotDataCdf(runid='129016W35', yname='xkiw20', xname='xb', zval=0.3, legend=r'$v_{\rm \phi}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W35', yname='vtor', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='m/s',
        # PlotDataCdf(runid='129016W35', yname='xdew20', xname='xb', zval=0.3, legend=r'$v_{\rm \theta}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W35', yname='vpol', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='m/s', title_override=r'$\ \ $vpol'
        # PlotDataCdf(runid='129016W35', yname='xkew20', xname='xb', zval=0.3, legend=r'$v_{\rm \parallel}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W35', yname='vpar', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='m/s', title_override=r'$\ \ $vpar'

        ## 6:
        # PlotDataCdf(runid='129016W36', yname='xkiw20', xname='xb', zval=0.3, legend=r'$g_{\rm v\phi}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W36', yname='gvtor', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016W36', yname='xdew20', xname='xb', zval=0.3, legend=r'$g_{\rm v\theta}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W36', yname='gvpol', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016W36', yname='xkew20', xname='xb', zval=0.3, legend=r'$g_{\rm v\parallel}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W36', yname='gvpar', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',

        ## 7:
        # PlotDataCdf(runid='129016W37', yname='xkiw20', xname='xb', zval=0.3, ymult=1e-4, legend=r'$\omega_{\rm E \times B}$ MMM Input', source=r'raw'),
        # PlotDataCdf(runid='129016W37', yname='wexb', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='1/s',
        # PlotDataCdf(runid='129016W37', yname='xdew20', xname='xb', zval=0.3, legend=r'$Z_{\rm eff}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W37', yname='zeff', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016W37', yname='xkew20', xname='xb', zval=0.4, legend=r'$\nabla\hat{\rho}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W37', yname='gxi', xname='xb', zval=0.4, source=r'mmm'),
        # ylabel_override='1/m', title_override='gxi',

        ## 8:
        # PlotDataCdf(runid='129016Q92', yname='xkiw20', xname='xb', zval=0.3, legend=r'$g_\kappa$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016Q92', yname='gelong', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=' ',
        # PlotDataCdf(runid='129016W38', yname='xdew20', xname='xb', zval=0.3, legend=r'$\overline{M}_\mathrm{h}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W38', yname='ah', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override='u',
        # PlotDataCdf(runid='129016W38', yname='xkew20', xname='xb', zval=0.4, legend=r'$\overline{M}_\mathrm{z}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W38', yname='az', xname='xb', zval=0.4, source=r'mmm'),
        # ylabel_override='u',

        ## 9:
        # PlotDataCdf(runid='129016W39', yname='xkiw20', xname='xb', zval=0.3, legend=r'$\eta_\mathrm{NC}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W39', yname='etanc', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'$\Omega {\rm m}$',
        # PlotDataCdf(runid='129016W39', yname='xdew20', xname='xb', zval=0.3, legend=r'$n_\mathrm{e}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W43', yname='xdew20', xname='xb', zval=0.3, legend=r'$n_\mathrm{e}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W39', yname='ne', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'1/m$^3$', title_override=r'$\ \ $ne',
        # PlotDataCdf(runid='129016W39', yname='xkew20', xname='xb', zval=0.3, legend=r'$n_\mathrm{h}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W39', yname='nh', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'1/m$^3$',

        ## 10:
        # PlotDataCdf(runid='129016W40', yname='xkiw20', xname='xb', zval=0.3, legend=r'$n_\mathrm{z}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W40', yname='nz', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'1/m$^3$',
        # PlotDataCdf(runid='129016W40', yname='xdew20', xname='xb', zval=0.3, legend=r'$n_\mathrm{f}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W40', yname='nf', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'1/m$^3$', title_override=r'$\ \ $nf',
        # PlotDataCdf(runid='129016W40', yname='xkew20', xname='xb', zval=0.3, legend=r'$\kappa$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W40', yname='elong', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'1/m$^3$',

        ## 11:
        # PlotDataCdf(runid='129016W41', yname='xkiw20', xname='xb', zval=0.3, legend=r'$r$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W41', yname='rmin', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'm',
        # PlotDataCdf(runid='129016W41', yname='xdew20', xname='xb', zval=0.3, legend=r'$R$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W41', yname='rmaj', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'm',
        # PlotDataCdf(runid='129016W41', yname='xkew20', xname='xb', zval=0.3, legend=r'$R_0$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W41', yname='rmaj', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'm',

        ## 12:
        # PlotDataCdf(runid='129016W42', yname='xkiw20', xname='xb', zval=0.3, legend=r'$r$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W42', yname='rmin', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'm',
        # PlotDataCdf(runid='129016W42', yname='xdew20', xname='xb', zval=0.3, legend=r'$T_{\rm e}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W42', yname='te', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'keV',
        # PlotDataCdf(runid='129016W42', yname='xkew20', xname='xb', zval=0.3, legend=r'$T_{\rm i}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016W42', yname='ti', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'keV',


        ## 13: ne test
        # PlotDataCdf(runid='129016Q89', yname='xkiw20', xname='xb', zval=0.3, legend=r'$n_\mathrm{e}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016Q89', yname='ne', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'1/m$^3$', title_override=r'$\ \ $ne, raw',
        # PlotDataCdf(runid='129016Q89', yname='xdew20', xname='xb', zval=0.3, legend=r'$n_\mathrm{e}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016Q89', yname='ne', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'1/m$^3$', title_override=r'$\qquad\quad$ne/10$^{19}$',
        # PlotDataCdf(runid='129016Q89', yname='xkew20', xname='xb', zval=0.3, legend=r'$n_\mathrm{e}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016Q89', yname='ne', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'1/m$^3$', title_override=r'$\ \ $ne, verified',


        ## 14:
        # PlotDataCdf(runid='129016Q91', yname='xkiw20', xname='xb', zval=0.3, legend=r'$n_\mathrm{h}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016Q91', yname='nh', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'1/m$^3$', title_override=r'$\quad$nh',
        # PlotDataCdf(runid='129016Q91', yname='xdew20', xname='xb', zval=0.3, legend=r'$n_\mathrm{e}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016Q91', yname='ne', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'1/m$^3$', title_override=r'$\quad$ne',
        # PlotDataCdf(runid='129016Q91', yname='xkew20', xname='xb', zval=0.3, legend=r'$n_\mathrm{z}$ MMM Input', source=r'cdf'),
        # PlotDataCdf(runid='129016Q91', yname='nz', xname='xb', zval=0.3, source=r'mmm'),
        # ylabel_override=r'1/m$^3$', title_override=r'$\quad$nz',

    )

    main(fig_data, savefig=False, savedata=False)


