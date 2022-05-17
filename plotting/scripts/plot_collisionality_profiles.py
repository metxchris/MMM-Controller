#!/usr/bin/python3

"""Automated plotting using the plot_variables module

This script automatically creates and save plots of the specified file type,
with the option to save CSV corresponding to the data of each plot as well.
"""

# Standard Packages
import sys; sys.path.insert(0, '../'); sys.path.insert(0, '../../')
import logging

# Third Party Packages
import matplotlib.pyplot as plt

# Local Packages
import modules.utils as utils
from plotting.modules.plotstyles import PlotStyles, StyleType
from plotting.plot_variables import AllPlotData, PlotDataCdf, PlotDataCsv, main


_log = logging.getLogger(__name__)


def plot_profiles(profile_list, all_data_list, savefig=True, savedata=True):
    print('plotting profiles...')
    for profiles, all_data in zip(profile_list, all_data_list):
        for p in profiles:
            all_data.set(*p)
            main(all_data, savefig=savefig, savedata=savedata)


def get_plot_data(yname):
    return (
        PlotDataCdf(r1, z1, yname, legend=l1),
        PlotDataCdf(r2, z2, yname, legend=l2),
        PlotDataCdf(r3, z3, yname, legend=l3),
    )

def get_csv_plot_data(scannum, yname):
    return (
        PlotDataCsv(r1, scannum, yname, legend=l1),
        PlotDataCsv(r2, scannum, yname, legend=l2),
        PlotDataCsv(r3, scannum, yname, legend=l3),
    )


if __name__ == '__main__':
    """Run this module directly to plot variable data stored in CDF or CSV files"""

    utils.init_logging()

    # Initialize visual styles for the plot
    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.AIP,
    )

    # Define settings for the plot
    all_data = AllPlotData(
        replace_offset_text=False,
        allow_title_runid=False,
        allow_title_time=False,
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
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
        # 'text.usetex': True,
    })

    # Discharges, time value
    r1, z1, l1 = '120968A02', 0.56, 'High'
    r2, z2, l2 = '129041A10', 0.49, 'Med'
    r3, z3, l3 = '120982A09', 0.62, 'Low'

    s = 1

    profiles = [
        # Input profiles
        # get_plot_data('te'),
        # get_plot_data('ti'),
        # get_plot_data('q'),
        # get_plot_data('gte'),
        # get_plot_data('gti'),
        # get_plot_data('gq'),
        # get_plot_data('bunit'),
        # get_plot_data('gbunit'),
        # get_plot_data('ne'),
        # get_plot_data('nf'),
        # get_plot_data('nz'),
        # get_plot_data('nh'),
        # get_plot_data('gne'),
        # get_plot_data('gni'),
        # get_plot_data('gnz'),
        # get_plot_data('gnh'),
        # get_plot_data('aimass'),
        # get_plot_data('aimp'),
        # get_plot_data('ahyd'),
        # get_plot_data('zimp'),
        # get_plot_data('zeff'),
        # get_plot_data('elong'),
        # get_plot_data('rmaj'),
        # get_plot_data('rmin'),
        # get_plot_data('gxi'),
        # get_plot_data('wexbs'),
        # get_plot_data('nuste'),
        # get_plot_data('nusti'),
        # get_plot_data('loge'),
        # get_plot_data('nuei'),
        # get_plot_data('tau'),
        # get_plot_data('betaeunit'),
        # get_plot_data('etae'),
        # get_plot_data('shat_gxi'),
        # get_plot_data('gyrfeunit'),
        # get_plot_data('gyrfiunit'),
        # get_plot_data('gyrfiunit'),
        # get_plot_data('lareunit'),
        # get_plot_data('rhosunit'),
        # get_plot_data('vthe'),
        # get_plot_data('vthi'),
        # get_plot_data('csound'),
        # get_plot_data('csound_a'),
        # get_plot_data('wtransit'),
        # get_plot_data('wbounce'),
        # get_plot_data('eps'),
        ## Output Profiles
        get_csv_plot_data(s, 'xte'),
        get_csv_plot_data(s, 'xti'),
        get_csv_plot_data(s, 'xdi'),
        get_csv_plot_data(s, 'xdz'),
        get_csv_plot_data(s, 'xvt'),
        get_csv_plot_data(s, 'xvp'),
        get_csv_plot_data(s, 'xteW20'),
        get_csv_plot_data(s, 'xteETGM'),
        get_csv_plot_data(s, 'xte2ETGM'),
        get_csv_plot_data(s, 'xteMTM'),
        get_csv_plot_data(s, 'xtiW20'),
        get_csv_plot_data(s, 'xdiW20'),
        get_csv_plot_data(s, 'xdiW20'),
        get_csv_plot_data(s, 'gmaW20ii'),
        get_csv_plot_data(s, 'gmaW20ie'),
        get_csv_plot_data(s, 'gmaW20ei'),
        get_csv_plot_data(s, 'gmaW20ee'),
        get_csv_plot_data(s, 'gmaMTM'),
        get_csv_plot_data(s, 'gmaETGM'),
        get_csv_plot_data(s, 'omgW20ii'),
        get_csv_plot_data(s, 'omgW20ie'),
        get_csv_plot_data(s, 'omgW20ei'),
        get_csv_plot_data(s, 'omgW20ee'),
        get_csv_plot_data(s, 'omgMTM'),
        get_csv_plot_data(s, 'omgETGM'),
    ]


    # profiles = [
    #     (PlotDataCdf(r1, z1, 'loge'), PlotDataCdf(r1, z1, 'loge', source='cdf'),),
    #     (PlotDataCdf(r2, z2, 'loge'), PlotDataCdf(r2, z2, 'loge', source='cdf'),),
    #     (PlotDataCdf(r3, z3, 'loge'), PlotDataCdf(r3, z3, 'loge', source='cdf'),),
    # ]

    # MAIN PROFILES
    plot_profiles([profiles], [all_data], savefig=1, savedata=1)
