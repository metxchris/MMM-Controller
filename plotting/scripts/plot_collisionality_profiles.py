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
from plotting.plot_variables import FigData, PlotDataCdf, PlotDataCsv, main


_log = logging.getLogger(__name__)


def plot_profiles(profile_list, fig_data_list, savename='', savefig=True, savedata=True):
    print('plotting profiles...')
    for profiles, fig_data in zip(profile_list, fig_data_list):
        for p in profiles:
            fig_data.set(*p, savename_append=savename)
            main(fig_data, savefig=savefig, savedata=savedata)


def get_plot_data(yname):
    return (
        PlotDataCdf(r1, z1, yname, xname='rmina', legend=l1, input_points=101),
        # PlotDataCdf(r2, z2, yname, xname='rmina', legend=l2, input_points=101),
        PlotDataCdf(r3, z3, yname, xname='rmina', legend=l3, input_points=101),
    )

def get_csv_plot_data(yname):
    return (
        PlotDataCsv(r1, s1, yname, xname='rmina', legend=l1),
        # PlotDataCsv(r2, s2, yname, xname='rmina', legend=l2),
        PlotDataCsv(r3, s3, yname, xname='rmina', legend=l3),
    )


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
        xmax=0.8,
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
        # 'text.usetex': True,
    })

    savename = ''

    # Discharges, time value
    r1, s1, z1, l1 = '120968A02', 13307, 0.56, 'High'
    r2, s2, z2, l2 = '129041A10', 13307, 0.49, 'Med'
    r3, s3, z3, l3 = '120982A09', 13307, 0.62, 'Low'

    savename = ''
    if s1 == 13301:
        savename = 'wexb'
    # elif s1 == 13302:
    #     savename = 'gte0'
    # elif s1 == 13303:
    #     savename = 'gti0'
    # elif s1 == 13304:
    #     savename = 'gne0'
    elif s1 == 13305:
        savename = 'ah1'
    elif s1 == 13306:
        savename = 'ah2'
    elif s1 == 13307:
        savename = 'ah3'

    profiles = [
        ### Input profiles
        # get_plot_data('te'),
        # get_plot_data('ti'),
        # get_plot_data('q'),
        # get_plot_data('gte'),
        # get_plot_data('gti'),
        # get_plot_data('gq'),
        # get_plot_data('btor'),
        # get_plot_data('bu'),
        # get_plot_data('gbu'),
        # get_plot_data('ne'),
        # get_plot_data('nf'),
        # get_plot_data('nz'),
        # get_plot_data('nh'),
        # get_plot_data('gne'),
        # get_plot_data('gni'),
        # get_plot_data('gnz'),
        # get_plot_data('gnh'),
        # get_plot_data('ai'),
        # get_plot_data('az'),
        # get_plot_data('ah'),
        # get_plot_data('zz'),
        # get_plot_data('zeff'),
        # get_plot_data('elong'),
        # get_plot_data('rmaj'),
        # get_plot_data('rmin'),
        # get_plot_data('gxi'),
        # get_plot_data('nuste'),
        # get_plot_data('loge'),
        # get_plot_data('nuei'),
        # get_plot_data('ti_te'),
        # get_plot_data('betaeu'),
        # get_plot_data('etae'),
        # get_plot_data('shear'),
        # get_plot_data('shat'),
        # get_plot_data('shat_gxi'),
        # get_plot_data('gyrfi'),
        # get_plot_data('gyrfiu'),
        # get_plot_data('rhos'),
        # get_plot_data('rhosu'),
        # get_plot_data('vthe'),
        # get_plot_data('vthi'),
        # get_plot_data('csound'),
        # get_plot_data('csound_a'),
        # get_plot_data('wte'),
        # get_plot_data('wbe'),
        # get_plot_data('eps'),
        # get_plot_data('wexb'),
        # get_plot_data('etanc'),
        ### Output Profiles
        get_csv_plot_data('xte'),
        get_csv_plot_data('xti'),
        get_csv_plot_data('xde'),
        get_csv_plot_data('xdz'),
        get_csv_plot_data('xvt'),
        get_csv_plot_data('xvp'),
        get_csv_plot_data('fte'),
        get_csv_plot_data('fti'),
        get_csv_plot_data('fde'),
        get_csv_plot_data('fdz'),
        get_csv_plot_data('fvt'),
        get_csv_plot_data('fvp'),
        get_csv_plot_data('xteW20'),
        get_csv_plot_data('xdeW20'),
        get_csv_plot_data('xtiW20'),
        get_csv_plot_data('xteDBM'),
        get_csv_plot_data('xte2DBM'),
        get_csv_plot_data('xdeDBM'),
        get_csv_plot_data('xde2DBM'),
        get_csv_plot_data('xtiDBM'),
        get_csv_plot_data('xti2DBM'),
        get_csv_plot_data('xteETGM'),
        get_csv_plot_data('xte2ETGM'),
        get_csv_plot_data('xteMTM'),
        get_csv_plot_data('xteETG'),
        get_csv_plot_data('gmaW20i'),
        get_csv_plot_data('gmaW20e'),
        get_csv_plot_data('gmaDBM'),
        get_csv_plot_data('gmaMTM'),
        get_csv_plot_data('gmaETGM'),
        get_csv_plot_data('omgW20i'),
        get_csv_plot_data('omgW20e'),
        get_csv_plot_data('omgDBM'),
        get_csv_plot_data('omgMTM'),
        get_csv_plot_data('omgETGM'),
        # get_csv_plot_data('vcz'),
        # get_csv_plot_data('vct'),
        # get_csv_plot_data('vcp'),
        ### Output Profiles
        # get_csv_plot_data('kyrhosDBM'),
        # get_csv_plot_data('gaveDBM'),
        # get_csv_plot_data('phi2DBM'),
        # get_csv_plot_data('Apara2DBM'),
        # get_csv_plot_data('xtiDBM'),
        # get_csv_plot_data('xteDBM'),
        # get_csv_plot_data('xdeDBM'),
        # get_csv_plot_data('fti'),
        # get_csv_plot_data('fte'),
        # get_csv_plot_data('fde'),
    ]


    # profiles = [
    #     (PlotDataCdf(r1, z1, 'loge'), PlotDataCdf(r1, z1, 'loge', source='cdf'),),
    #     (PlotDataCdf(r2, z2, 'loge'), PlotDataCdf(r2, z2, 'loge', source='cdf'),),
    #     (PlotDataCdf(r3, z3, 'loge'), PlotDataCdf(r3, z3, 'loge', source='cdf'),),
    # ]

    # MAIN PROFILES
    plot_profiles([profiles], [fig_data], savename=savename, savefig=1, savedata=1)
