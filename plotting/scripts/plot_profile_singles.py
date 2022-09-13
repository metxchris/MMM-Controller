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
from plotting.plot_variables import AllPlotData, PlotDataCsv, PlotDataCdf, main


_log = logging.getLogger(__name__)


def plot_profiles(profile_list, all_data_list, savefig=True, savedata=True):
    print('plotting profiles...')
    for profiles, all_data in zip(profile_list, all_data_list):
        for p in profiles:
            print(p[0].yname)
            all_data.set(*p)
            main(all_data, savefig=savefig, savedata=savedata)


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

    n = 458  # Scan numbers
    r = '138536A01'  # Discharge

    profiles = [
        # Input profiles
        # (PlotDataCsv(r, n, 'te'), PlotDataCsv(r, n, 'ti'), ),
        # (PlotDataCsv(r, n, 'q'), ),
        # (PlotDataCsv(r, n, 'wexb'), ),
        # (PlotDataCsv(r, n, 'gte'), PlotDataCsv(r, n, 'gti'), ),
        # (PlotDataCsv(r, n, 'gq'), ),
        # (PlotDataCsv(r, n, 'bu'), ),
        # (PlotDataCsv(r, n, 'btor'), ),
        # (PlotDataCsv(r, n, 'bu'), PlotDataCsv(r, n, 'btor'), ),
        # (PlotDataCsv(r, n, 'ne'), PlotDataCsv(r, n, 'nh'), ),
        # (PlotDataCsv(r, n, 'nf'), PlotDataCsv(r, n, 'nz'), ),
        # (PlotDataCsv(r, n, 'gne'), PlotDataCsv(r, n, 'gnh'), PlotDataCsv(r, n, 'gni'), ),
        # (PlotDataCsv(r, n, 'gnz'), ),
        # (PlotDataCsv(r, n, 'gbu'), ),
        # # (PlotDataCsv(r, n, 'vpol'), ),
        # # (PlotDataCsv(r, n, 'vtor'), PlotDataCsv(r, n, 'vpar'), ),
        # # (PlotDataCsv(r, n, 'gvpol'), ),
        # # (PlotDataCsv(r, n, 'gvtor'), PlotDataCsv(r, n, 'gvpar'), ),
        # (PlotDataCsv(r, n, 'ai'), ),
        # (PlotDataCsv(r, n, 'az'), ),
        # (PlotDataCsv(r, n, 'ah'), ),
        # (PlotDataCsv(r, n, 'zz'), ),
        # (PlotDataCsv(r, n, 'zeff'), ),
        # (PlotDataCsv(r, n, 'elong'), ),
        # (PlotDataCsv(r, n, 'rmaj'), ),
        # (PlotDataCsv(r, n, 'rmin'), ),
        # (PlotDataCsv(r, n, 'gxi'), ),
        # # Additional profiles
        # (PlotDataCsv(r, n, 'bu_btor'), ),
        # (PlotDataCsv(r, n, 'tau'), ),
        # (PlotDataCsv(r, n, 'tauh'), ),
        # (PlotDataCsv(r, n, 'loge'), ),
        # (PlotDataCsv(r, n, 'epsilonne'), ),
        # # (PlotDataCsv(r, n, 'beta'), PlotDataCsv(r, n, 'betae'), ),
        # (PlotDataCsv(r, n, 'betaeu'), ),
        # (PlotDataCsv(r, n, 'betapu'), ),
        # (PlotDataCsv(r, n, 'etae'), ),
        # (PlotDataCsv(r, n, 'eps'), ),
        # (PlotDataCsv(r, n, 'nuste'), PlotDataCsv(r, n, 'nusti'), ),
        (PlotDataCsv(r, n, 'shear'), PlotDataCsv(r, n, 'shat_gxi'), ),
        # (PlotDataCsv(r, n, 'alphamhd'), ),
        # (PlotDataCsv(r, n, 'alphamhdu'), ),
        # (PlotDataCsv(r, n, 'gyrfi'), ),
        # (PlotDataCsv(r, n, 'gyrfe'), ),
        # (PlotDataCsv(r, n, 'lare'), ),
        # (PlotDataCsv(r, n, 'gyrfiu'), ),
        # (PlotDataCsv(r, n, 'gyrfeu'), ),
        # (PlotDataCsv(r, n, 'lareu'), ),
        # (PlotDataCsv(r, n, 'rhosu'), ),
        # (PlotDataCsv(r, n, 'xetgm_const'), ),
        # (PlotDataCsv(r, n, 'vthe'), ),
        # (PlotDataCsv(r, n, 'vthi'), ),
        # (PlotDataCsv(r, n, 'csound'), ),
        # (PlotDataCsv(r, n, 'csound_a'), ),
        # (PlotDataCsv(r, n, 'wte'), ),
        # (PlotDataCsv(r, n, 'wbe'), ),
        # # Individual parameter scan profiles
        # (PlotDataCsv(r, n, 'betae'), ),
        # (PlotDataCsv(r, n, 'shat_gxi'), ),
        # (PlotDataCsv(r, n, 'gne'), ),
        # (PlotDataCsv(r, n, 'gte'), ),
        (PlotDataCsv(r, n, 'gte'), PlotDataCsv(r, n, 'gtecETG'), ),
        (PlotDataCsv(r, n, 'ne'), ),
        (PlotDataCsv(r, n, 'nuei'), ),
        (PlotDataCsv(r, n, 'te'), ),
        # Output Profiles
        # (PlotDataCsv(r, n, 'gmaETGM'), ),
        # (PlotDataCsv(r, n, 'omgETGM'), ),
        # (PlotDataCsv(r, n, 'xteETGM'), ),
        # (PlotDataCsv(r, n, 'xte2ETGM'), ),
        # (PlotDataCsv(r, n, 'kyrhoeETGM'), ),
        # (PlotDataCsv(r, n, 'kyrhosETGM'), ),
        (PlotDataCsv(r, n, 'gaveETGM'), ),
        # (PlotDataCsv(r, n, 'phiETGM'), ),
        # (PlotDataCsv(r, n, 'AparaETGM'), ),
        (PlotDataCsv(r, n, 'waETGM'), ),
        (PlotDataCsv(r, n, 'wseETGM'), ),
    ]

    # profiles = [
    #     (PlotDataCdf(r, 0.629, 'tauh'), ),
    # ]

    # MAIN PROFILES
    plot_profiles([profiles], [all_data], savefig=0, savedata=1)
