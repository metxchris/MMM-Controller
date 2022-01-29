#!/usr/bin/python3

# Standard Packages
import sys; sys.path.insert(0, '../')
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
from modules.enums import SaveType
from plotting.modules.plotstyles import PlotStyles, StyleType
from plotting.plot_variables import AllPlotData, PlotDataCsv, main


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
        title_override='',
        ylabel_override='',
        xlabel_override='',
    )

    n = 182
    r = '138536A01'
    profiles = [
        # Input profiles
        (PlotDataCsv(r, n, 'te'), PlotDataCsv(r, n, 'ti'), ),
        (PlotDataCsv(r, n, 'q'), ),
        (PlotDataCsv(r, n, 'wexbs'), ),
        (PlotDataCsv(r, n, 'gte'), PlotDataCsv(r, n, 'gti'), ),
        (PlotDataCsv(r, n, 'gq'), ),
        (PlotDataCsv(r, n, 'btor'), ),
        (PlotDataCsv(r, n, 'ne'), PlotDataCsv(r, n, 'nh'), ),
        (PlotDataCsv(r, n, 'nf'), PlotDataCsv(r, n, 'nz'), ),
        (PlotDataCsv(r, n, 'gne'), PlotDataCsv(r, n, 'gnh'), PlotDataCsv(r, n, 'gni'), ),
        (PlotDataCsv(r, n, 'gnz'), ),
        (PlotDataCsv(r, n, 'vpol'), ),
        (PlotDataCsv(r, n, 'vtor'), PlotDataCsv(r, n, 'vpar'), ),
        (PlotDataCsv(r, n, 'gvpol'), ),
        (PlotDataCsv(r, n, 'gvtor'), PlotDataCsv(r, n, 'gvpar'), ),
        (PlotDataCsv(r, n, 'aimass'), ),
        (PlotDataCsv(r, n, 'aimp'), ),
        (PlotDataCsv(r, n, 'ahyd'), ),
        (PlotDataCsv(r, n, 'zimp'), ),
        (PlotDataCsv(r, n, 'zeff'), ),
        (PlotDataCsv(r, n, 'elong'), ),
        (PlotDataCsv(r, n, 'rmaj'), ),
        # Additional profiles
        (PlotDataCsv(r, n, 'tau'), ),
        (PlotDataCsv(r, n, 'beta'), PlotDataCsv(r, n, 'betae'), ),
        (PlotDataCsv(r, n, 'etae'), ),
        (PlotDataCsv(r, n, 'nuste'), PlotDataCsv(r, n, 'nusti'), ),
        (PlotDataCsv(r, n, 'shear'), PlotDataCsv(r, n, 'shat'), PlotDataCsv(r, n, 'shat_gxi'), ),
        (PlotDataCsv(r, n, 'alphamhd'), ),
        (PlotDataCsv(r, n, 'gave'), PlotDataCsv(r, n, 'gave_shat'), PlotDataCsv(r, n, 'gave_shat_gxi'), ),
        (PlotDataCsv(r, n, 'gmax'), ),
        (PlotDataCsv(r, n, 'gyrfi'), ),
        (PlotDataCsv(r, n, 'gyrfe'), ),
        (PlotDataCsv(r, n, 'lare'), ),
        (PlotDataCsv(r, n, 'xetgm_const'), ),
        (PlotDataCsv(r, n, 'vthe'), ),
        (PlotDataCsv(r, n, 'vthi'), ),
        # Individual parameter scan profiles
        (PlotDataCsv(r, n, 'betae'), ),
        (PlotDataCsv(r, n, 'shear'), ),
        (PlotDataCsv(r, n, 'gne'), ),
        (PlotDataCsv(r, n, 'gte'), ),
        (PlotDataCsv(r, n, 'ne'), ),
        (PlotDataCsv(r, n, 'nuei'), ),
        (PlotDataCsv(r, n, 'te'), ),
        (PlotDataCsv(r, n, 'ti'), ),
    ]

    for p in profiles:
        all_data.set(*p)
        main(all_data, autosave=True)

