#!/usr/bin/python3

"""Plot growth rate thresholds from threshold scans

This script is a mess.
"""

# Standard Packages
import sys; sys.path.insert(0, '../'); sys.path.insert(0, '../../')
import logging

# 3rd Party Packages
import numpy as np
from PyQt5.QtGui import QImage
from PyQt5.QtWidgets import QApplication

# Local Packages
import modules.options
import modules.constants as constants
import modules.datahelper as datahelper
import modules.variables as variables
import modules.utils as utils
from modules.enums import SaveType
from plotting.modules.plotstyles import PlotStyles, StyleType
from plotting.plot_variables import FigData, PlotDataCsv, PlotDataCdf, PlotData, main


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
    threshold_pow = 5
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
        # title_override=r'$k_y\rho_\mathrm{e}$ Scan Counts',
        title_override=fr'$\gamma = 10^{{{threshold_pow}}}s^{{-1}}$ Threshold',
        # title_override=r'Min Growth Rate',
        # title_override=r' ',
        ylabel_override='',
        xlabel_override='',
    )

    '''omgse_etgm'''
    # r = '138536A01'
    # n = 463
    # options = modules.options.Options(runid=r, scan_num=n)
    # input_vars, output_vars, __ = datahelper.get_data_objects(options)
    # kyrhos = output_vars.kyrhosETGM.values
    # csound = input_vars.csound.values
    # gne = input_vars.gne.values
    # rmaj = input_vars.rmaj.values

    # input_vars.omgse_etgm.values = kyrhos * csound * gne / rmaj

    # fig_data.set(
    #     PlotData(options, r, 'omgse_etgm', 'rho', input_vars.omgse_etgm, input_vars.rho, is_csv=True,
    #              legend='')
    # )

    '''Thresholds (444 = gne, 445 = gte)'''
    r = '138536A01'
    n_gne = 1803  # gne threshold
    n_gte = 1804  # gte threshold
    n_base = 1784
    threshold = 10**threshold_pow

    options_base = modules.options.Options().load(r, n_base)
    options_gne = modules.options.Options().load(r, n_gne)
    options_gte = modules.options.Options().load(r, n_gte)
    input_vars_base, output_vars_base, __ = datahelper.get_data_objects(options_base)
    input_vars_gne, output_vars_gne, __ = datahelper.get_data_objects(options_gne)
    input_vars_gte, output_vars_gte, __ = datahelper.get_data_objects(options_gte)
    input_vars_dict_gne, output_vars_dict_gne, __ = datahelper.get_all_rho_data(options_gne)
    input_vars_dict_gte, output_vars_dict_gte, __ = datahelper.get_all_rho_data(options_gte)
    num_points_gne = options_gne.input_points
    num_points_gte = options_gte.input_points
    num_scans_gne = options_gne.scan_range.size
    num_scans_gte = options_gte.scan_range.size

    gmaETGM_gne = np.zeros((num_points_gne, num_scans_gne))
    gmaETGM_gte = np.zeros((num_points_gte, num_scans_gte))
    gte = np.zeros_like(gmaETGM_gte)
    gne = np.zeros_like(gmaETGM_gne)
    gte_threshold = np.zeros((num_points_gte))
    gne_threshold = np.zeros((num_points_gne))

    rho_strs = input_vars_dict_gne.keys()

    for i, rho_str in enumerate(rho_strs):
        gmaETGM_gne[i, :] = getattr(output_vars_dict_gne[rho_str], 'gmaETGM').values
        gmaETGM_gte[i, :] = getattr(output_vars_dict_gte[rho_str], 'gmaETGM').values
        gte[i, :] = getattr(input_vars_dict_gte[rho_str], 'gte').values
        gne[i, :] = getattr(input_vars_dict_gne[rho_str], 'gne').values

    gmaETGM_gne_bool = gmaETGM_gne < threshold
    gmaETGM_gte_bool = gmaETGM_gte < threshold

    for i in range(num_points_gne):

        if gmaETGM_gne_bool[i, :].all():
            # Threshold not reached (no threshold)
            gne_threshold[i] = np.nan

        elif (~gmaETGM_gne_bool[i, :]).all():
            # All above threshold (no threshold)
            gne_threshold[i] = np.nan

        elif gmaETGM_gne_bool[i, 0]:
            # Threshold, First value is True
            idx = gmaETGM_gne_bool[i, :].argmin() - 1  # Find first False, then go 1 index lower
            gne_threshold[i] = gne[i, idx]

        else:
            # Threshold, First value is False
            idx = gmaETGM_gne_bool[i, :].argmax()  # Find first True
            gne_threshold[i] = gne[i, idx]

    for i in range(num_points_gte):

        if gmaETGM_gte_bool[i, :].all():
            # Threshold not reached (no threshold)
            gte_threshold[i] = np.nan

        elif (~gmaETGM_gte_bool[i, :]).all():
            # All above threshold (no threshold)
            gte_threshold[i] = np.nan

        elif gmaETGM_gte_bool[i, 0]:
            # Threshold, First value is True
            idx = gmaETGM_gte_bool[i, :].argmin() - 1  # Find first False, then go 1 index lower
            gte_threshold[i] = gte[i, idx]

        else:
            # Threshold, First value is False
            idx = gmaETGM_gte_bool[i, :].argmax()  # Find first True
            gte_threshold[i] = gte[i, idx]

    input_vars_gne.gne_threshold.values = gne_threshold
    input_vars_gte.gte_threshold.values = gte_threshold
    output_vars_gne.gmaETGM.values = gmaETGM_gne.min(axis=1)
    output_vars_gte.gmaETGM.values = gmaETGM_gte.min(axis=1)

    # fig_data.set(
    #     PlotData(options_gne, r, 'gne_threshold', 'rho', input_vars_gne.gne_threshold, input_vars_gne.rho, is_csv=True,
    #              legend=fr'$g_\mathrm{{ne}}$ ($\gamma = 10^{{{threshold_pow}}}s^{{-1}}$ Threshold)'),
    # )

    # fig_data.set(
    #     PlotData(options_gte, r, 'gte_threshold', 'rho', input_vars_gte.gte_threshold, input_vars_gte.rho, is_csv=True,
    #              legend=fr'$g_\mathrm{{Te}}$ ($\gamma = 10^{{{threshold_pow}}}s^{{-1}}$ Threshold)'),
    # )


    # input_vars_base.gne.values *= 0.8
    fig_data.set(
        PlotData(options_gne, r, 'gne_threshold', 'rho', input_vars_gne.gne_threshold, input_vars_gne.rho, is_csv=True,
                 legend=fr'$g_\mathrm{{ne}}$ ($\gamma = 10^{{{threshold_pow}}}s^{{-1}}$ Threshold)'),
        PlotData(options_gne, r, 'shat_gxi', 'rho', input_vars_gne.shat_gxi, input_vars_gne.rho, is_csv=True),
    )

    # fig_data.set(
    #     PlotData(options_gne, r, 'gne_threshold', 'rho', input_vars_gne.gne_threshold, input_vars_gne.rho, is_csv=True,
    #              legend=fr'$g_\mathrm{{ne}}$ ($\gamma = 10^{{{threshold_pow}}}s^{{-1}}$ Threshold)'),
    #     # PlotData(options_gte, r, 'gne', 'rho', input_vars_gte.gne, input_vars_gte.rho, is_csv=True),
    # )

    # fig_data.set(
    #     PlotData(options_gte, r, 'gte_threshold', 'rho', input_vars_gte.gte_threshold, input_vars_gte.rho, is_csv=True,
    #              legend=fr'$g_\mathrm{{Te}}$ ($\gamma = 10^{{{threshold_pow}}}s^{{-1}}$ Threshold)'),
    #     # PlotData(options_gne, r, 'gte', 'rho', input_vars_gne.gte, input_vars_gne.rho, is_csv=True),
    # )

    # fig_data.set(
    #     PlotData(options_gte, r, 'gte_threshold', 'rho', input_vars_gte.gte_threshold, input_vars_gte.rho, is_csv=True,
    #              legend=fr'$g_\mathrm{{Te}}$ ($\gamma = 10^{{{threshold_pow}}}s^{{-1}}$ Threshold)'),
    #     PlotData(options_gte, r, 'shat_gxi', 'rho', input_vars_gte.shat_gxi, input_vars_gte.rho, is_csv=True, legend=''),
    # )

    # fig_data.set(
    #     PlotData(options_gne, r, 'gne_threshold', 'rho', input_vars_gne.gne_threshold, input_vars_gne.rho, is_csv=True,
    #              legend=fr'$g_\mathrm{{ne, ETGM}}$'),
    #     PlotData(options_gte, r, 'gte_threshold', 'rho', input_vars_gte.gte_threshold, input_vars_gte.rho, is_csv=True,
    #              legend=fr'$g_\mathrm{{Te, ETGM}}$'),
    #     PlotData(options_base, r, 'gtecritETG', 'rho', output_vars_base.gtecritETG, output_vars_base.rho, is_csv=True,
    #              legend=fr'$g_\mathrm{{Te, ETG}}$'),
    #     PlotData(options_base, r, 'gne', 'rho', input_vars_base.gne, input_vars_base.rho, is_csv=True,
    #              legend=fr'0.8$g_{{ne}}$'),
    # )

    # fig_data.set(
    #     PlotData(options, r, 'gmaETGM', 'rho', output_vars.gmaETGM, input_vars.rho, is_csv=True,
    #              legend='')
    # )

    main(fig_data)
