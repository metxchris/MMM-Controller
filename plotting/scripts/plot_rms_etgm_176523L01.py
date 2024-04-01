#!/usr/bin/python3

# Standard Packages
import sys; sys.path.insert(0, '../'), sys.path.insert(0, '../../')

# 3rd Party Packages
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# Local Packages
import modules.options
import modules.datahelper as datahelper
import modules.variables as variables
import modules.utils as utils
from modules.enums import SaveType
from plotting.modules.plotstyles import PlotStyles, StyleType
from plotting.plot_variables import FigData, PlotData, main
from plot_rms_mtm import load_all_rho_rms, load_all_rho, get_plot_data


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
        'legend.fontsize': 8,
        # 'text.usetex': True,
    })

    # Define settings for the plot
    fig_data = FigData(
        replace_offset_text=0,
        allow_title_runid=1,
        allow_title_time=0,
        allow_title_factor=1,
        allow_title_rho=1,
        invert_y_axis=False,
        invert_x_axis=False,
        nomralize_y_axis=False,
        nomralize_x_axis=False,
        title_override='D3D',
        ylabel_override=r'RMS $(\Gamma_{\rm e})$ ',
        xlabel_override='',
        ymax=0.1,
        # ymax=3,
        ymin=0,
        savefig=False,
        savedata=False,
    )

    np.set_printoptions(precision=4)

    vname, runid = 'fte', '176523L01'
    title = f'{runid}'  # Plot title
    INCLUDE_AVG_RMS = False

    # Base variable to compare RMS against
    v0, options = load_all_rho(runid, vname, 10011, True)  # E 1E6
    # v0, options = load_all_rho(runid, 10011, True)  # L 1E7

    # Load time variable (x-axis)
    time = variables.InputVariables().time
    time.set(values=options.scan_range)

    # Load variables to compute RMS (y-axis)
    rms01 = load_all_rho_rms(runid, vname, v0, 10012, r'$n_{\rm E} = 50$')
    rms02 = load_all_rho_rms(runid, vname, v0, 10013, r'$n_{\rm E} = 100$')
    rms03 = load_all_rho_rms(runid, vname, v0, 10014, r'$n_{\rm E} = 200$')
    rms07 = load_all_rho_rms(runid, vname, v0, 10018, r'$n_{\rm E} = 400$')
    rms04 = load_all_rho_rms(runid, vname, v0, 10015, r'$\overline{n}_{\rm \mathcal{E}} = 28.0$')
    rms05 = load_all_rho_rms(runid, vname, v0, 10016, r'$n_{\rm L} = 1000$')
    rms06 = load_all_rho_rms(runid, vname, v0, 10017, r'$n_{\rm L} = 2000$')
    rms08 = load_all_rho_rms(runid, vname, v0, 10019, r'$n_{\rm L} = 4000$')
    rms09 = load_all_rho_rms(runid, vname, v0, 10020, r'$n_{\rm L} = 8000$')
    rms10 = load_all_rho_rms(runid, vname, v0, 10021, r'$n_{\rm E} = 800$')
    rms11 = load_all_rho_rms(runid, vname, v0, 10022, r'$n_{\rm E} = 1600$')


    # Exponential Counts (50, 100, 200)
    # data_list = [
    #     get_plot_data(options, time, rms03),
    #     get_plot_data(options, time, rms07),
    #     get_plot_data(options, time, rms10),
    #     get_plot_data(options, time, rms04),
    # ]
    # fig_data.set(ymax=0.02)

    # Linear Counts (1250, 2500, 5000)
    # data_list = [
    #     get_plot_data(options, time, rms05),
    #     get_plot_data(options, time, rms06),
    #     get_plot_data(options, time, rms08),
    #     get_plot_data(options, time, rms09),
    # ]
    # fig_data.set(ymax=0.02)

    ##Linear vs Exponential
    data_list = [
        get_plot_data(options, time, rms09),
        get_plot_data(options, time, rms11),
        get_plot_data(options, time, rms04),
    ]
    fig_data.set(ymax=0.002)

    # Nested Linear vs Exponential
    # data_list = [
    #     get_plot_data(options, time, rms13),
    #     get_plot_data(options, time, rms12),
    # ]
    # fig_data.set(ymax=0.0001)

    # Call plotting script
    fig_data.set(*data_list)
    main(fig_data)

