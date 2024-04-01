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
    v0, options = load_all_rho(runid, vname, 10000, True)  # E 1E6
    # v0, options = load_all_rho(runid, 10011, True)  # L 1E7

    # Load time variable (x-axis)
    time = variables.InputVariables().time
    time.set(values=options.scan_range)

    # Load variables to compute RMS (y-axis)
    rms01 = load_all_rho_rms(runid, vname, v0, 10001, r'$n_{\rm E} = 50$')
    rms02 = load_all_rho_rms(runid, vname, v0, 10002, r'$n_{\rm E} = 100$')
    rms03 = load_all_rho_rms(runid, vname, v0, 10003, r'$n_{\rm E} = 200$')
    rms25 = load_all_rho_rms(runid, vname, v0, 10040, r'$n_{\rm E} = 1000$')
    rms04 = load_all_rho_rms(runid, vname, v0, 10004, r'$n_{\rm L} = 2000$')
    rms05 = load_all_rho_rms(runid, vname, v0, 10005, r'$n_{\rm L} = 4000$')
    rms06 = load_all_rho_rms(runid, vname, v0, 10006, r'$n_{\rm L} = 8000$')
    rms10 = load_all_rho_rms(runid, vname, v0, 10010, r'$n_{\rm L} = 16000$')  #0.90
    rms07 = load_all_rho_rms(runid, vname, v0, 10007, r'$n_{\rm S1} = 200$')  #0.85
    rms08 = load_all_rho_rms(runid, vname, v0, 10008, r'$n_{\rm S2} = 200$')  #0.75
    rms09 = load_all_rho_rms(runid, vname, v0, 10009, r'$n_{\rm S2} = 200$')  #0.90
    # rms13 = load_all_rho_rms(runid, vname, v0, 10025, r'$n_{\rm \mathcal{E}} =  99$')
    rms11 = load_all_rho_rms(runid, vname, v0, 10023, r'$n_{\rm \mathcal{E}} = 149$')
    rms12 = load_all_rho_rms(runid, vname, v0, 10024, r'$n_{\rm \mathcal{E}} = 249$')
    rms15 = load_all_rho_rms(runid, vname, v0, 10030, r'$\overline{n}_{\rm \mathcal{E}} = 243$')  # 6x10  (1.270e-03) [opt enabled]
    rms14 = load_all_rho_rms(runid, vname, v0, 10027, r'$\overline{n}_{\rm \mathcal{E}} = 245$')  # 7x10  (1.219e-03)
    rms16 = load_all_rho_rms(runid, vname, v0, 10031, r'$\overline{n}_{\rm \mathcal{E}} = 246$')  # 8x10  (1.206e-03)
    rms17 = load_all_rho_rms(runid, vname, v0, 10032, r'$\overline{n}_{\rm \mathcal{E}} = 248$')  # 10x10 (1.223e-03)
    rms20 = load_all_rho_rms(runid, vname, v0, 10039, r'$\overline{n}_{\rm \mathcal{E}} = 5\times 10$')  # (1.293e-03) [opt disabled]
    rms21 = load_all_rho_rms(runid, vname, v0, 10035, r'$\overline{n}_{\rm \mathcal{E}} = 6\times 8$')   # (1.229e-03)
    rms22 = load_all_rho_rms(runid, vname, v0, 10036, r'$\overline{n}_{\rm \mathcal{E}} = 7\times 7$')   # (1.202e-03)
    rms23 = load_all_rho_rms(runid, vname, v0, 10037, r'$\overline{n}_{\rm \mathcal{E}} = 8\times 6$')   # (1.205e-03)
    rms24 = load_all_rho_rms(runid, vname, v0, 10038, r'$\overline{n}_{\rm \mathcal{E}} = 10\times 5$')  # (1.223e-03)

    # Exponential default vs special
    # data_list = [
    #     get_plot_data(options, time, rms03),
    #     get_plot_data(options, time, rms08),
    #     get_plot_data(options, time, rms07),
    #     get_plot_data(options, time, rms09),
    # ]
    # fig_data.set(ymax=0.02)

    # # RMS of different scans per layer with optimizations enabled
    # data_list = [
    #     get_plot_data(options, time, rms12),  # 7x7  No Optimizations
    #     get_plot_data(options, time, rms14),  # 7x10
    #     get_plot_data(options, time, rms15),  # 6x10
    #     get_plot_data(options, time, rms16),  # 8x10
    #     get_plot_data(options, time, rms17),  # 10x10
    # ]
    # fig_data.set(ymax=0.006)

    # # RMS of different scans per layer no optimizations
    # data_list = [
    #     get_plot_data(options, time, rms20),  # 5x10
    #     get_plot_data(options, time, rms21),  # 6x8
    #     get_plot_data(options, time, rms22),  # 7x7
    #     get_plot_data(options, time, rms23),  # 8x6
    #     get_plot_data(options, time, rms24),  # 10x5
    # ]
    # fig_data.set(ymax=0.006)

    # # Exponential Counts (50, 100, 200)
    # data_list = [
    #     get_plot_data(options, time, rms01),
    #     get_plot_data(options, time, rms02),
    #     get_plot_data(options, time, rms03),
    #     get_plot_data(options, time, rms15),
    # ]
    # fig_data.set(ymax=0.06)

    # Linear Counts (4000, 8000, 16000)
    # data_list = [
    #     # get_plot_data(options, time, rms04),
    #     get_plot_data(options, time, rms05),
    #     get_plot_data(options, time, rms06),
    #     get_plot_data(options, time, rms10),
    # ]
    # fig_data.set(ymax=0.06)

    # Linear vs Exponential
    data_list = [
        get_plot_data(options, time, rms10),
        get_plot_data(options, time, rms03),
        get_plot_data(options, time, rms25),
        # get_plot_data(options, time, rms12),
        get_plot_data(options, time, rms14),
    ]
    fig_data.set(ymax=0.06)

    # Call plotting script
    fig_data.set(*data_list)
    main(fig_data)
