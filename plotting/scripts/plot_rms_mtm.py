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


def get_rms(exp_values, mmm_values):
    '''Calculates RMS of mmm_values against exp_values'''
    summed_diff = np.sum((exp_values - mmm_values)**2, axis=0)
    summed_exp = np.sum((exp_values)**2, axis=0)
    return (summed_diff / summed_exp)**(0.5)


def load_vars(runid, scan_num):
    '''Loads output variables for specified runid, scan_num'''
    vars0 = variables.OutputVariables(modules.options.Options(runid=runid, scan_num=scan_num))
    vars0.load_from_csv(SaveType.OUTPUT)
    return vars0


def load_all_rho(runid, var_to_plot, scan_num, base_var=False):
    '''Loads all rho data for variable var_to_plot'''
    options = modules.options.Options().load(runid, scan_num)

    __, output_vars, __ = datahelper.get_all_rho_data(options, get_input=0, get_control=0)

    rho_strs = list(output_vars.keys())
    var = np.zeros((len(rho_strs), len(getattr(output_vars[rho_strs[0]], var_to_plot).values)))

    for i, rho_str in enumerate(rho_strs):
        var[i, :] = getattr(output_vars[rho_str], var_to_plot).values

    if base_var:
        return var, options

    return var


def load_all_rho_rms(runid, var_to_plot, v0, scan_num, name):
    '''load rms for all rho data of variable var_to_plot'''
    v1 = load_all_rho(runid, var_to_plot, scan_num)
    return get_rms_var(name, v0, v1)


def get_plot_data(options, time, rms):
    print(f'{rms.label:<18} ({np.average(rms.values):.3e})')
    return PlotData(
        options=options,
        runid=options.runid,
        yname=rms.name,
        xname=time.name,
        yvar=rms,
        xvar=time,
        zidx=0,
        zval=0,
    )


def get_rms_var(name, v0, v1):
    rms_values = get_rms(v0, v1)
    # rms_avg_str = utils.get_sci_notation(np.average(rms_values), 2)
    rms = variables.Variable(
        'RMS',
        cdfvar='',
        label=name,
        units='',
    )
    rms.set(values=rms_values, label=f'{name}')
    return rms


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
        title_override='NSTX',
        ylabel_override=r'RMS $(\Gamma_{\rm e})$ ',
        xlabel_override='',
        ymax=0.1,
        # ymax=3,
        ymin=0,
        savefig=False,
        savedata=False,
    )

    np.set_printoptions(precision=4)

    vname, runid = 'fte', '138536A01'
    title = f'{runid}'  # Plot title

    # Base variable to compare RMS against
    v0, options = load_all_rho(runid, vname, 10000, True)  # E 1E6
    # v0, options = load_all_rho(runid, 10011, True)  # L 1E7

    # Load time variable (x-axis)
    time = variables.InputVariables().time
    time.set(values=options.scan_range)

    # Load variables to compute RMS (y-axis)
    rms1 = load_all_rho_rms(runid, vname, v0, 10001, r'$n_{\rm E} = 200$')
    rms2 = load_all_rho_rms(runid, vname, v0, 10009, r'$n_{\rm L} = 8000$')
    rms3 = load_all_rho_rms(runid, vname, v0, 10006, r'$n_{\rm L} = 10000$')
    rms4 = load_all_rho_rms(runid, vname, v0, 10008, r'$n_{\rm S} = 200$')  #mod_exp2 = 1.5
    rms5 = load_all_rho_rms(runid, vname, v0, 10007, r'$n_{\rm S} = 200$')  #mod_exp2 = 2
    rms6 = load_all_rho_rms(runid, vname, v0, 10004, r'$n_{\rm L} = 4000$')
    rms7 = load_all_rho_rms(runid, vname, v0, 10012, r'$n_{\rm L} = 2000$')
    rms8 = load_all_rho_rms(runid, vname, v0, 10013, r'$n_{\rm L} = 16000$')
    rms9 = load_all_rho_rms(runid, vname, v0, 10014, r'$n_{\rm E} = 50$')
    rms10 = load_all_rho_rms(runid, vname, v0, 10015, r'$n_{\rm E} = 100$')
    rms11 = load_all_rho_rms(runid, vname, v0, 10016, r'$n_{\rm E} = 400$')
    rms12 = load_all_rho_rms(runid, vname, v0, 10017, r'$n_{\rm E} = 8000$')
    rms13 = load_all_rho_rms(runid, vname, v0, 10018, r'$n_{\rm M} = 200$')  #mod_min = 0.75, mod_exp = 4
    rms14 = load_all_rho_rms(runid, vname, v0, 10019, r'$n_{\rm M} = 200$')  #mod_min = 0.65, mod_exp = 4
    rms15 = load_all_rho_rms(runid, vname, v0, 10020, r'$n_{\rm M} = 200$')  #mod_min = 0.85, mod_exp = 4

    # Exponential default vs special
    data_list = [
        get_plot_data(options, time, rms1),
        # get_plot_data(options, time, rms4),
        # get_plot_data(options, time, rms5),
        # get_plot_data(options, time, rms13),
        # get_plot_data(options, time, rms14),
        get_plot_data(options, time, rms15),
    ]

    # # Exponential Counts (50, 100, 200)
    # data_list = [
    #     get_plot_data(options, time, rms9),
    #     get_plot_data(options, time, rms10),
    #     get_plot_data(options, time, rms1),
    #     # get_plot_data(options, time, rms11),
    # ]
    # fig_data.set(ymax=0.15)

    # Linear Counts (2000, 4000, 8000)
    # data_list = [
    #     get_plot_data(options, time, rms7),
    #     get_plot_data(options, time, rms6),
    #     get_plot_data(options, time, rms2),
    # ]
    # fig_data.set(ymax=0.15)

    ## Linear vs Exponential
    # data_list = [
    #     get_plot_data(options, time, rms1),
    #     get_plot_data(options, time, rms2),
    # ]
    # fig_data.set(ymax=0.15)

    # Call plotting script
    fig_data.set(*data_list)
    main(fig_data)

