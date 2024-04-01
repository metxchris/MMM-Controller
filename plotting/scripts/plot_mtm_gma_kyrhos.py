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
from plotting.modules.plotstyles import PlotStyles, StyleType
from plotting.plot_variables import FigData, PlotData, main


_TIME_CUTOFF = 0.0005


class PlotID:
    def __init__(self, _id, _label=''):
        self.id = _id
        self.label = _label


def get_rms(exp_values, mmm_values):
    summed_diff = np.sum((exp_values - mmm_values)**2, axis=0)
    summed_exp = np.sum((exp_values)**2, axis=0)
    return (summed_diff / summed_exp)**(0.5)


def get_rms_data(shot, id_exp, id_mmm, yname):

    if id_mmm is None:
        return None

    # Load data from CDF
    __, exp_vars = datahelper.initialize_variables(
        modules.options.Options(runid=f'{shot}{id_exp.id}'), 2)
    __, mmm_vars = datahelper.initialize_variables(
        modules.options.Options(runid=f'{shot}{id_mmm.id}'), 2)

    # Get variable corresponding to yname
    mmm_var = getattr(mmm_vars, yname)
    exp_var = getattr(exp_vars, yname)

    # Interpolate experimental radial grid to predicted radial grid, if needed
    if len(exp_var.values[:, 0]) != len(mmm_var.values[:, 0]):
        set_interp = interp1d(exp_vars.xb.values[:, 0], exp_var.values,
                              kind='quadratic', fill_value="extrapolate", axis=0)
        exp_var.set(values=set_interp(mmm_vars.rho.values[:, 0]))

    # Track indices where experimental and predicted time values overlap
    exp_times = []
    mmm_times = []

    # Get matching time indices
    for i, v in enumerate(mmm_vars.time.values):
        if np.abs(exp_vars.time.values - v).min() < _TIME_CUTOFF:
            mmm_times.append(i)
            exp_times.append(np.abs(exp_vars.time.values - v).argmin())

    # Calculate RMS and set values for plotting
    mmm_vars.rms.values = get_rms(exp_var.values[:, exp_times], mmm_var.values[:, mmm_times])
    mmm_vars.time.values = mmm_vars.time.values[mmm_times]
    mmm_vars.rms.label = id_mmm.label

    return mmm_vars


def set_fig_data(fig_data, yname, id0, *id_list):

    data_list = []

    for idn in id_list:
        if idn is not None:
            rms_data = get_rms_data(shot, id0, idn, yname)

            data_list.append(
                PlotData(rms_data.options, rms_data.options.runid, f'{yname}RMS', '',
                         rms_data.rms, rms_data.time, 0, 0, runname=''),
            )

            # Plot title
            title = f'{shot} {getattr(rms_data, yname).label} RMS'

            print(f'Avg RMS ({idn.label}): {np.average(rms_data.rms.values):.4f}')

    # Set RMS plotting data
    fig_data.set(*data_list, title_override=title)


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
        allow_title_runid=0,
        allow_title_time=0,
        allow_title_factor=1,
        allow_title_rho=1,
        invert_y_axis=False,
        invert_x_axis=False,
        nomralize_y_axis=False,
        nomralize_x_axis=False,
        title_override=' ',
        ylabel_override='',
        xlabel_override='',
        # ymax=1.25,
        # ymax=3,
        # ymin=0,
        savefig=False,
        savedata=False,
    )

    np.set_printoptions(precision=4)

    ynames = ['te', 'ti']

    ## Define discharges for plotting.  id0 should be the experimental or base data
    id0 = id1 = id2 = id3 = id4 = None

    # shot = '129016'
    # id0 = PlotID('A03')
    # id1 = PlotID('Z36', '9.0.7') # +ETGM
    # # id2 = PlotID('Z44', '9.0.7 +Horton')
    # # id3 = PlotID('Z33', '9.0.7')
    # id4 = PlotID('Z29', '8.2.1')

    kyrhos = variables.Variable('kyrhos', label=r'$k_y\rho_s$')
    kyrhos.values = np.genfromtxt('data/z_mtm_kyrhos.dat', delimiter=',', dtype=float)[:, :-1].T

    gma = variables.Variable('gma', label=r'$\gamma\ ({\rm 1/s})$')
    gma.values = np.genfromtxt('data/z_mtm_gma.dat', delimiter=',', dtype=float)[:, :-1].T


        # def __init__(self, options, runid, yname, xname, yvar, xvar, zidx, zval,
        #          timeplot=False, yval_base=None, xval_base=None, source='mmm',
        #          is_cdf=False, is_csv=False, factor_symbol=None, scan_factor=None,
        #          ymult=1, xmult=1, rho_value=None, runname=None, legend=''):

    options = modules.options.Options()
    plot_data = PlotData(
        options=options,
        runid=options.runid,
        yname=gma.name,
        xname=kyrhos.name,
        yvar=gma,
        xvar=kyrhos,
        zidx=30,
        zval=0,
    )

    print(plot_data.xvals)


                         # rms_data.rms, rms_data.time, 0, 0, runname='')

    # Plot title
    title = f'NSTX-Like Discharge'

    fig_data.set(plot_data, title_override=title)

    main(fig_data)

    # for yname in ynames:
    #     set_fig_data(fig_data, yname, id0, id1, id2, id3, id4)
    #     main(fig_data)
