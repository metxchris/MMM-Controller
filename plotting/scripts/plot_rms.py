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
        layout=StyleType.Layout.AIP2,
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
        ymax=1.25,
        # ymax=3,
        ymin=0,
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

    ## 0.15s Start
    # shot = '129016'
    # id0 = PlotID('A03')
    # id1 = PlotID('W03', '9.0.7') # +ETGM
    # # id2 = PlotID('W04', '9.0.7 +Horton')
    # # id3 = PlotID('W02', '9.0.7')
    # id4 = PlotID('W01', '8.2.1')
    # #id4 = PlotID('W05', 'None')

    ## 0.15s Start
    # shot = '129016'
    # id0 = PlotID('W13')
    # id1 = PlotID('W03', r'9.0.7 $+\chi$')
    # id2 = PlotID('W12', r'9.0.7 $\pm\chi$')
    # # id3 = PlotID('W15', r'9.0.7 $\pm\chi$, 0 PT')
    # id4 = PlotID('W16', r'9.0.9 $\pm\chi$, pphi')

    # Region compare
    # shot = '129016'
    # id0 = PlotID('W13')
    # id1 = PlotID('W25', r'9.0.7 Only Confinement')
    # id2 = PlotID('W26', r'9.0.7 All three regions')

    # 0.15s Start
    # shot = '129016'
    # id0 = PlotID('W13')
    # id1 = PlotID('W12', r'9.0.7 Default')
    # id2 = PlotID('W15', r'9.0.7 Default, No Smoothing')
    # id3 = PlotID('W20', r'9.0.9 pphi')
    # id4 = PlotID('W19', r'9.0.9 pphi, No Smoothing')

    # 0.15s Start
    # shot = '129016'
    # id0 = PlotID('W13')
    # id1 = PlotID('W12', r'9.0.7 Default')
    # id2 = PlotID('W15', r'9.0.7 No Smoothing')

    # 0.15s Start
    # shot = '129016'
    # id0 = PlotID('W13')
    # id2 = PlotID('W12', r'9.0.7 Default')
    # id3 = PlotID('W21', r'9.0.7 $g \geq 1$')
    # id4 = PlotID('W22', r'9.0.7 $g \geq 1$, No Smoothing')

    # # 0.15s Start
    # shot = '129016'
    # id0 = PlotID('W13')
    # id1 = PlotID('W12', r'9.0.7')
    # id2 = PlotID('W24', r'9.0.9')
    # id3 = PlotID('W23', r'9.0.9 $g \geq 1$')
    

    ## TE only
    # shot = '129016'
    # id0 = PlotID('A03')
    # id1 = PlotID('W09', '9.0.7 +ETGM')
    # id2 = PlotID('W10', '9.0.7 +Horton')
    # id3 = PlotID('W08', '9.0.7')
    # id4 = PlotID('W07', '8.2.1')

    ## MTM Kyrhos scans
    # shot = '129016'
    # id0 = PlotID('A03')
    # id1 = PlotID('Z37', '100')
    # id2 = PlotID('Z38', '200')
    # id3 = PlotID('Z39', '400')
    # id4 = None

    shot = '120982'
    id0 = PlotID('A09')
    id1 = PlotID('W31', '9.0.7') # +ETGM
    # id2 = PlotID('W33', '9.0.7 +Horton')
    # id3 = PlotID('W32', '9.0.7')
    id4 = PlotID('W30', '8.2.1')

    # shot = '120982'
    # id0 = PlotID('A09')
    # id1 = PlotID('W31', 'ETGM Default')
    # id2 = PlotID('W05', 'ETGM Alternate')
    # id3 = PlotID('W06', 'ETGM Default Sum')
    # id4 = PlotID('W07', 'ETGM Alternate Sum')

    # shot = '120968'
    # id0 = PlotID('W02')
    # id1 = PlotID('W03', '9.0.7 +ETGM')
    # id2 = PlotID('W35', '9.0.7 +Horton')
    # id3 = PlotID('W33', '9.0.7')
    # id4 = PlotID('W32', '8.2.1')

    # shot = '141716'
    # id0 = PlotID('W01')
    # id1 = PlotID('W04', '9.0.7 +ETGM')
    # id2 = PlotID('W05', '9.0.7 +Horton')
    # id3 = PlotID('W03', '9.0.7')
    # id4 = PlotID('W02', '8.2.1')

    # shot = '138536'
    # id0 = PlotID('A01')
    # id1 = PlotID('W03', '9.0.7 +ETGM')
    # id2 = PlotID('W04', '9.0.7 +Horton')
    # id3 = PlotID('W02', '9.0.7')
    # id4 = PlotID('W01', '8.2.1')

    # shot = '129017'
    # id0 = PlotID('W01')
    # id1 = PlotID('W04', '9.0.7 +ETGM')
    # id2 = PlotID('W05', '9.0.7 +Horton')
    # id3 = PlotID('W03', '9.0.7')
    # id4 = PlotID('W02', '8.2.1')

    ## ----------------------------------------------------------------------    
    ## Different run starts
    # shot = '129016'
    # id0 = PlotID('A03')
    # id1 = PlotID('Z36', 'Start = 0.30s')
    # id2 = PlotID('W03', 'Start = 0.15s')
    # id1 = PlotID('Z44', 'Start = 0.30s')
    # id2 = PlotID('W04', 'Start = 0.15s')

    ## Particle counts
    # shot = '129016'
    # id0 = PlotID('A03')
    # id1 = PlotID('W02', 'N =  8000')
    # id2 = PlotID('W06', 'N = 32000')

    ## TE Only
    # shot = '129016'
    # id0 = PlotID('A03')
    # id1 = PlotID('W09', '9.0.7 +ETGM')
    # id2 = PlotID('W10', '9.0.7 +Horton')
    # id3 = PlotID('W08', '9.0.7')
    # id4 = PlotID('W07', '8.2.1')


    for yname in ynames:
        set_fig_data(fig_data, yname, id0, id1, id2, id3, id4)
        main(fig_data)
