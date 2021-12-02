# Standard Packages
from dataclasses import dataclass
from cycler import cycler
import copy
import sys
sys.path.insert(0, '../')
# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt
# Local Packages
from main import constants, utils, calculate_inputs
from plots.styles import standard as ps
import settings

# Subplot row and column counts
ROWS, COLS = ps.ROWS, ps.COLS

@dataclass
class PlotData:
    title: str
    xvar: np.ndarray
    yvars: list

class PlotType:
    input = 'Input'
    output = 'Output'
    compared = 'Compared'
    additional = 'Additional'

# Initializes the figure and subplots
def init_figure(input_options, profile_type):
    runid = input_options.runid
    shot_type = input_options.shot_type
    time = input_options.time
    points = input_options.interp_points

    # Init figure and subplots
    fig, axs = plt.subplots(ROWS, COLS)
    
    # Set figure title and subtitle
    modifier = 'Smoothed' if settings.APPLY_SMOOTHING else 'Unaltered'
    plt.figtext(*ps.TITLEPOS, f'MMM {profile_type} Profiles Using {modifier} Input Profiles',
        fontsize=15, ha='center')
    plt.figtext(*ps.SUBTITLEPOS, f'{shot_type} Shot {runid}, Measurement Time {time}s, {points} Input Points', 
        fontsize=10, ha='center')

    return fig, axs

# Creates one individual plot on the specified axis
def make_plot(ax, data, plot_type, time_idx=None):
    for i, yvar in enumerate(data.yvars):
        xvals = data.xvar.values if time_idx is None else data.xvar.values[:, time_idx]
        yvals = yvar.values if time_idx is None else yvar.values[:, time_idx]
        ax.plot(xvals, yvals, label=yvar.label)

    ax.set(title=data.title, xlabel=data.xvar.label, ylabel=data.yvars[0].units_label, xlim=(0, 1))
    ax.axis('on')

    # Check for ylim adjustment (needed when y-values are nearly constant)
    ymax = yvar.values[:, time_idx].max()
    ymin = yvar.values[:, time_idx].min()
    if round(ymax - ymin, 3) == 0:
        ax.set(ylim=(ymin - 1, ymax + 1))

    if plot_type != PlotType().output:
        ax.legend() # Legend disabled for output type plots

# Creates plots for each PlotData() object defined in the input plotdata list
def run_plotting_loop(plotdata, input_options, plot_type):
    print(f'Creating {plot_type.lower()} profile figures...')

    for i, data in enumerate(plotdata):

        # Logic to count (row, col) by col first, then by row; (0, 0), (0, 1), (0, 2), (1, 0), etc.
        row = int(i / COLS) % ROWS
        col = i % COLS

        # Create a new figure when we're on the first subplot
        if row == 0 and col == 0:
            fig, axs = init_figure(input_options, plot_type)

            # Disable all subplot axes until they are used
            for sub_axs in axs:
                for ax in sub_axs:
                    ax.axis('off')

        # Create subplot and enable axis.  Setting data to None will leave the subplot position empty
        if data is not None:
            if plot_type in [PlotType().input, PlotType().compared, PlotType().additional]:
                make_plot(axs[row, col], data, plot_type, input_options.time_idx)
            elif plot_type == PlotType().output:
                make_plot(axs[row, col], data, plot_type)

        # Figure is full of subplots, so save the sheet
        if (i + 1) % (ROWS * COLS) == 0:
            fig.savefig(utils.get_temp_path(f'{plot_type.lower()}_profiles_{int((i + 1) / 6)}.pdf'))


    # Save any remaining subplots to one final sheet
    if (i + 1) % (ROWS * COLS) != 0:
       fig.savefig(utils.get_temp_path(f'{plot_type.lower()}_profiles_{int((i + 1) / 6) + 1}.pdf'))

    # Merge individual pdf sheets with pdftk, then open file (may only open on Windows OS)
    utils.open_file(utils.merge_profile_sheets(input_options, plot_type))

    # Clear plots from memory
    plt.close('all')

def plot_input_profiles(vars, input_options):
    plotdata = [
        PlotData('Temperatures, Safety Factor', vars.rho, [vars.te, vars.ti, vars.q]),
        PlotData(r'$T, q$ Gradients', vars.rho, [vars.gte, vars.gti, vars.gq]),
        PlotData(vars.zeff.name, vars.rho, [vars.zeff]),
        PlotData(vars.btor.name, vars.rho, [vars.btor]),
        PlotData(vars.elong.name, vars.rho, [vars.elong]),
        PlotData(vars.rmaj.name, vars.rho, [vars.rmaj]),
        PlotData('Densities', vars.rho, [vars.ne, vars.ni, vars.nf, vars.nd]),
        PlotData(vars.nz.name, vars.rho, [vars.nz]),
        PlotData('Density Gradients', vars.rho, [vars.gne, vars.gni, vars.gnz]),
        PlotData(vars.nh.name, vars.rho, [vars.nh]),
        PlotData(vars.gnh.name, vars.rho, [vars.gnh]),
        PlotData(vars.wexbs.name, vars.rho, [vars.wexbs]),
        PlotData(vars.vpol.name, vars.rho, [vars.vpol]),
        PlotData(vars.vtor.name, vars.rho, [vars.vtor]),
        PlotData(vars.vpar.name, vars.rho, [vars.vpar]),
        PlotData(vars.gvpol.name, vars.rho, [vars.gvpol]),
        PlotData(vars.gvtor.name, vars.rho, [vars.gvtor]),
        PlotData(vars.gvpar.name, vars.rho, [vars.gvpar]),
        PlotData(vars.ahyd.name, vars.rho, [vars.ahyd]),
        PlotData(vars.aimass.name, vars.rho, [vars.aimass]),
        PlotData(vars.aimp.name, vars.rho, [vars.aimp]),
        PlotData(vars.zimp.name, vars.rho, [vars.zimp]),]

    run_plotting_loop(plotdata, input_options, PlotType().input)

def plot_additional_profiles(vars, input_options):
    plotdata = [
        PlotData(vars.tau.name, vars.rho, [vars.tau]),
        PlotData(vars.beta.name, vars.rho, [vars.beta, vars.betae]),
        PlotData('Gradient Ratios', vars.rho, [vars.etae, vars.etai]),
        PlotData('Collisionalities', vars.rho, [vars.nuste, vars.nusti]),
        PlotData('Magnetic Shear', vars.rho, [vars.shear, vars.shat]),
        PlotData(vars.alphamhd.name, vars.rho, [vars.alphamhd])]

    run_plotting_loop(plotdata, input_options, PlotType().additional)

def plot_output_profiles(vars, input_options):
    plotdata = [
        PlotData(vars.xti.name, vars.rho, [vars.xti]),
        PlotData(vars.xdi.name, vars.rho, [vars.xdi]),
        PlotData(vars.xte.name, vars.rho, [vars.xte]),
        PlotData(vars.xdz.name, vars.rho, [vars.xdz]),
        PlotData(vars.xvt.name, vars.rho, [vars.xvt]),
        PlotData(vars.xvp.name, vars.rho, [vars.xvp]),
        PlotData(vars.xtiW20.name, vars.rho, [vars.xtiW20]),
        PlotData(vars.xdiW20.name, vars.rho, [vars.xdiW20]),
        PlotData(vars.xteW20.name, vars.rho, [vars.xteW20]),
        PlotData(vars.xtiDBM.name, vars.rho, [vars.xtiDBM]),
        PlotData(vars.xdiDBM.name, vars.rho, [vars.xdiDBM]),
        PlotData(vars.xtiDBM.name, vars.rho, [vars.xtiDBM]),
        PlotData(vars.xteDBM.name, vars.rho, [vars.xteDBM]),
        PlotData(vars.xteETG.name, vars.rho, [vars.xteETG]),
        PlotData(vars.xteMTM.name, vars.rho, [vars.xteMTM]),
        PlotData(vars.xteETGM.name, vars.rho, [vars.xteETGM]),
        PlotData(vars.xdiETGM.name, vars.rho, [vars.xdiETGM]),
        PlotData(vars.gmaW20ii.name, vars.rho, [vars.gmaW20ii]),
        PlotData(vars.gmaW20ie.name, vars.rho, [vars.gmaW20ie]),
        PlotData(vars.gmaW20ei.name, vars.rho, [vars.gmaW20ei]),
        PlotData(vars.gmaW20ee.name, vars.rho, [vars.gmaW20ee]),
        PlotData(vars.omgW20ii.name, vars.rho, [vars.omgW20ii]),
        PlotData(vars.omgW20ie.name, vars.rho, [vars.omgW20ie]),
        PlotData(vars.omgW20ei.name, vars.rho, [vars.omgW20ei]),
        PlotData(vars.omgW20ee.name, vars.rho, [vars.omgW20ee]),
        PlotData(vars.gmaDBM.name, vars.rho, [vars.gmaDBM]),
        PlotData(vars.omgDBM.name, vars.rho, [vars.omgDBM]),
        PlotData(vars.gmaMTM.name, vars.rho, [vars.gmaMTM]),
        PlotData(vars.omgMTM.name, vars.rho, [vars.omgMTM]),
        PlotData(vars.gmaETGM.name, vars.rho, [vars.gmaETGM]),
        PlotData(vars.omgETGM.name, vars.rho, [vars.omgETGM]),
        PlotData(vars.dbsqprf.name, vars.rho, [vars.dbsqprf])]

    run_plotting_loop(plotdata, input_options, PlotType().output)

# Compares profiles of calculated values with values found in the CDF
def plot_profile_comparison(cdf_vars, input_vars, input_options):
    # Get list of variables that were both calculated and found in the CDF
    calculated_vars_list = calculate_inputs.get_calculated_vars()
    cdf_var_list = cdf_vars.get_cdf_variables()
    compare_list = [var for var in calculated_vars_list if var in cdf_var_list]

    plotdata = []

    # Automatically build plotdata list with variables that we want to compare
    for var_name in compare_list:

        # Make deep copies since we are modifying the labels below
        cdf_var = copy.deepcopy(getattr(cdf_vars, var_name))
        calc_var = copy.deepcopy(getattr(input_vars, var_name))

        # Skip this variable if there are any issues
        if cdf_var.values is None or cdf_var.values.ndim != calc_var.values.ndim:
            continue

        cdf_var.label += f' ({cdf_var.cdfvar})'
        calc_var.label += ' (Calc)'

        plotdata.append(PlotData(cdf_var.name, input_vars.rho, [cdf_var, calc_var]))

    run_plotting_loop(plotdata, input_options, PlotType().compared)

if __name__ == '__main__':
    # For testing purposes
    print(plt.rcParams.keys())
    pass
