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
from main.enums import PlotType, ShotType
from main.options import Options
from plots.styles import standard as ps


# Subplot row and column counts
ROWS, COLS = ps.ROWS, ps.COLS

@dataclass
class PlotData:
    title: str
    xvar: np.ndarray
    yvars: list


# Initializes the figure and subplots
def init_figure(profile_type, xvar_points):
    runid = Options.instance.runid
    shot_type = Options.instance.shot_type
    time = Options.instance.time_str
    points = xvar_points

    # Init figure and subplots
    fig, axs = plt.subplots(ROWS, COLS)
    
    # Set figure title and subtitle
    modifier = 'Smoothed' if Options.instance.apply_smoothing else 'Unsmoothed'
    title_txt = f'MMM {profile_type.name.capitalize()} Profiles Using {modifier} Input Profiles'
    subtitle_txt = f'{shot_type.name} Shot {runid}, Measurement Time {time}s, {points} Radial Points'
    plt.figtext(*ps.TITLEPOS, title_txt, fontsize=15, ha='center')
    plt.figtext(*ps.SUBTITLEPOS, subtitle_txt, fontsize=10, ha='center')

    return fig, axs

# Creates one individual plot on the specified axis
def make_plot(ax, data, plot_type, time_idx=None):

    xvals = data.xvar.values if time_idx is None else data.xvar.values[:, time_idx]

    for i, yvar in enumerate(data.yvars):
        yvals = yvar.values if time_idx is None else yvar.values[:, time_idx]
        ax.plot(xvals, yvals, label=yvar.label)

    ax.set(title=data.title, xlabel=data.xvar.label, ylabel=data.yvars[0].units_label, xlim=(xvals.min(), xvals.max()))
    ax.axis('on')

    # Check for ylim adjustment (needed when y-values are nearly constant and not nearly 0)
    ymax, ymin = yvals.max(), yvals.min()
    if round(ymax - ymin, 3) == 0 and round(ymax, 3) > 0:
        ax.set(ylim=(ymin - 5, ymax + 5))

    # Legend disabled for output type plots
    if plot_type != PlotType.OUTPUT:
        ax.legend() 

# Creates plots for each PlotData() object defined in the input plotdata list
def run_plotting_loop(plotdata, plot_type):
    from plots.styles import standard as ps
    from plots.colors import mmm

    opts = Options.instance

    print(f'Creating {plot_type.name.lower()} profile figures...')

    for i, data in enumerate(plotdata):

        # Logic to count (row, col) by col first, then by row; (0, 0), (0, 1), (0, 2), (1, 0), etc.
        row = int(i / COLS) % ROWS
        col = i % COLS

        # Create a new figure when we're on the first subplot
        if row == 0 and col == 0:
            fig, axs = init_figure(plot_type, data.xvar.values.shape[0])

            # Disable all subplot axes until they are used
            for sub_axs in axs:
                for ax in sub_axs:
                    ax.axis('off')

        # Create subplot and enable axis.  Setting data to None will leave the subplot position empty
        if data is not None:
            if plot_type in [PlotType.INPUT, PlotType.COMPARED, PlotType.ADDITIONAL]:
                make_plot(axs[row, col], data, plot_type, opts.time_idx)
            elif plot_type == PlotType.OUTPUT:
                make_plot(axs[row, col], data, plot_type)

        # Figure is full of subplots, so save the sheet
        if (i + 1) % (ROWS * COLS) == 0:
            fig.savefig(utils.get_temp_path(f'{plot_type.name.lower()}_profiles_{int((i + 1) / 6)}.pdf'))


    # Save any remaining subplots to one final sheet
    if (i + 1) % (ROWS * COLS) != 0:
       fig.savefig(utils.get_temp_path(f'{plot_type.name.lower()}_profiles_{int((i + 1) / 6) + 1}.pdf'))

    # Merge individual pdf sheets with pdftk, then open file (may only open on Windows OS)
    utils.open_file(utils.merge_profile_sheets(opts.runid, opts.scan_num, plot_type.name.capitalize()))

    # Clear plots from memory
    plt.close('all')

def plot_input_profiles(vars):
    plotdata = [
        PlotData('Temperatures', vars.rho, [vars.te, vars.ti]),
        PlotData(vars.q.name, vars.rho, [vars.q]),
        PlotData(vars.wexbs.name, vars.rho, [vars.wexbs]),
        PlotData(r'Temperature Gradients', vars.rho, [vars.gte, vars.gti]),
        PlotData(vars.gq.name, vars.rho, [vars.gq]),
        PlotData(vars.btor.name, vars.rho, [vars.btor]),
        PlotData('Densities', vars.rho, [vars.ne, vars.ni, vars.nf, vars.nd]),
        PlotData(vars.nz.name, vars.rho, [vars.nz]),
        PlotData(vars.nh.name, vars.rho, [vars.nh]),
        PlotData('Density Gradients', vars.rho, [vars.gne, vars.gni]),
        PlotData(vars.gnz.name, vars.rho, [vars.gnz]),
        PlotData(vars.gnh.name, vars.rho, [vars.gnh]),
        PlotData(vars.vpol.name, vars.rho, [vars.vpol]),
        PlotData(vars.vtor.name, vars.rho, [vars.vtor]),
        PlotData(vars.vpar.name, vars.rho, [vars.vpar]),
        PlotData(vars.gvpol.name, vars.rho, [vars.gvpol]),
        PlotData(vars.gvtor.name, vars.rho, [vars.gvtor]),
        PlotData(vars.gvpar.name, vars.rho, [vars.gvpar]),
        PlotData(vars.ahyd.name, vars.rho, [vars.ahyd]),
        PlotData(vars.aimass.name, vars.rho, [vars.aimass]),
        PlotData(vars.aimp.name, vars.rho, [vars.aimp]),
        PlotData(vars.zimp.name, vars.rho, [vars.zimp]),
        PlotData(vars.zeff.name, vars.rho, [vars.zeff]),
        PlotData(vars.elong.name, vars.rho, [vars.elong]),
        PlotData(vars.rmaj.name, vars.rho, [vars.rmaj])]

    run_plotting_loop(plotdata, PlotType.INPUT)

def plot_additional_profiles(vars):
    plotdata = [
        PlotData(vars.tau.name, vars.rho, [vars.tau]),
        PlotData(vars.beta.name, vars.rho, [vars.beta, vars.betae]),
        PlotData('Gradient Ratios', vars.rho, [vars.etae, vars.etai]),
        PlotData(vars.nuei.name, vars.rho, [vars.nuei]),
        PlotData('Collisionalities', vars.rho, [vars.nuste, vars.nusti]),
        PlotData('Magnetic Shear', vars.rho, [vars.shear, vars.shat]),
        PlotData(vars.alphamhd.name, vars.rho, [vars.alphamhd]),
        PlotData(vars.gave.name, vars.rho, [vars.gave]),
        PlotData(vars.gmax.name, vars.rho, [vars.gmax]),
        PlotData(vars.gyrfi.name, vars.rho, [vars.gyrfi]),
        PlotData(vars.vthe.name, vars.rho, [vars.vthe]),
        PlotData(vars.vthi.name, vars.rho, [vars.vthi])]

    run_plotting_loop(plotdata, PlotType.ADDITIONAL)

def plot_output_profiles(vars):
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

    run_plotting_loop(plotdata, PlotType.OUTPUT)

# Compares profiles of calculated values with values found in the CDF
def plot_profile_comparison(cdf_vars, input_vars):
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

    run_plotting_loop(plotdata, PlotType.COMPARED)


if __name__ == '__main__':
    # For testing purposes
    print(plt.rcParams.keys())
