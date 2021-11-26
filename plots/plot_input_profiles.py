# Standard Packages
from os.path import exists, dirname
import os
import sys
sys.path.insert(0, '../')
# 3rd Party Packages
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.ticker import ScalarFormatter
# Local Packages
from mmm_package import constants, utils

# Formatting list for plot lines
LINE_FORMATS = [{'color': constants.BLUE, 'ls': '-', 'lw': 1.5},
                {'color': constants.RED, 'ls': '-.', 'lw': 1.5},
                {'color': constants.GREEN, 'ls': '--', 'lw': 1.5},
                {'color': constants.PURPLE, 'ls': ':', 'lw': 2},
                {'color': constants.ORANGE, 'ls': '-.', 'lw': 1.75},
                {'color': constants.YELLOW, 'ls': '--', 'lw': 1.75}]

def init_subplots(input_options):
    runid = input_options.runid
    shot_type = input_options.shot_type
    time = input_options.time
    points = input_options.interp_points

    # Init figure and subplots
    fig, axs = plt.subplots(3, 2)
    fig.set_size_inches(8.5, 11)

    # Set plot layout properties (property values chosen for the saved PDF and not for the shown figure)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.2, hspace=0.38, left=0.1, right=0.9, bottom=0.08, top=0.82)
    
    # Set figure title and subtitle
    plt.figtext(0.5, 0.92, 'MMM Input Profiles Using Smoothed Input Parameters', fontsize=14, ha='center')
    plt.figtext(0.5, 0.90, '{0} Shot {1}, Measurement Time {2}s, {3} Input Points'
        .format(shot_type, runid, time, points), fontsize=10, ha='center')

    return fig, axs

def set_axes_style(ax, title, xlabel, ylabel, ylim=None):
    # font sizes for plot titles, labels, and tick values
    titlesize, labelsize, ticksize = 11, 10, 9

    # Set plot font properties, cmr10 = Computer Modern
    ax.set_title(title, fontsize=titlesize)
    ax.set_xlabel(xlabel, labelpad=2, fontsize=labelsize)
    ax.set_ylabel(ylabel, labelpad=2, fontsize=labelsize)
    ax.tick_params(direction='in', labelsize=ticksize)
    ax.xaxis.get_offset_text().set_fontsize(ticksize)
    ax.yaxis.get_offset_text().set_fontsize(ticksize)

    # ScalarFormatter used to modify tick value exponent display
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits(lims=(-2, 3))
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)

    # Set plot limits TODO: add option to adjust ylim range using yvar indices
    ax.set(xlim=(0, 1))

    # Font properties needed to set tick and legend labels
    tick_props = fm.FontProperties(size=ticksize)
    legend_props = fm.FontProperties(size=labelsize)

    # Spines are the frame surrounding each subplot
    for spine in ax.spines:
        ax.spines[spine].set_linewidth(0.5)
    for label in ax.get_xticklabels() :
        label.set_fontproperties(tick_props)
    for label in ax.get_yticklabels() :
        label.set_fontproperties(tick_props)

    ax.legend(borderpad=0, labelspacing=0, frameon=False, prop=legend_props)

def set_axes_plots(ax, time_idx, xvar, *yvars):
    for i, yvar in enumerate(yvars):
        ax.plot(xvar, yvar.values[:, time_idx], **LINE_FORMATS[i % len(LINE_FORMATS)], label=yvar.label)

def make_plots(vars, input_options):
    # Get the index of the measurement time
    time_idx = input_options.time_idx

    # x-axis parameter
    rho = vars.rho.values[:, time_idx]

    # Set plot parameters (unicode_minus does not work in Computer Modern font)
    plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'cm',
        'font.serif': 'cmr10',
        'axes.unicode_minus': False
    })

    # First figure
    fig, axs = init_subplots(input_options)

    set_axes_plots(axs[0, 0], time_idx, rho, vars.te, vars.ti, vars.q)
    set_axes_style(axs[0, 0], r'Temperatures, Safety Factor', r'$\rho$', r'$q, T$ (keV)')

    set_axes_plots(axs[0, 1], time_idx, rho, vars.ne, vars.ni, vars.nf, vars.nz)
    set_axes_style(axs[0, 1], r'Densities', r'$\rho$', r'$n$ $\left(\mathrm{N/m}^3\right)$')

    set_axes_plots(axs[1, 0], time_idx, rho, vars.gte, vars.gti, vars.gq)
    set_axes_style(axs[1, 0], r'Temperatures, Safety Factor Gradients', r'$\rho$', r'$g$')

    set_axes_plots(axs[1, 1], time_idx, rho, vars.gne, vars.gni, vars.gnz)
    set_axes_style(axs[1, 1], r'Density Gradients', r'$\rho$', r'$g$')

    set_axes_plots(axs[2, 0], time_idx, rho, vars.wexbs)
    set_axes_style(axs[2, 0], r'$E \times B$ Shear Rate', r'$\rho$', r'(rad/s)')

    set_axes_plots(axs[2, 1], time_idx, rho, vars.elong)
    set_axes_style(axs[2, 1], r'Elongation', r'$\rho$', '')

    fig.savefig(utils.get_temp_path("input_profiles_1.pdf"))

    # Second figure
    fig, axs = init_subplots(input_options)

    set_axes_plots(axs[0, 0], time_idx, rho, vars.tau)
    set_axes_style(axs[0, 0], r'Temperatures Ratio', r'$\rho$', r'$\tau$')

    set_axes_plots(axs[0, 1], time_idx, rho, vars.beta, vars.betae)
    set_axes_style(axs[0, 1], r'Plasma Betas', r'$\rho$', r'$\beta$')

    set_axes_plots(axs[1, 0], time_idx, rho, vars.etae, vars.etai)
    set_axes_style(axs[1, 0], r'Gradient Ratios', r'$\rho$', r'$\eta$')

    set_axes_plots(axs[1, 1], time_idx, rho, vars.nuste, vars.nusti)
    set_axes_style(axs[1, 1], r'Collisionalities', r'$\rho$', r'$\nu$')

    set_axes_plots(axs[2, 0], time_idx, rho, vars.shear, vars.shat)
    set_axes_style(axs[2, 0], r'Magnetic Shear', r'$\rho$', r'$s$')

    set_axes_plots(axs[2, 1], time_idx, rho, vars.alphamhd)
    set_axes_style(axs[2, 1], r'MHD Alpha', r'$\rho$', r'$\alpha_\mathrm{MHD}$')

    fig.savefig(utils.get_temp_path("input_profiles_2.pdf"))

    # Third figure
    fig, axs = init_subplots(input_options)

    set_axes_plots(axs[0, 0], time_idx, rho, vars.vpol)
    set_axes_style(axs[0, 0], r'Poloidal Velocity', r'$\rho$', r'$\nu$ (m/s)')

    set_axes_plots(axs[0, 1], time_idx, rho, vars.vtor)
    set_axes_style(axs[0, 1], r'Toroidal Velocity', r'$\rho$', r'$\nu$ (m/s)')

    set_axes_plots(axs[1, 0], time_idx, rho, vars.gvpol, vars.gvtor)
    set_axes_style(axs[1, 0], r'Velocity Gradients', r'$\rho$', r'$g$')

    axs[1, 1].axis('off')
    axs[2, 0].axis('off')
    axs[2, 1].axis('off')

    fig.savefig(utils.get_temp_path("input_profiles_3.pdf"))

    # Show figures
    # plt.show()

    utils.merge_input_profile_sheets(input_options)
    # TODO: use pdftek to merge individual profile pdfs 
    # "C:\Users\metxc\Documents\MMM-Package Old\CDF\..\MMM Package\pdftk\pdftk" *.pdf cat output "132017T01 Profiles, 2.100s (08).pdf"


if __name__ == '__main__':
    pass
