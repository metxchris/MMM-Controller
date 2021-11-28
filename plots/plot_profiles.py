# Standard Packages
import sys
sys.path.insert(0, '../')
# 3rd Party Packages
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.ticker import ScalarFormatter
# Local Packages
from main import constants, utils
import settings

# Formatting list for plot lines
LINE_FORMATS = [{'color': constants.BLUE, 'ls': '-', 'lw': 1.5},
                {'color': constants.RED, 'ls': '-.', 'lw': 1.5},
                {'color': constants.GREEN, 'ls': '--', 'lw': 1.5},
                {'color': constants.PURPLE, 'ls': ':', 'lw': 2},
                {'color': constants.ORANGE, 'ls': '-.', 'lw': 1.75},
                {'color': constants.YELLOW, 'ls': '--', 'lw': 1.75}]

def init_subplots(input_options, profile_type):
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
    modifier = 'Smoothed' if settings.APPLY_SMOOTHING else 'Unaltered'
    plt.figtext(0.5, 0.92, 'MMM {0} Profiles Using {1} Input Profiles'
        .format(profile_type, modifier), fontsize=14, ha='center')
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

def set_axes_input_plots(ax, time_idx, xvar, *yvars):
    for i, yvar in enumerate(yvars):
        ax.plot(xvar, yvar.values[:, time_idx], **LINE_FORMATS[i % len(LINE_FORMATS)], label=yvar.label)

def set_axes_output_plots(ax, xvar, *yvars):
    for i, yvar in enumerate(yvars):
        ax.plot(xvar.values, yvar.values, **LINE_FORMATS[i % len(LINE_FORMATS)], label=yvar.label)

# Set plot parameters 
def set_rcparams():
    plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'cm',
        'font.serif': 'cmr10',
        'axes.unicode_minus': False # (unicode_minus does not work in Computer Modern font)
    })

def plot_input_profiles(vars, input_options):
    set_rcparams()

    print('Creating input profile figures...')

    # Get the index of the measurement time
    time_idx = input_options.time_idx

    # x-axis parameter
    rho = vars.rho.values[:, time_idx]

    # First figure
    fig, axs = init_subplots(input_options, 'Input')

    set_axes_input_plots(axs[0, 0], time_idx, rho, vars.te, vars.ti, vars.q)
    set_axes_style(axs[0, 0], r'Temperatures, Safety Factor', r'$\rho$', r'$q, T$ (keV)')

    set_axes_input_plots(axs[0, 1], time_idx, rho, vars.ne, vars.ni, vars.nf, vars.nz)
    set_axes_style(axs[0, 1], r'Densities', r'$\rho$', r'$n$ $\left(\mathrm{N/m}^3\right)$')

    set_axes_input_plots(axs[1, 0], time_idx, rho, vars.gte, vars.gti, vars.gq)
    set_axes_style(axs[1, 0], r'Temperatures, Safety Factor Gradients', r'$\rho$', r'$g$')

    set_axes_input_plots(axs[1, 1], time_idx, rho, vars.gne, vars.gni, vars.gnz)
    set_axes_style(axs[1, 1], r'Density Gradients', r'$\rho$', r'$g$')

    set_axes_input_plots(axs[2, 0], time_idx, rho, vars.wexbs)
    set_axes_style(axs[2, 0], r'$E \times B$ Shear Rate', r'$\rho$', r'(rad/s)')

    set_axes_input_plots(axs[2, 1], time_idx, rho, vars.elong)
    set_axes_style(axs[2, 1], r'Elongation', r'$\rho$', '')

    fig.savefig(utils.get_temp_path("input_profiles_1.pdf"))

    # Second figure
    fig, axs = init_subplots(input_options, 'Input')

    set_axes_input_plots(axs[0, 0], time_idx, rho, vars.tau)
    set_axes_style(axs[0, 0], r'Temperatures Ratio', r'$\rho$', r'$\tau$')

    set_axes_input_plots(axs[0, 1], time_idx, rho, vars.beta, vars.betae)
    set_axes_style(axs[0, 1], r'Plasma Betas', r'$\rho$', r'$\beta$')

    set_axes_input_plots(axs[1, 0], time_idx, rho, vars.etae, vars.etai)
    set_axes_style(axs[1, 0], r'Gradient Ratios', r'$\rho$', r'$\eta$')

    set_axes_input_plots(axs[1, 1], time_idx, rho, vars.nuste, vars.nusti)
    set_axes_style(axs[1, 1], r'Collisionalities', r'$\rho$', r'$\nu$')

    set_axes_input_plots(axs[2, 0], time_idx, rho, vars.shear, vars.shat)
    set_axes_style(axs[2, 0], r'Magnetic Shear', r'$\rho$', r'$s$')

    set_axes_input_plots(axs[2, 1], time_idx, rho, vars.alphamhd)
    set_axes_style(axs[2, 1], r'MHD Alpha', r'$\rho$', r'$\alpha_\mathrm{MHD}$')

    fig.savefig(utils.get_temp_path("input_profiles_2.pdf"))

    # Third figure
    fig, axs = init_subplots(input_options, 'Input')

    set_axes_input_plots(axs[0, 0], time_idx, rho, vars.vpol)
    set_axes_style(axs[0, 0], r'Poloidal Velocity', r'$\rho$', r'$\nu$ (m/s)')

    set_axes_input_plots(axs[0, 1], time_idx, rho, vars.vtor)
    set_axes_style(axs[0, 1], r'Toroidal Velocity', r'$\rho$', r'$\nu$ (m/s)')

    set_axes_input_plots(axs[1, 0], time_idx, rho, vars.gvpol, vars.gvtor)
    set_axes_style(axs[1, 0], r'Velocity Gradients', r'$\rho$', r'$g$')

    axs[1, 1].axis('off')
    axs[2, 0].axis('off')
    axs[2, 1].axis('off')

    fig.savefig(utils.get_temp_path("input_profiles_3.pdf"))

    # Merge individual pdf sheets with pdftk
    merged_pdf = utils.merge_profile_sheets(input_options, 'Input')

    # Open merged pdf (May only work on Windows)
    utils.open_file(merged_pdf)

    """
    plt.show() can be used to view the individual MatPlotLib figures, 
    but doing so pauses code execution until the figures are closed
    """
    # plt.show()

def plot_output_profiles(vars, input_options):
    set_rcparams()

    print('Creating output profile figures...')

    # First figure
    fig, axs = init_subplots(input_options, 'Output')

    set_axes_output_plots(axs[0, 0], vars.rho, vars.xti, vars.xdi, vars.xte, vars.xdz)
    set_axes_style(axs[0, 0], r'xti, xdi, xte, xdz', r'$\rho$', r'(m/s)')

    set_axes_output_plots(axs[0, 1], vars.rho, vars.xvt, vars.xvp)
    set_axes_style(axs[0, 1], r'xvt, xvp', r'$\rho$', r'$\left(\mathrm{m}^2/\mathrm{s}\right)$')

    set_axes_output_plots(axs[1, 0], vars.rho, vars.xtiW20, vars.xdiW20, vars.xteW20)
    set_axes_style(axs[1, 0], r'xti, xdi, xte (W20)', r'$\rho$', r'$\left(\mathrm{m}^2/\mathrm{s}\right)$')

    set_axes_output_plots(axs[1, 1], vars.rho, vars.xtiDBM, vars.xdiDBM, vars.xteDBM)
    set_axes_style(axs[1, 1], r'xti, xdi, xte (DBM)', r'$\rho$', r'$\left(\mathrm{m}^2/\mathrm{s}\right)$')

    set_axes_output_plots(axs[2, 0], vars.rho, vars.xteETG, vars.xteMTM)
    set_axes_style(axs[2, 0], r'xte, (ETG, MTM)', r'$\rho$', r'$\left(\mathrm{m}^2/\mathrm{s}\right)$')

    set_axes_output_plots(axs[2, 1], vars.rho, vars.xteETGM, vars.xdiETGM)
    set_axes_style(axs[2, 1], r'xte, xdi (ETGM)', r'$\rho$', r'$\left(\mathrm{m}^2/\mathrm{s}\right)$')

    fig.savefig(utils.get_temp_path("output_profiles_1.pdf"))

    # Second figure
    fig, axs = init_subplots(input_options, 'Output')

    set_axes_output_plots(axs[0, 0], vars.rho, vars.gmaW20ii, vars.gmaW20ie, vars.gmaW20ei, vars.gmaW20ee)
    set_axes_style(axs[0, 0], r'gmaW20 (ii, ie, ei, ee)', r'$\rho$', r'$\left(\mathrm{s}^{-1}\right)$')

    set_axes_output_plots(axs[0, 1], vars.rho, vars.omgW20ii, vars.omgW20ie, vars.omgW20ei, vars.omgW20ee)
    set_axes_style(axs[0, 1], r'omgW20 (ii, ie, ei, ee)', r'$\rho$', r'$(\mathrm{rad/s})$')

    set_axes_output_plots(axs[1, 0], vars.rho, vars.gmaDBM, vars.omgDBM)
    set_axes_style(axs[1, 0], r'gma, omg (DBM)', r'$\rho$', r'$\left(\mathrm{s}^{-1}\right), (\mathrm{rad/s})$')

    set_axes_output_plots(axs[1, 1], vars.rho, vars.gmaMTM, vars.omgMTM)
    set_axes_style(axs[1, 1], r'gma, omg (MTM)', r'$\rho$', r'$\left(\mathrm{s}^{-1}\right), (\mathrm{rad/s})$')

    set_axes_output_plots(axs[2, 0], vars.rho, vars.gmaETGM, vars.omgETGM)
    set_axes_style(axs[2, 0], r'gma, omg (ETGM)', r'$\rho$', r'$\left(\mathrm{s}^{-1}\right), (\mathrm{rad/s})$')

    set_axes_output_plots(axs[2, 1], vars.rho, vars.dbsqprf)
    set_axes_style(axs[2, 1], r'dbsqprf', r'$\rho$', r'$\left(\mathrm{s}^{-1}\right)$')

    fig.savefig(utils.get_temp_path("output_profiles_2.pdf"))

    # Merge individual pdf sheets with pdftk
    merged_pdf = utils.merge_profile_sheets(input_options, 'Output')

    # Open merged pdf (May only work on Windows)
    utils.open_file(merged_pdf)

    """
    plt.show() can be used to view the individual MatPlotLib figures, 
    but doing so pauses code execution until the figures are closed
    """
    # plt.show()

if __name__ == '__main__':
    pass
