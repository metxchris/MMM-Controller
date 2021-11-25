# 3rd Party Packages
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

def init_subplots():
    # Init figure and subplots
    fig, axs = plt.subplots(3, 2)
    fig.set_size_inches(8.5, 11)

    # Set plot layout properties (property values chosen for the saved PDF and not for the shown figure)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.2, hspace=0.38, left=0.1, right=0.9, bottom=0.08, top=0.82)

    # Set fonts (unicode_minus does not work in Computer Modern font)
    plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'cm',
        'font.serif': 'cmr10',
        'axes.unicode_minus': False
    })
    
    # Set figure title and subtitle TODO: Read CDF info from variables
    plt.figtext(0.5, 0.92, 'MMM Input Profiles Using Smoothed Input Parameters', fontsize=14, ha='center')
    plt.figtext(0.5, 0.90, 'DIII-D Shot 132017T01, Measurement Time 2.100s, 200 Input Points', fontsize=10, ha='center')

    return fig, axs

def set_axes_style(ax, title, xlabel, ylabel):
    # font sizes for plot titles, labels, and tick values
    titlesize, labelsize, ticksize = 11, 10, 9

    # Set plot font properties, cmr10 = Computer Modern
    ax.set_title(title, fontsize=titlesize, fontname='cmr10')
    ax.set_xlabel(xlabel, labelpad=2, fontsize=labelsize, fontname='cmr10')
    ax.set_ylabel(ylabel, labelpad=2, fontsize=labelsize, fontname='cmr10')
    ax.tick_params(direction='in', labelsize=ticksize)
    ax.ticklabel_format(style='sci', scilimits=(-2, 3))
    ax.xaxis.get_offset_text().set_size(labelsize)
    ax.yaxis.get_offset_text().set_size(labelsize)
    ax.xaxis.get_offset_text().set_font('cmr10')
    ax.yaxis.get_offset_text().set_font('cmr10')

    # Set plot limits
    ax.set(xlim=(0, 1))

    # Font properties needed to set tick and legend labels
    font_prop = fm.FontProperties(family='cmr10', size=ticksize)

    # Spines are the frame surrounding each subplot
    for spine in ax.spines:
        ax.spines[spine].set_linewidth(0.5)
    for label in ax.get_xticklabels() :
        label.set_fontproperties(font_prop)
    for label in ax.get_yticklabels() :
        label.set_fontproperties(font_prop)

    ax.legend(borderpad=0, labelspacing=0, frameon=False, prop=font_prop)

def set_axes_plots(ax, t_idx, xvar, *yvars):
    # TODO: update format_str with better colors
    format_str = ['b-', 'r-.', 'g--', 'm:', 'c-.', 'y--']
    for i, yvar in enumerate(yvars):
        ax.plot(xvar, yvar.values[:, t_idx], format_str[i % len(format_str)], label=yvar.label)

def make_plots(vars, input_time):
    # x-axis parameter
    rho = vars.rho.values[:, 0]

    # Get the array index and measurement time value corresponding to the input time
    t_idx = vars.get_measurement_time_idx(input_time)
    t = vars.get_measurement_time(t_idx)

    # First figure
    fig, axs = init_subplots()

    set_axes_plots(axs[0, 0], t_idx, rho, vars.te, vars.ti, vars.q)
    set_axes_style(axs[0, 0], r'Temperatures, Safety Factor', r'$\rho$', r'$q, T$ (keV)')

    set_axes_plots(axs[0, 1], t_idx, rho, vars.ne, vars.nf, vars.ni, vars.nz)
    set_axes_style(axs[0, 1], r'Densities', r'$\rho$', r'$n$ $\left(\mathrm{N/m}^3\right)$')

    set_axes_plots(axs[1, 0], t_idx, rho, vars.gte, vars.gti, vars.gq)
    set_axes_style(axs[1, 0], r'Temperatures, Safety Factor Gradients', r'$\rho$', r'$g$')

    set_axes_plots(axs[1, 1], t_idx, rho, vars.gne, vars.gni, vars.gnz)
    set_axes_style(axs[1, 1], r'Density Gradients', r'$\rho$', r'$g$')

    set_axes_plots(axs[2, 0], t_idx, rho, vars.wexbs)
    set_axes_style(axs[2, 0], r'$E \times B$ Shear Rate', r'$\rho$', r'(rad/s)')

    set_axes_plots(axs[2, 1], t_idx, rho, vars.elong)
    set_axes_style(axs[2, 1], r'Elongation', r'$\rho$', '')

    fig.savefig("input_profiles_1.pdf")

    # Second figure
    fig, axs = init_subplots()

    set_axes_plots(axs[0, 0], t_idx, rho, vars.tau)
    set_axes_style(axs[0, 0], r'Temperatures Ratio', r'$\rho$', r'$\tau$')

    set_axes_plots(axs[0, 1], t_idx, rho, vars.beta, vars.betae)
    set_axes_style(axs[0, 1], r'Plasma Betas', r'$\rho$', r'$\beta$')

    set_axes_plots(axs[1, 0], t_idx, rho, vars.etae, vars.etai)
    set_axes_style(axs[1, 0], r'Gradient Ratios', r'$\rho$', r'$\eta$')

    set_axes_plots(axs[1, 1], t_idx, rho, vars.nuste, vars.nusti)
    set_axes_style(axs[1, 1], r'Collisionalities', r'$\rho$', r'$\nu$')

    set_axes_plots(axs[2, 0], t_idx, rho, vars.shear, vars.shat)
    set_axes_style(axs[2, 0], r'Magnetic Shear', r'$\rho$', r'$s$')

    set_axes_plots(axs[2, 1], t_idx, rho, vars.alphamhd)
    set_axes_style(axs[2, 1], r'MHD Alpha', r'$\rho$', r'$\alpha_\mathrm{MHD}$')

    fig.savefig("input_profiles_2.pdf")

    # Third figure
    fig, axs = init_subplots()

    set_axes_plots(axs[0, 0], t_idx, rho, vars.vpol)
    set_axes_style(axs[0, 0], r'Poloidal Velocity', r'$\rho$', r'$\nu$ (m/s)')

    set_axes_plots(axs[0, 1], t_idx, rho, vars.vtor)
    set_axes_style(axs[0, 1], r'Toroidal Velocity', r'$\rho$', r'$\nu$ (m/s)')

    set_axes_plots(axs[1, 0], t_idx, rho, vars.gvpol, vars.gvtor)
    set_axes_style(axs[1, 0], r'Velocity Gradients', r'$\rho$', r'$g$')

    axs[1, 1].axis('off')
    axs[2, 0].axis('off')
    axs[2, 1].axis('off')

    fig.savefig("input_profiles_3.pdf")

    # Show figures
    plt.show()

    # TODO: use pdftek to merge individual profile pdfs 

if __name__ == '__main__':
    pass
