# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

def plot2d(x1, y1, x2=None, y2=None):
    plt.figure()
    plt.plot(x1,y1,lw=2)
    if x2 is not None and y2 is not None:
        plt.plot(x2, y2, 'o', mew=1.5, fillstyle='none', markersize=4)
    plt.show()

def set_axes_style(ax, title, xlabel, ylabel):
    titlesize, labelsize, ticksize = 11, 10, 9
    ax.set_title(title, fontsize=titlesize, fontname='cmr10')
    ax.set_xlabel(xlabel, labelpad=2, fontsize=labelsize, fontname='cmr10')
    ax.set_ylabel(ylabel, labelpad=2, fontsize=labelsize, fontname='cmr10')
    ax.tick_params(direction='in', labelsize=ticksize)
    ax.xaxis.get_offset_text().set_size(labelsize)
    ax.yaxis.get_offset_text().set_size(labelsize)
    ax.xaxis.get_offset_text().set_font('cmr10')
    ax.yaxis.get_offset_text().set_font('cmr10')
    ax.set(xlim=(0, 1))

    font_prop = fm.FontProperties(family='cmr10', size=ticksize)

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

def plot_input_profiles(vars, input_time):
    # x-axis parameter
    xb = vars.xb.values[:, 0]

    # Get the array index and measurement time value corresponding to the input time
    t_idx = vars.get_measurement_time_idx(input_time)
    t = vars.get_measurement_time(t_idx)

    # Init figure and subplots
    fig, axs = plt.subplots(3, 2)
    fig.set_size_inches(8.5, 11)

    # Set plot layout properties
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.2, hspace=0.38, left=0.1, right=0.9, bottom=0.08, top=0.82)

    # Set fonts
    plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'cm',
        'font.serif': 'cmr10',
        'axes.unicode_minus': False
    })
    
    # Set figure title and subtitle TODO: Read info from variables
    plt.figtext(0.5, 0.92, 'MMM Input Profiles Using Smoothed Input Parameters', fontsize=14, ha='center')
    plt.figtext(0.5, 0.90, 'DIII-D Shot 132017T01, Measurement Time 2.100s, 200 Input Points', fontsize=10, ha='center')

    # Set values for each subplot
    set_axes_plots(axs[0, 0], t_idx, xb, vars.te, vars.ti, vars.q)
    set_axes_style(axs[0, 0], r'Temperatures, Safety Factor', r'$\rho$', r'$q, T$ (kEV)')

    set_axes_plots(axs[0, 1], t_idx, xb, vars.ne, vars.nf, vars.ni, vars.nz)
    set_axes_style(axs[0, 1], r'Densities', r'$\rho$', r'$n$ (N/m$^3$)')

    set_axes_plots(axs[1, 0], t_idx, xb, vars.gte, vars.gti, vars.gq)
    set_axes_style(axs[1, 0], r'Temperatures, Safety Factor Gradients', r'$\rho$', r'$g$')

    set_axes_plots(axs[1, 1], t_idx, xb, vars.gne, vars.gni, vars.gnz)
    set_axes_style(axs[1, 1], r'Density Gradients', r'$\rho$', r'$g$')

    set_axes_plots(axs[2, 0], t_idx, xb, vars.wexbs)
    set_axes_style(axs[2, 0], r'$E \times B$ Shear Rate', r'$\rho$', r'rad/s')

    set_axes_plots(axs[2, 1], t_idx, xb, vars.elong)
    set_axes_style(axs[2, 1], r'Elongation', r'$\rho$', '')

    # Save figures as pdf
    fig.savefig("foo.pdf")
    plt.show()

if __name__ == '__main__':
    pass
