# 3rd Party Packages
import matplotlib.pyplot as plt


def init():
    # Show all rcParams keys using: print(plt.rcParams.keys())
    plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'cm',
        'mathtext.bf': 'cmb10',
        'font.serif': 'cmr10',
        'axes.unicode_minus': False,  # unicode_minus does not work in Computer Modern font (cm, cmr10)
        'axes.formatter.limits': [-2, 3],  # Forces exponent notation below 10**(-2) and above 10**3
        'axes.grid': True,
        'axes.labelpad': 2,
        'axes.labelsize': 10,
        'axes.titlesize': 10,
        'axes.linewidth': 0.5,
        'axes.facecolor': '#f8f8f8',
        'axes.formatter.use_mathtext': True,
        'figure.dpi': 150.0,
        'figure.figsize': [3.5, 3],
        'figure.subplot.bottom': 0.15,
        'figure.subplot.hspace': 0.38,
        'figure.subplot.left': 0.16,
        'figure.subplot.right': 0.9,
        'figure.subplot.top': 0.9,
        'figure.subplot.wspace': 0.22,
        'grid.alpha': 1.0,
        'grid.color': '#fff',
        'legend.fontsize': 9,
        'legend.frameon': False,
        'legend.borderpad': 0,
        'legend.labelspacing': 0.1,
        'legend.handlelength': 2.5,
        'lines.dash_joinstyle': 'round',
        'lines.dash_capstyle': 'butt',
        'lines.linewidth': 1.0,
        'ytick.direction': 'in',
        'xtick.direction': 'in',
        'xtick.labelsize': 9,
        'xtick.major.size': 0.0,
        'ytick.labelsize': 9,
        'ytick.major.size': 0.0,
    })
