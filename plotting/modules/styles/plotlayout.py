# 3rd Party Packages
from matplotlib.pyplot import rcParams

# Local Packages
import plotting.modules.plotstyles as ps


class Dimensions:
    text1_pos: tuple = None
    text2_pos: tuple = None
    text3_pos: tuple = None
    text4_pos: tuple = None
    rows: int = None
    cols: int = None


def init(style):

    # General rcParams
    rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'cm',
        'mathtext.bf': 'cmb10',
        'font.serif': 'cmr10',
        'figure.dpi': 150.0,
        'axes.linewidth': 0.5,
        'axes.unicode_minus': False,  # unicode_minus does not work in Computer Modern font (cm, cmr10)
        'axes.formatter.limits': [-2, 3],  # Forces exponent notation below 1e-2 and above 1e3
        'axes.formatter.use_mathtext': True,
        'legend.frameon': False,
        'legend.borderpad': 0,
        'legend.borderaxespad': 0.5,
        'legend.labelspacing': 0.1,
        'legend.handlelength': 1.9,
        'legend.handleheight': 0.7,
        'legend.handletextpad': 0.5,
        'lines.dash_joinstyle': 'round',
        'lines.dash_capstyle': 'butt',
    })

    if style is ps.StyleType.Layout.SINGLE:
        rcParams.update({
            'axes.labelpad': 2,
            'axes.labelsize': 9,
            'axes.titlesize': 8.5,
            'figure.figsize': [3.5, 3],
            'figure.subplot.bottom': 0.15,
            'figure.subplot.hspace': 0.38,
            'figure.subplot.left': 0.16,
            'figure.subplot.right': 0.9,
            'figure.subplot.top': 0.9,
            'figure.subplot.wspace': 0.22,
            'legend.fontsize': 8,
            'xtick.labelsize': 8,
            'ytick.labelsize': 8,
        })

    elif style is ps.StyleType.Layout.GRID3X2:
        rcParams.update({
            'axes.labelpad': 2,
            'axes.labelsize': 10,
            'axes.titlesize': 10,
            'figure.figsize': [11, 8.5],
            'figure.subplot.bottom': 0.10,
            'figure.subplot.hspace': 0.38,
            'figure.subplot.left': 0.08,
            'figure.subplot.right': 0.92,
            'figure.subplot.top': 0.78,
            'figure.subplot.wspace': 0.22,
            'legend.fontsize': 11,
            'xtick.labelsize': 9,
            'ytick.labelsize': 9,
        })

        Dimensions.text1_pos = (0.5, 0.905)
        Dimensions.text2_pos = (0.5, 0.88)
        Dimensions.text3_pos = (0.5, 0.861)
        Dimensions.text4_pos = (0.5, 0.842)
        Dimensions.rows = 2
        Dimensions.cols = 3
