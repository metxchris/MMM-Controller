# 3rd Party Packages
from matplotlib.pyplot import rcParams

# Local Packages
import plotting.modules.plotstyles


class Dimensions:
    text1_pos: tuple = None
    text2_pos: tuple = None
    text3_pos: tuple = None
    text4_pos: tuple = None
    rows: int = None
    cols: int = None


def init(style):
    Layout = plotting.modules.plotstyles.StyleType.Layout

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
        'legend.frameon': True,  # enables legend border and face
        'legend.loc': 'best',
        'legend.framealpha': 0.6,  # alpha of the border and face
        'legend.edgecolor': '1',  # border color
        'legend.borderpad': 0,  # inner padding of frame
        'legend.borderaxespad': 0.5,  # outer padding of frame
        'legend.labelspacing': 0.1,  # spacing between each label (typically vertical)
        'legend.handlelength': 1.9,  # length of lines
        'legend.handleheight': 1,  # height for each line (and label)
        'legend.handletextpad': 0.5,  # padding between lines and text
        'legend.facecolor': 'inherit',
        'legend.fancybox': False,  # rounded legend frame corners
        'lines.dash_joinstyle': 'round',
        'lines.dash_capstyle': 'butt',
        'patch.linewidth': 0,  # legend frame line width
        'savefig.format': 'pdf',
        'xtick.top': True,
        'ytick.right': True,
    })

    if style is Layout.SINGLE1:
        rcParams.update({
            'axes.formatter.limits': [-1, 2],  # Forces exponent notation below 1e-1 and above 1e2
            'axes.labelpad': 2,
            'axes.labelsize': 7.5,
            'axes.titlesize': 7.5,
            'axes.titlepad': 4,
            'figure.figsize': [2.5, 2.08],
            'figure.subplot.bottom': 0.16,
            'figure.subplot.hspace': 0.38,
            'figure.subplot.left': 0.18,
            'figure.subplot.right': 0.94,
            'figure.subplot.top': 0.9,
            'figure.subplot.wspace': 0.22,
            'legend.fontsize': 7,
            'xtick.labelsize': 7,
            'ytick.labelsize': 7,
        })

    if style is Layout.SINGLE1B:
        rcParams.update({
            'axes.formatter.limits': [-2, 2],  # Forces exponent notation below 1e-1 and above 1e2
            'axes.labelpad': 2,
            'axes.labelsize': 9,
            'axes.titlesize': 8,
            'axes.titlepad': 4,
            'figure.figsize': [2.5, 2.08],
            'figure.subplot.bottom': 0.16,
            'figure.subplot.hspace': 0.38,
            'figure.subplot.left': 0.2,
            'figure.subplot.right': 0.94,
            'figure.subplot.top': 0.9,
            'figure.subplot.wspace': 0.22,
            'legend.fontsize': 8,
            'xtick.labelsize': 7,
            'ytick.labelsize': 7,
        })

    elif style is Layout.SINGLE2:
        rcParams.update({
            'axes.labelpad': 2,
            'axes.labelsize': 8,
            'axes.titlesize': 8,
            'axes.titlepad': 4,
            'figure.figsize': [3, 2.5],
            'figure.subplot.bottom': 0.15,
            'figure.subplot.hspace': 0.38,
            'figure.subplot.left': 0.16,
            'figure.subplot.right': 0.94,
            'figure.subplot.top': 0.9,
            'figure.subplot.wspace': 0.22,
            'legend.fontsize': 7.5,
            'xtick.labelsize': 7.5,
            'ytick.labelsize': 7.5,
        })

    elif style is Layout.SINGLE3:
        rcParams.update({
            'axes.labelpad': 2,
            'axes.labelsize': 8.5,
            'axes.titlesize': 8.5,
            'axes.titlepad': 4,
            'figure.figsize': [3.6, 3],
            'figure.subplot.bottom': 0.15,
            'figure.subplot.hspace': 0.38,
            'figure.subplot.left': 0.16,
            'figure.subplot.right': 0.94,
            'figure.subplot.top': 0.9,
            'figure.subplot.wspace': 0.22,
            'legend.fontsize': 8,
            'xtick.labelsize': 8,
            'ytick.labelsize': 8,
        })

    elif style is Layout.GRID3X2:
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

    elif style is Layout.AIP:
        rcParams.update({
            'axes.formatter.limits': [-2, 2],  # Forces exponent notation below 1e-1 and above 1e2
            'axes.labelpad': 2,
            'axes.labelsize': 10,
            'axes.titlesize': 10,
            'axes.titlepad': 4,
            'figure.figsize': [2.5, 2.08],
            'figure.subplot.bottom': 0.16,
            'figure.subplot.hspace': 0.38,
            'figure.subplot.left': 0.2,
            'figure.subplot.right': 0.94,
            'figure.subplot.top': 0.9,
            'figure.subplot.wspace': 0.22,
            'legend.fontsize': 10,
            'xtick.labelsize': 8,
            'ytick.labelsize': 8,
        })

        Dimensions.text1_pos = (0.5, 0.905)
        Dimensions.text2_pos = (0.5, 0.88)
        Dimensions.text3_pos = (0.5, 0.861)
        Dimensions.text4_pos = (0.5, 0.842)
        Dimensions.rows = 2
        Dimensions.cols = 3

    elif style is Layout.AIP2:
        rcParams.update({
            'axes.formatter.limits': [-2, 2],  # Forces exponent notation below 1e-1 and above 1e2
            'axes.labelpad': 2,
            'axes.labelsize': 10,
            'axes.titlesize': 8.5,
            'axes.titlepad': 3.5,
            'figure.figsize': [2.5, 2.08],
            'figure.subplot.bottom': 0.17,
            'figure.subplot.hspace': 0.38,
            'figure.subplot.left': 0.2,
            'figure.subplot.right': 0.94,
            'figure.subplot.top': 0.9,
            'figure.subplot.wspace': 0.22,
            'legend.fontsize': 10,
            'xtick.labelsize': 8,
            'ytick.labelsize': 8,
        })

        Dimensions.text1_pos = (0.5, 0.905)
        Dimensions.text2_pos = (0.5, 0.88)
        Dimensions.text3_pos = (0.5, 0.861)
        Dimensions.text4_pos = (0.5, 0.842)
        Dimensions.rows = 2
        Dimensions.cols = 3
