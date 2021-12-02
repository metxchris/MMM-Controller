from cycler import cycler
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '../')
sys.path.insert(0, '../../')
from main import constants

line_cycle = (  cycler(color=[constants.BLUE, constants.RED, constants.ORANGE, constants.GREEN, constants.PURPLE, constants.YELLOW])
          # cycler(color=['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c'])
        + cycler(dashes=[(1, 0), (8, 2, 1, 2), (5, 2), (1.5, 0.75), (2, 1), (4, 1)])
        + cycler(lw=[1.75, 1.75, 1.75, 2.0, 2.0, 2.0]))

plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'cm',
        'mathtext.bf': 'cmb10',
        'font.serif': 'cmr10',
        'axes.prop_cycle': line_cycle,
        'axes.unicode_minus': False, # unicode_minus does not work in Computer Modern font (cm, cmr10)
        'axes.formatter.limits': [-2, 3], # Forces exponent notation below 10**(-2) and above 10**3
        'axes.grid': True,
        'axes.labelpad': 2,
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'axes.linewidth': 0.5,
        'axes.facecolor': '#f8f8f8',
        'figure.dpi': 150.0,
        'figure.subplot.bottom': 0.10,
        'figure.subplot.hspace': 0.38,
        'figure.subplot.left': 0.08,
        'figure.subplot.right': 0.92,
        'figure.subplot.top': 0.78,
        'figure.subplot.wspace': 0.22,
        'grid.alpha': 1.0,
        'grid.color': '#fff',
        'ytick.direction': 'in',
        'xtick.direction': 'in',
        'xtick.labelsize': 9,
        'xtick.major.size': 0.0,
        'ytick.labelsize': 9,
        'ytick.major.size': 0.0,
        'axes.formatter.use_mathtext': True,
        'figure.figsize': [11, 8.5],
        'legend.fontsize': 10,
        'legend.frameon': False,
        'legend.borderpad': 0, 
        'legend.labelspacing': 0,
        'legend.handlelength': 2.5,
        'lines.dash_joinstyle': 'round',
        'lines.dash_capstyle': 'butt',
        'lines.linewidth': 1.0,
    })

TITLEPOS = (0.5, 0.9)
SUBTITLEPOS = (0.5, 0.88)
ROWS = 2
COLS = 3