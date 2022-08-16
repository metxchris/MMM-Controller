# Standard Packages
import sys; sys.path.insert(0, '../../')

import matplotlib.pyplot as plt
from matplotlib.pyplot import rcParams
import numpy as np

from plotting.modules.plotstyles import PlotStyles, StyleType

# Initialize visual styles for the plot
PlotStyles(
    axes=StyleType.Axes.WHITE,
    lines=StyleType.Lines.RHO_MMM,
    # layout=StyleType.Layout.SINGLE1B,
)

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
    'patch.linewidth': 0.5,  # legend/bar plot frame line width
    'savefig.format': 'pdf',
    'xtick.top': False,
    'xtick.bottom': False,
    'ytick.right': True,
    'axes.formatter.limits': [-2, 2],  # Forces exponent notation below 1e-1 and above 1e2
    'axes.labelpad': 2,
    'axes.labelsize': 8,
    'axes.titlesize': 8,
    'axes.titlepad': 4,
    'figure.figsize': [5, 1.5],
    'figure.subplot.bottom': 0.16,
    'figure.subplot.hspace': 0.38,
    'figure.subplot.left': 0.1,
    'figure.subplot.right': 0.96,
    'figure.subplot.top': 0.84,
    'figure.subplot.wspace': 0.22,
    'legend.fontsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'ytick.minor.left': True,
    'ytick.minor.pad': 0,
    'ytick.minor.right': True,
    'ytick.minor.size': 1.25,
    'ytick.minor.visible': True,
    'ytick.minor.width': 0.5,
})
 
#data
#x-axis
discharges = ['ITB', r'High-$\beta_{\rm\, N}$',r'High-$q$',r'High-$\beta_{\rm\, P}$']
#y-axis
te_rms = [8.54, 4.825, 11.45, 7]
ti_rms = [14.8, 8.3, 15.13333333, 6.5]
te_off = [5.66, -0.5, -6.1, -7]
ti_off = [14.3, -0.45, 12.13333333, 4.583]
 
#bar chart properties
x = np.arange(len(discharges))
width = 0.3
spacing = 1.05
 
# draw grouped bar chart
fig, ax = plt.subplots()
ax.axhline(0, color='#ccc', linewidth=0.5, dashes=(7,2), zorder=-10)
bar1 = ax.bar(x - spacing*width/2, te_off, width, label=r'$T_{\rm e}$',
        edgecolor=(0.275, 0.482, 0.730), color=(0.094, 0.353, 0.663), hatch='')
bar2 = ax.bar(x + spacing*width/2, ti_off, width, label=r'$T_{\rm i}$',
        edgecolor=(0.994, 0.344, 0.347), color=(0.933, 0.180, 0.184), hatch='/////')

plt.yticks([-10, 0, 10])
 
# ax.set_ylabel('RMS Deviation %')
ax.set_ylabel('Offset %')
# ax.set_title('18 KSTAR Discharges: RMS Using MMM')
ax.set_title('18 KSTAR Discharges: Offset Using MMM')
ax.set_xticks(x, discharges)
ax.legend()

plt.show()

print(rcParams.keys())