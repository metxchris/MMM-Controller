# 3rd Party Packages
from matplotlib.pyplot import rcParams

# Local Packages
import plotting.modules.plotstyles as ps


def init(style):
    if style is ps.StyleType.Axes.WHITE:
        rcParams.update({
            'axes.grid': False,
            'axes.facecolor': '#fff',
            'ytick.direction': 'in',
            'xtick.direction': 'in',
            'xtick.major.size': 3.0,
            'ytick.major.size': 3.0,
            'xtick.major.width': 0.5,
            'ytick.major.width': 0.5,
        })

    elif style is ps.StyleType.Axes.GRAY:
        rcParams.update({
            'axes.grid': True,
            'axes.facecolor': '#f8f8f8',
            'grid.alpha': 1.0,
            'grid.color': '#fff',
            'xtick.major.size': 0.0,
            'ytick.major.size': 0.0,
        })
