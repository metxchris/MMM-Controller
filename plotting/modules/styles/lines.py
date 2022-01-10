# Standard Packages
from cycler import cycler

# 3rd Party Packages
from matplotlib.pyplot import rcParams

# Local Packages
import plotting.modules.plotstyles


def init(style):
    Lines = plotting.modules.plotstyles.StyleType.Lines

    if style is Lines.MMM:
        prop_cycle = cycler(
            color=[
                (0.094, 0.353, 0.663),  # Blue
                (0.933, 0.180, 0.184),  # Red
                (0.000, 0.549, 0.282),  # Green
                (0.400, 0.173, 0.569),  # Purple
                (0.957, 0.490, 0.137),  # Orange
                (0.984, 0.722, 0.153),  # Yellow
            ],
            dashes=[
                (1, 0),  # Solid line
                (6.5, 1.1, 1.25, 1.1),  # Dash-dot
                (4.5, 1.5),  # Dashed line
                (3, 0.8, 1, 0.5, 1, 0.8),  # Dash-dot-dot
                (1, 0.25, 0.5, 0.5),  # Dashed line (short dashes)
                (1.5, 0.5, 1.5, 0.5, 1.5, 2.0),  # Dashes as Dot-dot-dot
            ],
            linewidth=[
                1.7,
                1.6,
                1.5,
                1.5,
                1.4,
                1.4,
            ],
        )

    elif style is Lines.RHO_MMM:
        '''
        There's a bug in MatPlotLib when saving PDF and EPS when specifying
        dashes as a blank line using (0, 1), where the legend label of a
        marker contains line artifacts only when the file is saved (it looks
        fine using plt.show()).  Furthermore, if we try and set the linewidth
        to 0 to fix this, the PDF fails to save correctly.

        Consequently, we are setting the lines to white here so that any
        artifacts blend into the background (the artifacts should only be a
        pixel per legend marker), and the corresponding line widths to 0.1.
        '''
        prop_cycle = cycler(
            color=[
                (0.094, 0.353, 0.663),  # Blue
                (1, 1, 1),  # White
                (0.933, 0.180, 0.184),  # Red
                (1, 1, 1),  # White
                (0.000, 0.549, 0.282),  # Green
                (1, 1, 1),  # White
                (0.400, 0.173, 0.569),  # Purple
                (1, 1, 1),  # White
                (0.957, 0.490, 0.137),  # Orange
                (1, 1, 1),  # White
                (0.984, 0.722, 0.153),  # Yellow
                (1, 1, 1),  # White
            ],
            markeredgecolor=[
                (1, 1, 1),  # White
                (0.094, 0.353, 0.663),  # Blue
                (1, 1, 1),  # White
                (0.933, 0.180, 0.184),  # Red
                (1, 1, 1),  # White
                (0.000, 0.549, 0.282),  # Green
                (1, 1, 1),  # White
                (0.400, 0.173, 0.569),  # Purple
                (1, 1, 1),  # White
                (0.957, 0.490, 0.137),  # Orange
                (1, 1, 1),  # White
                (0.984, 0.722, 0.153),  # Yellow
            ],
            dashes=[
                (1, 0),  # Solid line
                (0, 100),  # Blank
                (6.5, 1.1, 1.25, 1.1),  # Dash-dot
                (0, 100),  # Blank
                (4.5, 1.5),  # Dashed line
                (0, 100),  # Blank
                (3, 0.8, 1, 0.5, 1, 0.8),  # Dash-dot-dot
                (0, 100),  # Blank
                (1, 0.25, 0.5, 0.5),  # Dashed line (short dashes)
                (0, 100),  # Blank
                (1.25, 0.5, 1.25, 0.5, 1.25, 1.5),  # Dashes as Dot-dot-dot
                (0, 100),  # Blank
            ],
            marker=[
                '',
                'o',
                '',
                'o',
                '',
                'o',
                '',
                'o',
                '',
                'o',
                '',
                'o',
            ],
            linewidth=[
                1.7,
                0.1,
                1.6,
                0.1,
                1.5,
                0.1,
                1.5,
                0.1,
                1.4,
                0.1,
                1.4,
                0.1,
            ],
        )

        rcParams.update({
            'lines.markeredgecolor': 'auto',
            'lines.markerfacecolor': '#fff',
            'lines.markeredgewidth': 1.25,
            'lines.markersize': 4,
            'lines.dashdot_pattern': [6.5, 1.1, 1.25, 1.1],
            'lines.dashed_pattern': [4.5, 1.5],
            'lines.dotted_pattern': [3, 0.8, 1, 0.5, 1, 0.8],
        })

    elif style is Lines.FTE:
        prop_cycle = cycler(
            color=[
                '#008fd5',
                '#fc4f30',
                '#e5ae38',
                '#6d904f',
                '#8b8b8b',
                '#810f7c',
            ],
            dashes=[
                (1, 0),  # Solid line
                (6.5, 1.5, 1, 1.5),  # Dash-dot
                (4.5, 1.5),  # Dashed line
                (3, 1, 1, 0.5, 1, 1),  # Dash-dot-dot
                (1.5, 1),  # Dashed line (short dashes)
                (1.5, 0.5, 1.5, 0.5, 1.5, 2.0),  # Dashes as Dot-dot-dot
            ],
            linewidth=[
                1.7,
                1.6,
                1.5,
                1.5,
                1.4,
                1.4,
            ],
        )

    rcParams.update({'axes.prop_cycle': prop_cycle})
