# Standard Packages
from cycler import cycler

# 3rd Party Packages
from matplotlib.pyplot import rcParams

# Local Packages
import plotting.modules.plotstyles as ps


def init(style):
    if style is ps.StyleType.Lines.MMM:
        line_cycle = cycler(
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
                (1.0, 0.9),  # Dashed line (short dashes)
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

    elif style is ps.StyleType.Lines.MMM_RHO:
        line_cycle = cycler(
            color=[
                (0.094, 0.353, 0.663),  # Blue
                (0.933, 0.180, 0.184),  # Red
                (0.000, 0.549, 0.282),  # Green
                (0.400, 0.173, 0.569),  # Purple
                (0.957, 0.490, 0.137),  # Orange
                (0.984, 0.722, 0.153),  # Yellow
            ],
            # Defining blank lines using dashes produces artifacts when saving
            # the image as a PDF or EPS
            linestyle=[
                '-',
                '',
                '-.',
                '',
                ':',
                '',
            ],
            marker=[
                '',
                'o',
                '',
                'o',
                '',
                'o',
            ],
            linewidth=[
                1.7,
                0,
                1.6,
                0,
                1.5,
                0,
            ],
            alpha=[
                1.0,
                0.8,
                1.0,
                0.8,
                1.0,
                0.8,
            ]
        )

        rcParams.update({
            'lines.markeredgewidth': 0,
            'lines.markersize': 4.25,
            'lines.dashdot_pattern': [6.5, 1.1, 1.25, 1.1],
            'lines.dashed_pattern': [4.5, 1.5],
            'lines.dotted_pattern': [3, 0.8, 1, 0.5, 1, 0.8],
        })

    elif style is ps.StyleType.Lines.FTE:
        line_cycle = cycler(
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

    rcParams.update({'axes.prop_cycle': line_cycle})
