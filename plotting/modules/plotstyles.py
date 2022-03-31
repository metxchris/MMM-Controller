"""Loads styles used for plotting with MatPlotLib

Each style that is loaded is located in the styles subfolder.  See Enum
docstrings below for more information about the different styles that are
loaded.

Dev Note: We are not using custom .mplstyle sheets because the axes.prop_cycle
needs to all be specified on one line, and the prop_cycle's defined in
lines.py would be much harder to maintain if all on one line
"""

# Standard Packages
from enum import Enum

# Local Packages
import plotting.modules.styles.axes
import plotting.modules.styles.layout
import plotting.modules.styles.lines


class StyleType:
    '''
    Defines style type Enums for plotting

    Enums:
    * Axes: Specifies the axes (background, grid) properties of a plot
    * Lines: Specifies the line properties of a plot (colors, dash types, thickness)
    * Layout: Specifies the general layout of a plot
    '''

    class Axes(Enum):
        '''
        Specifies the axes (background, grid) properties of a plot

        Members:
        * WHITE: White background, no grid, with tickmarks
        * GRAY: Gray background, white grid, no tickmarks
        '''

        NONE = 0
        WHITE = 1
        WHITEGRID = 2
        GRAY = 3

    class Lines(Enum):
        '''
        Specifies the line properties of a plot (colors, dash types, thickness)

        Members:
        * MMM: Default theme for MMM Explorer
        * FTE: FiveThirtyEight.com styled lines
        * RHO_MMM: Default theme for MMM Explorer with support for rho plots
        '''

        NONE = 0
        MMM = 1
        FTE = 2
        MAGMA = 3
        RHO_MMM = 4
        RHO_MAGMA = 5

    class Layout(Enum):
        '''
        Specifies the general layout of a plot

        Members:
        * GRID3X2: Layout for a grid of plots with 3 columns and 2 rows
        * SINGLE1: Layout for a single plot (small)
        * SINGLE2: Layout for a single plot (medium)
        * SINGLE3: Layout for a single plot (large)
        '''

        NONE = 0
        GRID3X2 = 1
        SINGLE1 = 2
        SINGLE1B = 3
        SINGLE2 = 4
        SINGLE3 = 5


class PlotStyles:
    '''
    Styles to control the visual elements of the plot

    Style values all correspond to Plot Enums defined in enums.py.  See
    docstrings on each enum type for more information about the different
    values that can be specified for each parameter here.

    Parameters:
    * axes (AxesType): The style of the plot axes
    * layout (LayoutType): The general layout of the plot
    * lines (LineType): The style of the plot lines
    '''

    def __init__(self, axes, lines, layout):
        # Initialize styles
        plotting.modules.styles.axes.init(axes)
        plotting.modules.styles.lines.init(lines)
        plotting.modules.styles.layout.init(layout)

        # Members
        self.dimensions = plotting.modules.styles.layout.Dimensions
