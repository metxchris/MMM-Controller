"""Loads styles used for plotting with MatPlotLib

Each style that is loaded is located in the styles subfolder.  See Enum
docstrings below for more information about the different styles that are
loaded.
"""

# Standard Packages
from enum import Enum

# Local Packages
import plotting.modules.styles.plotaxes as plotaxes
import plotting.modules.styles.plotlayout as plotlayout
import plotting.modules.styles.plotlines as plotlines


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
        GRAY = 2

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
        RHO_MMM = 3

    class Layout(Enum):
        '''
        Specifies the general layout of a plot

        Members:
        * SINGLE: Layout for a single plot
        * GRID3X2: Layout for a grid of plots with 3 columns and 2 rows
        * SINGLE_SMALL: Layout for a small single plot
        '''

        NONE = 0
        SINGLE = 1
        GRID3X2 = 2
        SINGLESMALL = 3


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
        self.axes: StyleType.Axes = axes
        self.lines: StyleType.Lines = lines
        self.layout: StyleType.Layout = layout

        self.dimensions = plotlayout.Dimensions  # Reference, No values until init_styles is called
        self.init_styles()

    def init_styles(self):
        plotaxes.init(self.axes)
        plotlines.init(self.lines)
        plotlayout.init(self.layout)
