# Standard Packages
import sys; sys.path.insert(0, '../')
from dataclasses import dataclass

# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt

# Local Packages
import modules.options
import modules.datahelper as datahelper
from plotting.modules.styles import single as plotlayout
from plotting.modules.colors import mmm as plotcolors


@dataclass
class PlotSettings:
    '''
    Settings to control behavior of the plot, set directly by the user.
    '''

    allow_title_runid: bool = True
    allow_title_time: bool = True
    title_override: str = ''
    ylabel_override: str = ''
    xlabel_override: str = ''
    xpad: int = 0
    ypad: int = 0


class PlotData:
    '''
    Makes a PlotData object from user input parameters

    The names of the x-variable and y-variable must match member definitions
    in the InputVariables class.

    Init Parameters:
    * runid (str): The runid of the CDF
    * time (float): The time value to plot
    * yname (str): The name of the y-variable to plot
    * xname (str): The name of the x-variable to plot (Optional)
    * runname (str): A string to replace the runid that shows in plot legends or titles (Optional)
    * legend_override (str): A string to completely replace the legend label of a y-variable (Optional)
    * use_cdf_vars (bool): Uses uncalculated CDF variables instead of calculated MMM variables (Optional)
    '''

    def __init__(self, runid, time, yname, xname='rho', runname='', legend_override='', use_cdf_vars=False):
        options = modules.options.Options(runid=runid, input_time=time)
        mmm_vars, cdf_vars, __ = datahelper.initialize_variables(options)
        plot_vars = mmm_vars if not use_cdf_vars else cdf_vars
        xvar = getattr(plot_vars, xname)
        yvar = getattr(plot_vars, yname)

        self.xvals: np.ndarray = xvar.values[:, plot_vars.options.time_idx]
        self.yvals: np.ndarray = yvar.values[:, plot_vars.options.time_idx]
        self.xsymbol: str = xvar.label
        self.ysymbol: str = yvar.label
        self.xunits: str = xvar.units_label
        self.yunits: str = yvar.units_label
        self.xname: str = xvar.name
        self.yname: str = yvar.name
        self.time: float = plot_vars.options.time_str
        self.runid: str = runid
        self.runname: str = runname
        self.legend_override: str = legend_override

    def get_legend_label(self, legend_attrs):
        '''
        Gets the y-variable label for the plot legend

        If the legend override is set, then that value is used for the legend
        label.  Otherwise, either the y-variable symbol, the runname (if set) or
        runid, and the time value may be added to the legend.  Each of these
        attributes are only added to the legend if multiple y-variables have
        different values for a given attribute.

        Parameters:
        * all_data (list[PlotData]): List of PlotData objects
        * attr_name (str): The name of a LegendAttributes member

        Returns:
        * (str): The y-variable label for the legend
        '''

        # TODO: Show ysymbol vs xsymbol when xsymbols are different
        if self.legend_override:
            return self.legend_override

        legend_items = []
        if legend_attrs.show_ysymbol:
            legend_items.append(self.ysymbol)
        if legend_attrs.show_runid or legend_attrs.show_runname:
            legend_items.append(self.runname or self.runid)
        if legend_attrs.show_time:
            legend_items.append(f'{self.time}s')

        return' '.join(legend_items)


class AllPlotData:

    def __init__(self, *args):
        self.data = [arg for arg in args if isinstance(arg, PlotData)]


class LegendAttributes:
    '''
    Stores a bool for if the corresponding variable attribute should appear in the legend.
    '''

    def __init__(self, ysymbol, runid, runname, time, override):
        self.show_ysymbol: bool = ysymbol
        self.show_runid: bool = runid
        self.show_runname: bool = runname
        self.show_time: bool = time
        self.show_override: bool = override

    def show_legend(self):
        '''Returns (bool): True if the legend should be shown'''

        show_legend = False
        attributes = self.get_attributes()
        for a in attributes:
            if getattr(self, a):
                show_legend = True
                break
        return show_legend

    def get_attributes(self):
        '''Returns (list[str]): all boolean legend attributes'''
        return [a for a in dir(self) if isinstance(getattr(self, a), bool)]


def legend_include_attr(all_pdata, attr_name):
    '''
    Determines if a LegendAttribute should be included in the legend

    The rule for adding an attribute to the legend is if multiple plotted
    variables have different values for a checked attribute.  For example, if
    var A and var B both have different runid's, then their runid's are added
    to the legend.  If there is only one variable to plot, then no attributes
    will be added to the legend.

    Parameters:
    * all_data (list[PlotData]): List of PlotData objects
    * attr_name (str): The name of a LegendAttributes member

    Returns:
    * (bool) True if the attribute should be added to the legend
    '''

    attrs = set()
    for pdata in all_pdata:
        attrs.add(getattr(pdata, attr_name))

    return len(attrs) > 1


def get_plot_limits(psettings, all_pdata):
    '''
    Gets the limits of the plot, adjusted by padding parameters set in psettings

    The axes limits are adjusted by the percentage given for xpad and ypad in
    the psettings object, in each direction for either axis.  For example, if
    ypad = 1, then the limits of the yaxis are both increased and decreased
    by 1%.  Note that MatPlotLib uses default padding of something like 5% on
    each axis, so setting padding values smaller than this value will produce
    tighter plots than what would be created if the padding wasn't adjusted.

    Parameters:
    * psettings (PlotSettings): Plot settings object
    * all_data (list[PlotData]): List of PlotData objects
    '''

    xmin = ymin = float("inf")
    xmax = ymax = -float("inf")
    for pdata in all_pdata:
        xmin = min(xmin, pdata.xvals.min())
        xmax = max(xmin, pdata.xvals.max())
        ymin = min(ymin, pdata.yvals.min())
        ymax = max(ymax, pdata.yvals.max())

    xoffset = (xmax - xmin) * psettings.xpad / 100
    yoffset = (ymax - ymin) * psettings.ypad / 100

    return (xmin - xoffset, xmax + xoffset), (ymin - yoffset, ymax + yoffset)


def get_plot_title(psettings, all_pdata, legend_attrs):
    '''
    Gets the title for the plot

    If the title override is set, then that is used for the title of the plot.
    Otherwise, the name of the first y-variable defined in the all_data list
    is used as the title of the plot.  Additionally, details are added to the
    title in parenthesis if the title details switches are enabled.
    Specifically, if all plotted variables share the same runid, then the
    runid is added to the title.  Similarly, if all plotted variables share
    the same time value, then the time value is added to the title as well.

    Parameters:
    * psettings (PlotSettings): Plot settings object
    * all_data (list[PlotData]): List of PlotData objects
    * legend_attrs (LegendAttributes): Legend data object

    Returns:
    * (str): The title for the plot
    '''

    base_title = psettings.title_override or all_pdata[0].yname

    title_details = ''
    if psettings.allow_title_runid or psettings.allow_title_time:
        title_details_list = []
        if psettings.allow_title_runid and not (legend_attrs.show_runid or legend_attrs.show_runname):  # all lines have same runid or runname
            title_details_list.append(all_pdata[0].runname or all_pdata[0].runid)
        if psettings.allow_title_time and not legend_attrs.show_time:  # all lines have same time
            title_details_list.append(f'{all_pdata[0].time}s')
        if title_details_list:
            title_details = f' ({", ".join(title_details_list)})'

    return f'{base_title}{title_details}'


def get_plot_ylabel(psettings, all_pdata, legend_attrs):
    '''
    Gets the yaxis label for the plot

    If the ylabel override is set, then that is used for the ylabel of the
    plot.  Otherwise, the unique units of each y-variable are added to the
    ylabel (the same units aren't repeated).  Additionally, if the y-variable
    symbols do not appear in the legend, then these unique symbols are also
    added to the ylabel.

    Parameters:
    * psettings (PlotSettings): Plot settings object
    * all_data (list[PlotData]): List of PlotData objects
    * legend_attrs (LegendAttributes): Legend data object

    Returns:
    * (str): The ylabel for the plot
    '''

    if psettings.ylabel_override:
        return psettings.ylabel_override

    ylabels = []  # Not using a set to preserve order
    for pd in all_pdata:
        if not legend_attrs.show_ysymbol:  # ysymbol is not in the legend
            ystr = f'{pd.ysymbol} {pd.yunits}'
            if ystr not in ylabels:
                ylabels.append(ystr)
        elif pd.yunits not in ylabels:  # ysymbol is in the legend
            ylabels.append(pd.yunits)

    return ', '.join(ylabels)


def get_plot_xlabel(psettings, all_pdata, legend_attrs):
    '''
    Gets the xaxis label for the plot

    If the xlabel override is set, then that is used for the xlabel of the
    plot.  Otherwise, the unique symbols and units of each x-variable are
    added to the ylabel (the same units aren't repeated).

    Parameters:
    * psettings (PlotSettings): Plot settings object
    * all_data (list[PlotData]): List of PlotData objects
    * legend_attrs (LegendAttributes): Legend data object

    Returns:
    * (str): The xlabel for the plot
    '''

    if psettings.xlabel_override:
        return psettings.xlabel_override

    xlabels = []  # Not using a set to preserve order
    for pd in all_pdata:
        xstr = f'{pd.xsymbol} {pd.xunits}'
        if xstr not in xlabels:
            xlabels.append(xstr)

    return ', '.join(xlabels)


def main(psettings, all_pdata):
    '''
    Create a plot using CDF data

    Parameters:
    * psettings (PlotSettings): Plot settings object
    * all_data (list[PlotData]): List of PlotData objects
    '''

    plotlayout.init()
    plotcolors.init()
    ax = plt.gca()

    legend_attrs = LegendAttributes(
        ysymbol=legend_include_attr(all_pdata, 'ysymbol'),
        runid=legend_include_attr(all_pdata, 'runid'),
        runname=legend_include_attr(all_pdata, 'runname'),
        time=legend_include_attr(all_pdata, 'time'),
        override=legend_include_attr(all_pdata, 'legend_override')
    )

    for pdata in all_pdata:
        ax.plot(pdata.xvals, pdata.yvals, label=pdata.get_legend_label(legend_attrs))

    xlims, ylims = get_plot_limits(psettings, all_pdata)

    ax.set(
        title=get_plot_title(psettings, all_pdata, legend_attrs),
        xlabel=get_plot_xlabel(psettings, all_pdata, legend_attrs),
        ylabel=get_plot_ylabel(psettings, all_pdata, legend_attrs),
        xlim=xlims,
        ylim=ylims,
    )

    if legend_attrs.show_legend():
        ax.legend()

    plt.show()


# Run this file directly to make a simple plot of variable profiles
if __name__ == '__main__':

    psettings = PlotSettings(
        allow_title_runid=True,
        allow_title_time=True,
        title_override='',
        ylabel_override='',
        xlabel_override='',
        ypad=1,
        xpad=0,
    )

    all_pdata = [
        PlotData(runid='129041A10', yname='gne', xname='rho', time=0.50, runname=''),
        PlotData(runid='120982A09', yname='gne', xname='rho', time=0.50, runname=''),
        # PlotData(runid='120968A02', yname='te', xname='rho', time=0.50, runname=r''),
        # PlotData(runid='138536A01', yname='te', xname='rho', time=0.50, runname=r''),
    ]

    # all_pdata = [
    #     PlotData(runid='129041A10', yname='te', xname='ti', time=0.50, runname=''),
    #     PlotData(runid='129041A10', yname='tau', xname='rho', time=0.50, runname=''),
    # ]

    main(psettings, all_pdata)
