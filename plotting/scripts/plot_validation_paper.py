#!/usr/bin/python3

"""General non-automated plotting for the ETGM paper

Discharges listed below

EAST:
* 80208: T04, R02 xo, q
* 85122: T04, W57, q
* 85124: T02, W02
* 85126: T02, W05, q
* 85610: T01, W02
* 90328: T01, W23

KSTR:
* 15334: T03, P01
* 16295: T10, P03
* 16296: T10, P02 xx
* 16297: T01, P02 xx
* 16299: T10, P02
* 16325: T10, P01
* 16901: T01, P03
* 16949: T02, P01
* 18399: T05, P01
* 18400: T01, P01
* 18402: T01, P01
* 18404: T05, P01
* 18476: T02, P01
* 18492: T01, P02 x?
* 18495: T01, P02
* 18602: T01, P03 xo

DIII-D:
* 118341: T54, V01, T55
* 144449: T54, T52, T53
* 150840: T03, T02, B10
* 153283: T51, T52, T50 xxo

JET:
* Need specific discharges



"""


# Standard Packages
import sys; sys.path.insert(0, '../'), sys.path.insert(0, '../../')
import copy

# 3rd Party Packages
import matplotlib.pyplot as plt

# Local Packages
import modules.options
import modules.utils as utils
from plotting.modules.plotstyles import PlotStyles, StyleType
from plotting.plot_variables import AllPlotData, PlotDataCsv, PlotDataCdf, main


def make_plot(*args, time=None, title=None, **kwargs):
    """Calls plotting.plot_variables.main() to make a plot"""
    plots = [PlotDataCdf(time=time, **arg) for arg in args]
    all_data = copy.deepcopy(base_data)
    all_data.set(*plots, title_override=fr'{title}, $t={time}$s', **kwargs)
    main(all_data, savefig, savedata)


if __name__ == '__main__':
    """Run this module directly to plot variable data stored in CDF or CSV files"""

    utils.init_logging()

    # Initialize visual styles for the plot
    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.AIP,
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
        # 'text.usetex': True,
    })

    # Define settings for the plot
    base_data = AllPlotData(
        replace_offset_text=0,
        allow_title_name=0,
        allow_title_runid=0,
        allow_title_time=0,
        allow_title_factor=True,
        allow_title_rho=True,
        invert_y_axis=False,
        invert_x_axis=False,
        nomralize_y_axis=False,
        nomralize_x_axis=False,
        title_override=r'',
        ylabel_override=r'',
        xlabel_override=r'',
    )

    savefig = 0
    savedata = 1

    # make_plot(
    #     {'runid': '85126T02', 'yname': 'xke', 'legend': 'XKEMMM07'},
    #     {'runid': '85126T02', 'yname': 'conde', 'legend': 'CONDE'},
    #     {'runid': '85126T02', 'yname': 'condepr', 'legend': 'CONDEPR'},
    #     time=2.175,
    #     title='85126:T02',
    #     ymin=0,
    #     ymax=4,
    #     xmax=0.91,
    # )
    # make_plot(
    #     {'runid': '85126T02', 'yname': 'xki', 'legend': 'xki'},
    #     {'runid': '85126T02', 'yname': 'xkimmm07', 'legend': 'XKIMMM07'},
    #     {'runid': '85126T02', 'yname': 'condiwnc', 'legend': 'CONDIWNC'},
    #     {'runid': '85126T02', 'yname': 'condipr', 'legend': 'CONDIPR'},
    #     time=2.175,
    #     title='85126:T02',
    #     ymin=0,
    #     ymax=4,
    #     xmax=0.91,
    # )
    # make_plot(
    #     {'runid': '85610T01', 'yname': 'xke', 'legend': 'XKEMMM07'},
    #     {'runid': '85610T01', 'yname': 'conde', 'legend': 'CONDE'},
    #     {'runid': '85610T01', 'yname': 'condepr', 'legend': 'CONDEPR'},
    #     time=2.118,
    #     title='85610:T01',
    #     ymin=0,
    #     ymax=7,
    #     xmax=0.91,
    # )
    # make_plot(
    #     {'runid': '85610T01', 'yname': 'xki', 'legend': 'XKIMMM07'},
    #     {'runid': '85610T01', 'yname': 'condi', 'legend': 'CONDI'},
    #     {'runid': '85610T01', 'yname': 'condipr', 'legend': 'CONDIPR'},
    #     time=2.118,
    #     title='85610:T01',
    #     ymin=0,
    #     ymax=7,
    #     xmax=0.91,
    # )


    """EAST PLOTS"""

    """85126"""
    # make_plot(
    #     {'runid': '85126T02', 'yname': 'te', 'legend': 'Prediction'},
    #     {'runid': '85126W05', 'yname': 'te', 'legend': 'Experiment'},
    #     title='85126:T02,W05', time=2.175,
    # )
    # make_plot(
    #     {'runid': '85126T02', 'yname': 'ti', 'legend': 'Prediction'},
    #     # {'runid': '85126W05', 'yname': 'ti', 'legend': 'Experiment'},
    #     title='85126:T02,W05', time=2.175,
    # )
    # make_plot(
    #     {'runid': '85126T02', 'yname': 'q', 'legend': 'Prediction'},
    #     {'runid': '85126W05', 'yname': 'q', 'legend': 'Analysis'},
    #     title='85126:T02,W05', time=3.5,
    # )
    # make_plot(
    #     {'runid': '85126T02', 'yname': 'xke'},
    #     {'runid': '85126T02', 'yname': 'xki'},
    #     title='85126:T02', time=2.175,
    #     ymin=0, ymax=4, xmax=0.91,
    # )


    """85610"""
    # make_plot(
    #     {'runid': '85610T01', 'yname': 'te', 'legend': 'Prediction'},
    #     {'runid': '85610W02', 'yname': 'te', 'legend': 'Experiment'},
    #     time=2.118,
    #     title='85610:T01,W02',
    # )
    # make_plot(
    #     {'runid': '85610T01', 'yname': 'ti', 'legend': 'Prediction'},
    #     {'runid': '85610W02', 'yname': 'ti', 'legend': 'Experiment'},
    #     time=2.118,
    #     title='85610:T01,W02',
    # )
    # make_plot(
    #     {'runid': '85610T01', 'yname': 'q', 'legend': 'Prediction'},
    #     {'runid': '85610W02', 'yname': 'q', 'legend': 'Analysis'},
    #     time=3.02,
    #     title='85610:T01,W02',
    # )
    # make_plot(
    #     {'runid': '85610T01', 'yname': 'xke'},
    #     {'runid': '85610T01', 'yname': 'xki'},
    #     title='85610:T01', time=2.118,
    #     ymin=0, ymax=7, xmax=0.91,
    # )


    """85122"""
    # make_plot(
    #     {'runid': '85122T04', 'yname': 'te', 'legend': 'Prediction'},
    #     {'runid': '85122W57', 'yname': 'te', 'legend': 'Experiment'},
    #     time=2.17,
    #     title='85122:T04,W57',
    # )
    # make_plot(
    #     {'runid': '85122T04', 'yname': 'ti', 'legend': 'Prediction'},
    #     {'runid': '85122W57', 'yname': 'ti', 'legend': 'Experiment'},
    #     time=2.17,
    #     title='85122:T04,W57',
    # )
    # make_plot(
    #     {'runid': '85122T04', 'yname': 'q', 'legend': 'Prediction'},
    #     {'runid': '85122W57', 'yname': 'q', 'legend': 'Analysis'},
    #     time=4.38,
    #     title='85122:T04,W57',
    # )
    # make_plot(
    #     {'runid': '85122T04', 'yname': 'xke'},
    #     {'runid': '85122T04', 'yname': 'xki'},
    #     time=2.17,
    #     title='85122:T04',
    #     ymin=0, ymax=4, xmax=0.91,
    # )


    """85124"""
    # make_plot(
    #     {'runid': '85124T02', 'yname': 'te', 'legend': 'Prediction'},
    #     {'runid': '85124W02', 'yname': 'te', 'legend': 'Experiment'},
    #     time=2.076,
    #     title='85124:T02,W02',
    # )
    # make_plot(
    #     {'runid': '85124T02', 'yname': 'ti', 'legend': 'Prediction'},
    #     {'runid': '85124W02', 'yname': 'ti', 'legend': 'Experiment'},
    #     time=2.076,
    #     title='85124:T02,W02',
    # )
    # make_plot(
    #     {'runid': '85124T02', 'yname': 'q', 'legend': 'Prediction'},
    #     {'runid': '85124W02', 'yname': 'q', 'legend': 'Analysis'},
    #     time=3.5,
    #     title='85124:T02,W02',
    # )
    # make_plot(
    #     {'runid': '85124T02', 'yname': 'xke'},
    #     {'runid': '85124T02', 'yname': 'xki'},
    #     time=2.076,
    #     title='85124:T02',
    #     ymin=0, ymax=4, xmax=0.91,
    # )


    """90328"""
    # make_plot(
    #     {'runid': '90328T01', 'yname': 'te', 'legend': 'Prediction'},
    #     {'runid': '90328W23', 'yname': 'te', 'legend': 'Experiment'},
    #     time=1.997,
    #     title='90328:T01,W23',
    # )
    # make_plot(
    #     {'runid': '90328T01', 'yname': 'ti', 'legend': 'Prediction'},
    #     {'runid': '90328W23', 'yname': 'ti', 'legend': 'Experiment'},
    #     time=1.997,
    #     title='90328:T01,W23',
    # )
    # make_plot(
    #     {'runid': '90328T01', 'yname': 'q', 'legend': 'Prediction'},
    #     {'runid': '90328W23', 'yname': 'q', 'legend': 'Analysis'},
    #     time=3.02,
    #     title='90328:T01,W23',
    # )
    make_plot(
        {'runid': '90328T01', 'yname': 'xke'},
        {'runid': '90328T01', 'yname': 'xki'},
        time=1.997,
        title='90328:T01',
        ymin=0, ymax=6, xmax=0.91,
    )


    """80208: CDF MISSING"""
    # make_plot(
    #     {'runid': '80208T04', 'yname': 'te', 'legend': 'Prediction'},
    #     {'runid': '80208R02', 'yname': 'te', 'legend': 'Experiment'},
    #     time=4,
    #     title='80208:T04,R02',
    # )
    # make_plot(
    #     {'runid': '80208T04', 'yname': 'te', 'legend': 'Prediction'},
    #     {'runid': '80208R02', 'yname': 'te', 'legend': 'Experiment'},
    #     time=3.5,
    #     title='80208:T04,R02',
    # )
    # make_plot(
    #     {'runid': '80208T04', 'yname': 'ti', 'legend': 'Prediction'},
    #     {'runid': '80208R02', 'yname': 'ti', 'legend': 'Experiment'},
    #     time=3.5,
    #     title='80208:T04,R02',
    # )
    # make_plot(
    #     {'runid': '80208T04', 'yname': 'ti', 'legend': 'Prediction'},
    #     {'runid': '80208R02', 'yname': 'ti', 'legend': 'Experiment'},
    #     time=4,
    #     title='80208:T04,R02',
    # )
    # make_plot(
    #     {'runid': '80208T04', 'yname': 'q', 'legend': 'Prediction'},
    #     {'runid': '80208R02', 'yname': 'q', 'legend': 'Analysis'},
    #     time=3.0,
    #     title='80208:T04,R02',
    # )
    make_plot(
        {'runid': '80208T04', 'yname': 'xke'},
        {'runid': '80208T04', 'yname': 'xki'},
        time=3.5,
        title='80208:T04',
        ymin=0,
        ymax=5,
        xmax=0.91,
    )



    """KSTR PLOTS"""
