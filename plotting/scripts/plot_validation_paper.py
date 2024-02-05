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
* 16295: T10, P03 xx
* 16296: T10, P02 xx
* 16297: T01, P02 xx
* 16299: T01, P02 xx
* 16325: T10, P01 xx
* 16901: T01, P03
* 16949: T02, P01 xx
* 18399: T05, P01 xx
* 18400: T01, P01 xx
* 18402: T01, P01 xx
* 18404: T05, P01 xx
* 18476: T02, P01 xx
* 18492: T01, P02 xx
* 18495: T01, P02 xx
* 18602: T01, P03 xo

DIII-D:
* 118341: T54, V01, T55 xxx
* 144449: T54, T52, T53 xxx
* 150840: T03, T02, B10 xxo
* 153283: T51, T52, T50 xxx

JET:
* Need specific discharges



"""


# Standard Packages
import sys; sys.path.insert(0, '../'), sys.path.insert(0, '../../')
import copy

# 3rd Party Packages
import matplotlib.pyplot as plt
import numpy as np

# Local Packages
import modules.options
import modules.conversions as conversions
import modules.variables as variables
import modules.utils as utils
from plotting.modules.plotstyles import PlotStyles, StyleType
from plotting.plot_variables import FigData, PlotDataCsv, PlotDataCdf, PlotData, main


def make_plot(*args, time=None, timeplot=False, title=None, extradata=None, **kwargs):
    """Calls plotting.plot_variables.main() to make a plot"""
    print(f'Making plot: {title}')
    plots = [PlotDataCdf(zval=time, timeplot=timeplot, **arg) for arg in args]
    if extradata:
        for d in extradata:
            plots.append(d)
    fig_data = copy.deepcopy(base_data)
    title_override = f'{title}, ' + (fr'$t={time}$s' if not timeplot else fr'$\hat{{\rho}}={time}$')
    fig_data.set(*plots, title_override=title_override, **kwargs)

    main(fig_data)


def csv_to_plotdata(filename, runid, zval, yname, xname, yunits='', xunits='', timeplot=False, runname='', legend=''):
    delimiter = ''
    if '.csv' in filename:
        delimiter = ','

    csvdata = np.genfromtxt(f'data/{filename}', delimiter=delimiter)
    options = modules.options.Options(
        runid=runid, input_time=zval, ignore_exceptions=True
    )

    input_vars = variables.InputVariables(options)
    output_vars = variables.OutputVariables(options)
    # input_controls = controls.InputControls(options)

    if hasattr(input_vars, yname):
        yvar = getattr(input_vars, yname) if hasattr(input_vars, yname) else getattr(output_vars, yname)
    else:
        yvar = variables.Variable(yname)
    if hasattr(input_vars, xname):
        xvar = getattr(input_vars, xname) if hasattr(input_vars, xname) else getattr(output_vars, xname)
    else:
        xvar = variables.Variable(xname)

    yvar.values = csvdata[:, 1]
    xvar.values = csvdata[:, 0]

    if yunits:
        yvar.units = yunits
        conversions.convert_units(yvar)
    if xunits:
        xvar.units = xunits
        conversions.convert_units(xvar)

    zidx = 0

    return PlotData(
        options, runid, yname, xname, yvar, xvar, zidx, zval,
        timeplot=timeplot, runname=runname, legend=legend,
    )


if __name__ == '__main__':
    """Run this module directly to plot variable data stored in CDF or CSV files"""

    utils.init_logging()

    # Initialize visual styles for the plot
    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.AIP2,
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
        # 'text.usetex': True,
    })

    # Define settings for the plot
    base_data = FigData(
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
        savefig=True,
        savedata=True,
    )

    plot_east = 0
    plot_kstr = 0
    plot_d3d = 1
    plot_jet = 0

    
    # quit()


    """EAST PLOTS"""
    if plot_east:
        # """85126, 85610 Currents"""
        # make_plot(
        #     extradata=(
        #         csv_to_plotdata('85126LHC.txt', '85126', 4, 'curdlh', 'rho', yunits='AMPS/CM**2', legend='High Density'),
        #         csv_to_plotdata('85610LHC.txt', '85610', 4, 'curdlh', 'rho', yunits='AMPS/CM**2', legend='Low Density'),
        #     ),
        #     title='EAST 85126, 85610', time=4.0, ymax=1.4, xmin=0, xmax=1,
        # )
        # make_plot(
        #     extradata=(
        #         csv_to_plotdata('85126LHX.txt', '85126', 4, 'curlh', 'rho', yunits='AMPS', legend='High Density'),
        #         csv_to_plotdata('85610LHX.txt', '85610', 4, 'curlh', 'rho', yunits='AMPS', legend='Low Density'),
        #     ),
        #     title='EAST 85126, 85610', time=4.0, ymin=0, ymax=0.5, xmin=0, xmax=1,
        # )
        # make_plot(
        #     extradata=(
        #         csv_to_plotdata('85126LHT.txt', '85126', 1.0, 'curlh', 'time', yunits='AMPS', timeplot=True, legend='High Density'),
        #         csv_to_plotdata('85610LHT.txt', '85610', 1.0, 'curlh', 'time', yunits='AMPS', timeplot=True, legend='Low Density'),
        #     ),
        #     title='EAST 85126, 85610', time=1.0, ymin=0, ymax=0.5, xmax=6, xmin=1, timeplot=True, xticks=np.arange(1,7),
        # )
        # make_plot(
        #     extradata=(
        #         csv_to_plotdata('85126OHC.txt', '85126', 4, 'curdoh', 'rho', yunits='AMPS/CM**2', legend='High Density'),
        #         csv_to_plotdata('85610OHC.txt', '85610', 4, 'curdoh', 'rho', yunits='AMPS/CM**2', legend='Low Density'),
        #     ),
        #     title='EAST 85126, 85610', time=4.0, ymax=2.5, xmin=0, xmax=1,
        # )
        # make_plot(
        #     extradata=(
        #         csv_to_plotdata('85126OHX.txt', '85126', 4, 'curoh', 'rho', yunits='AMPS', legend='High Density'),
        #         csv_to_plotdata('85610OHX.txt', '85610', 4, 'curoh', 'rho', yunits='AMPS', legend='Low Density'),
        #     ),
        #     title='EAST 85126, 85610', time=4.0, ymin=0, ymax=0.5, xmin=0, xmax=1,
        # )
        # make_plot(
        #     extradata=(
        #         csv_to_plotdata('85126oHT.txt', '85126', 1.0, 'curoh', 'time', yunits='AMPS', timeplot=True, legend='High Density'),
        #         csv_to_plotdata('85610oHT.txt', '85610', 1.0, 'curoh', 'time', yunits='AMPS', timeplot=True, legend='Low Density'),
        #     ),
        #     title='EAST 85126, 85610', time=1.0, ymin=0, ymax=0.5, xmax=6, xmin=1, timeplot=True, xticks=np.arange(1,7),
        # )

        """85126"""
        make_plot(
            {'runid': '85126T02', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '85126W05', 'yname': 'te', 'legend': 'Experiment'},
            title='EAST 85126', time=2.175, ymin=0, ymax=1.5,
        )
        make_plot(
            {'runid': '85126T02', 'yname': 'ti', 'legend': 'Prediction'},
            # {'runid': '85126W05', 'yname': 'ti', 'legend': 'Experiment'},
            title='EAST 85126', time=2.175, ymin=0, ymax=1.5,
        )
        make_plot(
            extradata=(
                csv_to_plotdata('85126q.dat', '85126', 3.5, 'q', 'rho', legend='Prediction'),
                csv_to_plotdata('85126qexp.dat', '85126', 3.5, 'q', 'rho', legend='Experiment'),
                csv_to_plotdata('85126qa.dat', '85126', 3.5, 'q', 'rho', legend='Analysis'),
            ),
            title='EAST 85126', time=3.5, ymax=7, ymin=0, xmin=0, xmax=1, yticks=np.arange(0, 8),
        )
        make_plot(
            {'runid': '85126T02', 'yname': 'xke'},
            {'runid': '85126T02', 'yname': 'xki'},
            title='EAST 85126', time=2.175, ymin=0, ymax=4, xmax=0.91,
        )


        """85610"""
        make_plot(
            {'runid': '85610T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '85610W02', 'yname': 'te', 'legend': 'Experiment'},
            time=2.118, title='EAST 85610', ymax=1.25,
        )
        make_plot(
            {'runid': '85610T01', 'yname': 'ti', 'legend': 'Prediction'},
            # {'runid': '85610W02', 'yname': 'ti', 'legend': 'Experiment'},
            time=2.118, title='EAST 85610', ymax=1.25,
        )
        make_plot(
            extradata=(
                csv_to_plotdata('85610q.dat', '85610', 3.02, 'q', 'rho', legend='Prediction'),
                csv_to_plotdata('85610qexp.dat', '85610', 3.02, 'q', 'rho', legend='Experiment'),
                csv_to_plotdata('85610qa.dat', '85610', 3.02, 'q', 'rho', legend='Analysis'),
            ),
            time=3.02, title='EAST 85610', ymax=7, ymin=0, xmin=0, xmax=1, yticks=np.arange(0, 8),
        )
        make_plot(
            {'runid': '85610T01', 'yname': 'xke'},
            {'runid': '85610T01', 'yname': 'xki'},
            title='EAST 85610', time=2.118, ymin=0, ymax=7, xmax=0.91,
        )


        ##"""85122"""
        ##make_plot(
        ##    {'runid': '85122T04', 'yname': 'te', 'legend': 'Prediction'},
        ##    {'runid': '85122W57', 'yname': 'te', 'legend': 'Experiment'},
        ##    time=2.17, title='EAST 85122',
        ##)
        ##make_plot(
        ##    {'runid': '85122T04', 'yname': 'ti', 'legend': 'Prediction'},
        ##    # {'runid': '85122W57', 'yname': 'ti', 'legend': 'Experiment'},
        ##    time=2.17, title='EAST 85122',
        ##)
        ##make_plot(
        ##    {'runid': '85122T04', 'yname': 'q', 'legend': 'Prediction'},
        ##    {'runid': '85122W57', 'yname': 'q', 'legend': 'Analysis'},
        ##    extradata=(csv_to_plotdata('85122q.csv', '85122', 4.38, 'q', 'rho', legend='Experiment'),),
        ##    time=4.38, title='EAST 85122', ymax=7, ymin=0,
        ##)
        ##make_plot(
        ##    {'runid': '85122T04', 'yname': 'xke'},
        ##    {'runid': '85122T04', 'yname': 'xki'},
        ##    time=2.17, title='EAST 85122', ymin=0, ymax=4, xmax=0.91,
        ##)


        ##"""85124"""
        ##make_plot(
        ##    {'runid': '85124T02', 'yname': 'te', 'legend': 'Prediction'},
        ##    {'runid': '85124W02', 'yname': 'te', 'legend': 'Experiment'},
        ##    time=2.076, title='EAST 85124',
        ##)
        ##make_plot(
        ##    {'runid': '85124T02', 'yname': 'ti', 'legend': 'Prediction'},
        ##    # {'runid': '85124W02', 'yname': 'ti', 'legend': 'Experiment'},
        ##    time=2.076, title='EAST 85124',
        ##)
        ##make_plot(
        ##    {'runid': '85124T02', 'yname': 'q', 'legend': 'Prediction'},
        ##    {'runid': '85124W02', 'yname': 'q', 'legend': 'Analysis'},
        ##    extradata=(csv_to_plotdata('85124q.csv', '85124', 3.5, 'q', 'rho', legend='Experiment'),),
        ##    time=3.5, title='EAST 85124', ymax=7, ymin=0,
        ##)
        ##make_plot(
        ##    {'runid': '85124T02', 'yname': 'xke'},
        ##    {'runid': '85124T02', 'yname': 'xki'},
        ##    time=2.076, title='EAST 85124', ymin=0, ymax=4, xmax=0.91,
        ##)


        ##"""90328"""
        ##make_plot(
        ##    {'runid': '90328T01', 'yname': 'te', 'legend': 'Prediction'},
        ##    {'runid': '90328W23', 'yname': 'te', 'legend': 'Experiment'},
        ##    time=1.997, title='EAST 90328',
        ##)
        ##make_plot(
        ##    {'runid': '90328T01', 'yname': 'ti', 'legend': 'Prediction'},
        ##    # {'runid': '90328W23', 'yname': 'ti', 'legend': 'Experiment'},
        ##    time=1.997, title='EAST 90328',
        ##)
        ##make_plot(
        ##    {'runid': '90328T01', 'yname': 'q', 'legend': 'Prediction'},
        ##    {'runid': '90328W23', 'yname': 'q', 'legend': 'Analysis'},
        ##    extradata=(csv_to_plotdata('90328q.csv', '90328', 3.02, 'q', 'rho', legend='Experiment'),),
        ##    time=3.02, title='EAST 90328', ymax=7, ymin=0,
        ##)
        ##make_plot(
        ##    {'runid': '90328T01', 'yname': 'xke'},
        ##    {'runid': '90328T01', 'yname': 'xki'},
        ##    time=1.997, title='EAST 90328', ymin=0, ymax=6, xmax=0.91,
        ##)

        """80208"""
        make_plot(
            extradata=(
                csv_to_plotdata('80208TE4.dat', '80208', 4.0, 'te', 'rho', legend='Prediction'),
                csv_to_plotdata('80208TEexp4.dat', '80208', 4.0, 'te', 'rho', legend='Experiment'),
            ),
            time=4, title='EAST 80208', xmin=0, xmax=1, ymin=0, ymax=2, savename_append='4',
        )
        make_plot(
            extradata=(
                csv_to_plotdata('80208TE.dat', '80208', 3.5, 'te', 'rho', legend='Prediction'),
                csv_to_plotdata('80208TEexp.dat', '80208', 3.5, 'te', 'rho', legend='Experiment'),
            ),
            time=3.5, title='EAST 80208', xmin=0, xmax=1, ymin=0, ymax=2,
        )
        make_plot(
            extradata=(
                csv_to_plotdata('80208TI4.dat', '80208', 4.0, 'ti', 'rho', legend='Prediction'),
            ),
            time=4, title='EAST 80208', xmin=0, xmax=1, ymin=0, ymax=2, savename_append='4',
        )
        make_plot(
            extradata=(
                csv_to_plotdata('80208TI.dat', '80208', 3.5, 'ti', 'rho', legend='Prediction'),
            ),
            time=3.5, title='EAST 80208', xmin=0, xmax=1, ymin=0, ymax=2,
        )
        make_plot(
            extradata=(
                csv_to_plotdata('80208q.dat', '80208', 3.0, 'q', 'rho', legend='Prediction'),
                csv_to_plotdata('80208qexp.dat', '80208', 3.0, 'q', 'rho', legend='Experiment'),
                csv_to_plotdata('80208qa.dat', '80208', 3.0, 'q', 'rho', legend='Analysis'),
            ),
            time=3.0, title='EAST 80208', xmin=0, xmax=1, ymin=0,
            yticks=np.arange(0, 9),
        )
        make_plot(
            {'runid': '80208T04', 'yname': 'xke'},
            {'runid': '80208T04', 'yname': 'xki'},
            time=3.5, title='EAST 80208', ymin=0, ymax=5, xmax=0.91,
        )

    """KSTR PLOTS"""
    if plot_kstr:
        """Beta Plots"""
        # make_plot(
        #     extradata=(
        #         csv_to_plotdata('16295BTN.dat', '16295', 4, 'betant', 'time', yunits='', legend=''),
        #         csv_to_plotdata('16295BTH.dat', '16295', 4, 'betanh', 'time', yunits='', legend=''),
        #     ),
        #     title='KSTAR 16295', time=0.5, timeplot=True, xmin=2, ymin=0,
        # )
        # make_plot(
        #     extradata=(
        #         csv_to_plotdata('16296BTN.dat', '16296', 4, 'betant', 'time', yunits='', legend=''),
        #         csv_to_plotdata('16296BTH.dat', '16296', 4, 'betanh', 'time', yunits='', legend=''),
        #     ),
        #     title='KSTAR 16296', time=0.5, timeplot=True, xmin=1, ymin=0,
        # )
        # make_plot(
        #     extradata=(
        #         csv_to_plotdata('16297BTN.dat', '16297', 4, 'betant', 'time', yunits='', legend=''),
        #         csv_to_plotdata('16297BTH.dat', '16297', 4, 'betanh', 'time', yunits='', legend=''),
        #     ),
        #     title='KSTAR 16297', time=0.5, timeplot=True, xmin=1, ymin=0,
        # )
        # make_plot(
        #     extradata=(
        #         csv_to_plotdata('16299BTN.dat', '16299', 4, 'betant', 'time', yunits='', legend=''),
        #         csv_to_plotdata('16299BTH.dat', '16299', 4, 'betanh', 'time', yunits='', legend=''),
        #     ),
        #     title='KSTAR 16299', time=0.5, timeplot=True, ymax=5, ymin=0, xmin=1, 
        # )
        # make_plot(
        #     extradata=(
        #         csv_to_plotdata('16901BTN.dat', '16901', 4, 'betant', 'time', yunits='', legend=''),
        #         csv_to_plotdata('16901BTH.dat', '16901', 4, 'betanh', 'time', yunits='', legend=''),
        #     ),
        #     title='KSTAR 16901', time=0.5, timeplot=True, ymax=5, ymin=0, xmin=1, 
        # )


        """16901"""
        make_plot(
            {'runid': '16901T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '16901P03', 'yname': 'te', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='KSTAR 16901', xmin=1, xmax=9, ymin=0, ymax=2.5,
        )
        PlotStyles(lines=StyleType.Lines.COMPARE_MMM)
        make_plot(
            {'runid': '16901T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '16901P03', 'yname': 'ti', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='KSTAR 16901', xmin=1, xmax=9, ymin=0, ymax=2.5,
        )
        PlotStyles(lines=StyleType.Lines.RHO_MMM)
        make_plot(
            {'runid': '16901T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '16901P03', 'yname': 'te', 'legend': 'Experiment'},
            time=5.0, title='KSTAR 16901', ymin=0, ymax=3,
        )
        make_plot(
            {'runid': '16901T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '16901P03', 'yname': 'ti', 'legend': 'Experiment'},
            time=5.0, title='KSTAR 16901', ymin=0, ymax=3,
        )
        make_plot(
            {'runid': '16901T01', 'yname': 'xke', 'legend': 'Total'},
            {'runid': '16901T01', 'yname': 'xkemmm07', 'legend': 'MMM'},
            time=5.0, title='KSTAR 16901', xmax=0.8, ymin=0, ymax=5, 
        )
        make_plot(
            {'runid': '16901T01', 'yname': 'xki', 'legend': 'Total'},
            {'runid': '16901T01', 'yname': 'xkimmm07', 'legend': 'MMM'},
            time=5.0, title='KSTAR 16901', xmax=0.8, ymax=5, ymin=0,
        )

        """16295, 16296, 16297"""
        make_plot(
            {'runid': '16295T10', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '16295P03', 'yname': 'te', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='KSTAR 16295', xmin=2, ymax=2.5, ymin=0,
        )
        make_plot(
            {'runid': '16296T10', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '16296P02', 'yname': 'te', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='KSTAR 16296', xmin=0.8, xmax=3, ymax=4, ymin=0,
        )
        make_plot(
            {'runid': '16297T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '16297P02', 'yname': 'te', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='KSTAR 16297', xmin=0.8, ymin=0, ymax=2.5,
        )
        PlotStyles(lines=StyleType.Lines.COMPARE_MMM)
        make_plot(
            {'runid': '16295T10', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '16295P03', 'yname': 'ti', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='KSTAR 16295', xmin=2, ymin=0, ymax=2.5,
        )
        make_plot(
            {'runid': '16296T10', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '16296P02', 'yname': 'ti', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='KSTAR 16296', xmin=0.8, xmax=3, ymin=0, ymax=4,
        )
        make_plot(
            {'runid': '16297T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '16297P02', 'yname': 'ti', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='KSTAR 16297', xmin=0.8, ymin=0, ymax=2.5,
        )
        PlotStyles(lines=StyleType.Lines.RHO_MMM)

        ## """18399, 18400, 16325"""
        ## make_plot(
        ##     {'runid': '18399T05', 'yname': 'te', 'legend': 'Prediction'},
        ##     {'runid': '18399P01', 'yname': 'te', 'legend': 'Experiment'},
        ##     time=10.0, title='KSTAR 18399', ymin=0, ymax=3.5,
        ## )
        ## make_plot(
        ##     {'runid': '18399T05', 'yname': 'ti', 'legend': 'Prediction'},
        ##     {'runid': '18399P01', 'yname': 'ti', 'legend': 'Experiment'},
        ##     time=10.0, title='KSTAR 18399', ymin=0, ymax=3.5,
        ## )
        ## make_plot(
        ##     {'runid': '18400T01', 'yname': 'te', 'legend': 'Prediction'},
        ##     {'runid': '18400P01', 'yname': 'te', 'legend': 'Experiment'},
        ##     time=11.5, title='KSTAR 18400', ymin=0, ymax=3,
        ## )
        ## make_plot(
        ##     {'runid': '18400T01', 'yname': 'ti', 'legend': 'Prediction'},
        ##     {'runid': '18400P01', 'yname': 'ti', 'legend': 'Experiment'},
        ##     time=11.5, title='KSTAR 18400', ymin=0, ymax=3,
        ## )
        ## make_plot(
        ##     {'runid': '16325T10', 'yname': 'te', 'legend': 'Prediction'},
        ##     {'runid': '16325P01', 'yname': 'te', 'legend': 'Experiment'},
        ##     time=10.0, title='KSTAR 16325', ymin=0,
        ## )
        ## make_plot(
        ##     {'runid': '16325T10', 'yname': 'ti', 'legend': 'Prediction'},
        ##     {'runid': '16325P01', 'yname': 'ti', 'legend': 'Experiment'},
        ##     time=10.0, title='KSTAR 16325', ymin=0,
        ## )

        ##"""18404, 16949"""
        ##make_plot(
        ##    {'runid': '18404T05', 'yname': 'te', 'legend': 'Prediction'},
        ##    {'runid': '18404P01', 'yname': 'te', 'legend': 'Experiment'},
        ##    time=12.0, title='KSTAR 18404', ymin=0,
        ##)
        ##make_plot(
        ##    {'runid': '18404T05', 'yname': 'ti', 'legend': 'Prediction'},
        ##    {'runid': '18404P01', 'yname': 'ti', 'legend': 'Experiment'},
        ##    time=12.0, title='KSTAR 18404', ymin=0,
        ##)
        ##make_plot(
        ##    {'runid': '16949T02', 'yname': 'te', 'legend': 'Prediction'},
        ##    {'runid': '16949P01', 'yname': 'te', 'legend': 'Experiment'},
        ##    time=9.0, title='KSTAR 16949', ymin=0,
        ##)
        ##make_plot(
        ##    {'runid': '16949T02', 'yname': 'ti', 'legend': 'Prediction'},
        ##    {'runid': '16949P01', 'yname': 'ti', 'legend': 'Experiment'},
        ##    time=9.0, title='KSTAR 16949', ymin=0,
        ##)

        """18402"""
        make_plot(
            {'runid': '18402T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18402P01', 'yname': 'te', 'legend': 'Experiment'},
            time=11.5, title='KSTAR 18402', ymin=0, ymax=3.5,
        )
        make_plot(
            {'runid': '18402T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18402P01', 'yname': 'ti', 'legend': 'Experiment'},
            time=11.5, title='KSTAR 18402', ymin=0, ymax=3.5,
        )
        PlotStyles(lines=StyleType.Lines.COMPARE_MMM)
        make_plot(
            {'runid': '18402T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18402P01', 'yname': 'te', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='KSTAR 18402', xmin=1, xmax=12, ymin=0, ymax=2,
        )
        make_plot(
            {'runid': '18402T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18402P01', 'yname': 'ti', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='KSTAR 18402', xmin=1, xmax=12, ymin=0, ymax=2,
        )
        PlotStyles(lines=StyleType.Lines.RHO_MMM)
        make_plot(
            {'runid': '18402T01', 'yname': 'q', 'runname': ' '},
            {'runid': '18402T01', 'yname': 'xke', 'runname': r'(m$^2$/s)'},
            {'runid': '18402T01', 'yname': 'xki', 'runname': r'(m$^2$/s)'},
            time=10.0, title='KSTAR 18402', ymin=0, ylabel_override=' ',
        )
        make_plot(
            {'runid': '18402T01', 'yname': 'q', 'runname': ' '},
            {'runid': '18402T01', 'yname': 'xke', 'runname': r'(m$^2$/s)'},
            {'runid': '18402T01', 'yname': 'xki', 'runname': r'(m$^2$/s)'},
            time=0.5, timeplot=True, title='KSTAR 18402', ylabel_override=' ',
            xmin=1.1, xmax=12, ymin=0, ymax=4.5,
        )

        """18476"""
        make_plot(
            {'runid': '18476T02', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18476P01', 'yname': 'te', 'legend': 'Experiment'},
            time=7.0, title='KSTAR 18476', ymin=0, ymax=3.5,
        )
        make_plot(
            {'runid': '18476T02', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18476P01', 'yname': 'ti', 'legend': 'Experiment'},
            time=7.0, title='KSTAR 18476', ymin=0, ymax=3.5,
        )
        make_plot(
            {'runid': '18476T02', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18476P01', 'yname': 'te', 'legend': 'Experiment'},
            timeplot=True, time=0.4, title='KSTAR 18476', xmin=1, xmax=7, ymin=0, xticks=np.arange(1,8),
        )
        make_plot(
            {'runid': '18476T02', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18476P01', 'yname': 'ti', 'legend': 'Experiment'},
            timeplot=True, time=0.4, title='KSTAR 18476', xmin=1, xmax=7, ymin=0, xticks=np.arange(1,8),
        )
        make_plot(
            {'runid': '18476T02', 'yname': 'omega'},
            time=7.0, title='KSTAR 18476', ymin=0, ymax=8e4,
        )
        make_plot(
            {'runid': '18476T02', 'yname': 'xke', 'legend': 'Total'},
            {'runid': '18476T02', 'yname': 'xkemmm07', 'legend': 'MMM'},
            time=7.0, title='KSTAR 18476', ymin=0, xmax=0.8,
        )
        make_plot(
            {'runid': '18476T02', 'yname': 'xki', 'legend': 'Total'},
            {'runid': '18476T02', 'yname': 'xkimmm07', 'legend': 'MMM'},
            time=7.0, title='KSTAR 18476', ymin=0, ymax=6, xmax=0.8,
        )

        """18477"""
        make_plot(
            {'runid': '18477T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18477P02', 'yname': 'te', 'legend': 'Experiment'},
            time=6.0, title='KSTAR 18477', ymin=0,
        )
        make_plot(
            {'runid': '18477T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18477P02', 'yname': 'ti', 'legend': 'Experiment'},
            time=6.0, title='KSTAR 18477', ymin=0,
        )

        """18602"""
        make_plot(
            {'runid': '18602T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18602P03', 'yname': 'te', 'legend': 'Experiment'},
            time=6.9, title='KSTAR 18602', ymin=0, ymax=3,
        )
        make_plot(
            {'runid': '18602T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18602P03', 'yname': 'ti', 'legend': 'Experiment'},
            time=6.9, title='KSTAR 18602', ymin=0, ymax=3,
        )
        make_plot(
            {'runid': '18602T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18602P03', 'yname': 'te', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='KSTAR 18602', xmin=4, xmax=8, ymin=0, ymax=2.5,
        )
        PlotStyles(lines=StyleType.Lines.COMPARE_MMM)
        make_plot(
            {'runid': '18602T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18602P03', 'yname': 'ti', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='KSTAR 18602', xmin=4, xmax=8, ymin=0, ymax=2.5,
        )
        PlotStyles(lines=StyleType.Lines.RHO_MMM)

        # # """18492"""
        # # make_plot(
        # #     {'runid': '18492T01', 'yname': 'te', 'legend': 'Prediction'},
        # #     {'runid': '18492P02', 'yname': 'tepro', 'legend': 'Experiment'},
        # #     time=6.7, title='KSTAR 18492', ymin=0,
        # # )
        # # make_plot(
        # #     {'runid': '18492T01', 'yname': 'ti', 'legend': 'Prediction'},
        # #     {'runid': '18492P02', 'yname': 'ti', 'legend': 'Experiment'},
        # #     time=6.7, title='KSTAR 18492', ymin=0,
        # # )

        # # """18495"""
        # # make_plot(
        # #     {'runid': '18495T01', 'yname': 'te', 'legend': 'Prediction'},
        # #     {'runid': '18495P02', 'yname': 'te', 'legend': 'Experiment'},
        # #     time=6.1, title='KSTAR 18495', ymin=0,
        # # )
        # # make_plot(
        # #     {'runid': '18495T01', 'yname': 'ti', 'legend': 'Prediction'},
        # #     {'runid': '18495P02', 'yname': 'ti', 'legend': 'Experiment'},
        # #     time=6.1, title='KSTAR 18495', ymin=0,
        # # )

        # # """18499"""
        # # make_plot(
        # #     {'runid': '18499T01', 'yname': 'te', 'legend': 'Prediction'},
        # #     {'runid': '18499P02', 'yname': 'te', 'legend': 'Experiment'},
        # #     time=6.0, title='KSTAR 18499', ymin=0,
        # # )
        # # make_plot(
        # #     {'runid': '18499T01', 'yname': 'ti', 'legend': 'Prediction'},
        # #     {'runid': '18499P02', 'yname': 'ti', 'legend': 'Experiment'},
        # #     time=6.0, title='KSTAR 18499', ymin=0,
        # # )

    """D3D PLOTS"""
    if plot_d3d:
        """118341"""
        make_plot(
            {'runid': '118341T55', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '118341V01', 'yname': 'te', 'legend': 'Experiment'},
            time=5.85, title='DIII-D 118341', xmin=0, xmax=1, ymin=0, ymax=5,
        )
        make_plot(
            {'runid': '118341T55', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '118341V01', 'yname': 'ti', 'legend': 'Experiment'},
            time=5.85, title='DIII-D 118341', xmin=0, xmax=1, ymin=0, ymax=5,
        )
        make_plot(
            {'runid': '118341T55', 'yname': 'ne', 'legend': 'Prediction'},
            {'runid': '118341V01', 'yname': 'ne', 'legend': 'Experiment'},
            time=5.85, title=r'$\ \ $DIII-D 118341', xmin=0, xmax=1, ymin=0,
        )
        # make_plot(
        #     {'runid': '118341T55', 'yname': 'omega', 'legend': 'Prediction'},
        #     {'runid': '118341V01', 'yname': 'omegadata', 'legend': 'Experiment'},
        #     time=5.85, title=r'$\ \ $DIII-D 118341', xmin=0, xmax=1, ymin=0,
        # )

        """144449"""
        make_plot(
            {'runid': '144449T54', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '144449T52', 'yname': 'te', 'legend': 'Experiment'},
            time=3.0, title='DIII-D 144449', xmin=0, xmax=1, ymin=0, ymax=5,
        )
        make_plot(
            {'runid': '144449T54', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '144449T52', 'yname': 'ti', 'legend': 'Experiment'},
            time=3.0, title='DIII-D 144449', xmin=0, xmax=1, ymin=0, ymax=5,
        )
        make_plot(
            {'runid': '144449T54', 'yname': 'ne', 'legend': 'Prediction'},
            {'runid': '144449T52', 'yname': 'ne', 'legend': 'Experiment'},
            time=3.0, title=r'$\ \ $DIII-D 144449', xmin=0, xmax=1, ymin=0,
        )

        """153283"""
        make_plot(
            {'runid': '153283T51', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '153283T52', 'yname': 'te', 'legend': 'Experiment'},
            time=3.3, title='DIII-D 153283', xmin=0, xmax=1, ymin=0,
        )
        make_plot(
            {'runid': '153283T51', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '153283T52', 'yname': 'ti', 'legend': 'Experiment'},
            time=3.3, title='DIII-D 153283', xmin=0, xmax=1, ymin=0,
        )
        make_plot(
            {'runid': '153283T51', 'yname': 'ne', 'legend': 'Prediction'},
            {'runid': '153283T52', 'yname': 'ne', 'legend': 'Experiment'},
            time=3.3, title=r'$\ \ $DIII-D 153283', xmin=0, xmax=1, ymin=0,
        )
        # make_plot(
        #     {'runid': '153283T51', 'yname': 'omega', 'legend': 'Prediction'},
        #     {'runid': '153283T52', 'yname': 'omegadata', 'legend': 'Experiment'},
        #     time=3.3, title=r'$\ \ $DIII-D 153283', xmin=0, xmax=1,
        # )

        """150840"""
        make_plot(
            {'runid': '150840T03', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '150840T02', 'yname': 'te', 'legend': 'Experiment'},
            time=3.1, title='DIII-D 150840', xmin=0, xmax=1, ymin=0, ymax=3,
        )
        make_plot(
            {'runid': '150840T03', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '150840T02', 'yname': 'ti', 'legend': 'Experiment'},
            time=3.1, title='DIII-D 150840', xmin=0, xmax=1, ymin=0, ymax=3,
        )
        make_plot(
            {'runid': '150840T03', 'yname': 'ne', 'legend': 'Prediction'},
            {'runid': '150840T02', 'yname': 'ne', 'legend': 'Experiment'},
            time=3.1, title=r'$\ \ $DIII-D 150840', xmin=0, xmax=1, ymin=0,
        )
        # make_plot(
        #     {'runid': '150840T03', 'yname': 'omega', 'legend': 'Prediction'},
        #     {'runid': '150840T02', 'yname': 'omegadata', 'legend': 'Experiment'},
        #     time=3.1, title=r'$\ \ $DIII-D 150840', xmin=0, xmax=1, ymin=0,
        # )

    """JET PLOTS"""
    if plot_jet:
        make_plot(
            extradata=(
                csv_to_plotdata('84599TE.dat', '84599', 11.5, 'te', 'rho', yunits='', legend='Prediction'),
                csv_to_plotdata('84599TEexp.dat', '84599', 11.5, 'te', 'rho', yunits='', legend='Experiment'),
            ),
            title='JET 84599', time=11.5, timeplot=False, ymax=2.5, ymin=0, xmin=0, xmax=1,
        )

        make_plot(
            extradata=(
                csv_to_plotdata('84599TI.dat', '84599', 11.5, 'ti', 'rho', yunits='', legend='Prediction'),
                csv_to_plotdata('84599TIexp.dat', '84599', 11.5, 'ti', 'rho', yunits='', legend='Experiment'),
            ),
            title='JET 84599', time=11.5, timeplot=False, ymax=2.5, ymin=0, xmin=0, xmax=1,
        )

        make_plot(
            extradata=(
                csv_to_plotdata('86911TE.dat', '86911', 12.53, 'te', 'rho', yunits='', legend='Prediction'),
                csv_to_plotdata('86911TEexp.dat', '86911', 12.53, 'te', 'rho', yunits='', legend='Experiment'),
            ),
            title='JET 86911', time=12.53, timeplot=False, ymax=3, ymin=0, xmin=0, xmax=1,
        )

        make_plot(
            extradata=(
                csv_to_plotdata('86911TI.dat', '86911', 12.53, 'ti', 'rho', yunits='', legend='Prediction'),
                csv_to_plotdata('86911TIexp.dat', '86911', 12.53, 'ti', 'rho', yunits='', legend='Experiment'),
            ),
            title='JET 86911', time=12.53, timeplot=False, ymax=3, ymin=0, xmin=0, xmax=1,
        )

        make_plot(
            extradata=(
                csv_to_plotdata('87261TE.dat', '87261', 7.55, 'te', 'rho', yunits='', legend='Prediction'),
                csv_to_plotdata('87261TEexp.dat', '87261', 7.55, 'te', 'rho', yunits='', legend='Experiment'),
            ),
            title='JET 87261', time=7.55, timeplot=False, ymax=4, ymin=0, xmin=0, xmax=1,
        )

        make_plot(
            extradata=(
                csv_to_plotdata('87261TI.dat', '87261', 7.55, 'ti', 'rho', yunits='', legend='Prediction'),
                csv_to_plotdata('87261TIexp.dat', '87261', 7.55, 'ti', 'rho', yunits='', legend='Experiment'),
            ),
            title='JET 87261', time=7.55, timeplot=False, ymax=4, ymin=0, xmin=0, xmax=1,
        )

        make_plot(
            extradata=(
                csv_to_plotdata('87215TE.dat', '87215', 8.5, 'te', 'rho', yunits='', legend='Prediction'),
                csv_to_plotdata('87215TEexp.dat', '87215', 8.5, 'te', 'rho', yunits='', legend='Experiment'),
            ),
            title='JET 87215', time=8.5, timeplot=False, ymax=6, ymin=0, xmin=0, xmax=1,
        )

        make_plot(
            extradata=(
                csv_to_plotdata('87215TI.dat', '87215', 8.5, 'ti', 'rho', yunits='', legend='Prediction'),
                csv_to_plotdata('87215TIexp.dat', '87215', 8.5, 'ti', 'rho', yunits='', legend='Experiment'),
            ),
            title='JET 87215', time=8.5, timeplot=False, ymax=6, ymin=0, xmin=0, xmax=1,
        )


