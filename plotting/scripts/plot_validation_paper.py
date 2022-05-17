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
from plotting.plot_variables import AllPlotData, PlotDataCsv, PlotDataCdf, PlotData, main


def make_plot(*args, time=None, timeplot=False, title=None, extradata=None, **kwargs):
    """Calls plotting.plot_variables.main() to make a plot"""
    print(f'Making plot: {title}')
    plots = [PlotDataCdf(zval=time, timeplot=timeplot, **arg) for arg in args]
    if extradata:
        for d in extradata:
            plots.append(d)
    all_data = copy.deepcopy(base_data)
    title_override = f'{title}, ' + (fr'$t={time}$s' if not timeplot else fr'$\rho={time}$')
    all_data.set(*plots, title_override=title_override, **kwargs)

    main(all_data, savefig, savedata)


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

    savefig = 1
    savedata = 1

    plot_east = 0
    plot_kstr = 1
    plot_d3d = 0
    plot_jet = 0

    
    # quit()


    """EAST PLOTS"""
    if plot_east:
        """85126, 85610 Currents"""
        make_plot(
            extradata=(
                csv_to_plotdata('85126LHC.txt', '85126', 4, 'curdlh', 'rho', yunits='AMPS/CM**2', legend='High Density'),
                csv_to_plotdata('85610LHC.txt', '85610', 4, 'curdlh', 'rho', yunits='AMPS/CM**2', legend='Low Density'),
            ),
            title='85126, 85610', time=4.0, ymax=1.4,
        )
        make_plot(
            extradata=(
                csv_to_plotdata('85126LHX.txt', '85126', 4, 'curlh', 'rho', yunits='AMPS', legend='High Density'),
                csv_to_plotdata('85610LHX.txt', '85610', 4, 'curlh', 'rho', yunits='AMPS', legend='Low Density'),
            ),
            title='85126, 85610', time=4.0, ymax=0.4,
        )
        make_plot(
            extradata=(
                csv_to_plotdata('85126LHT.txt', '85126', 0.8, 'curlh', 'time', yunits='AMPS', timeplot=True, legend='High Density'),
                csv_to_plotdata('85610LHT.txt', '85610', 0.8, 'curlh', 'time', yunits='AMPS', timeplot=True, legend='Low Density'),
            ),
            title='85126, 85610', time=0.8, ymax=0.5, xmax=6, xmin=1, timeplot=True,
        )
        make_plot(
            extradata=(
                csv_to_plotdata('85126OHC.txt', '85126', 4, 'curdoh', 'rho', yunits='AMPS/CM**2', legend='High Density'),
                csv_to_plotdata('85610OHC.txt', '85610', 4, 'curdoh', 'rho', yunits='AMPS/CM**2', legend='Low Density'),
            ),
            title='85126, 85610', time=4.0, ymax=2.5,
        )
        make_plot(
            extradata=(
                csv_to_plotdata('85126OHX.txt', '85126', 4, 'curoh', 'rho', yunits='AMPS', legend='High Density'),
                csv_to_plotdata('85610OHX.txt', '85610', 4, 'curoh', 'rho', yunits='AMPS', legend='Low Density'),
            ),
            title='85126, 85610', time=4.0, ymax=0.4,
        )
        make_plot(
            extradata=(
                csv_to_plotdata('85126oHT.txt', '85126', 0.8, 'curoh', 'time', yunits='AMPS', timeplot=True, legend='High Density'),
                csv_to_plotdata('85610oHT.txt', '85610', 0.8, 'curoh', 'time', yunits='AMPS', timeplot=True, legend='Low Density'),
            ),
            title='85126, 85610', time=0.8, ymax=0.6, xmax=6, xmin=1, timeplot=True,
        )

        """85126"""
        make_plot(
            {'runid': '85126T02', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '85126W05', 'yname': 'te', 'legend': 'Experiment'},
            title='85126', time=2.175,
        )
        make_plot(
            {'runid': '85126T02', 'yname': 'ti', 'legend': 'Prediction'},
            # {'runid': '85126W05', 'yname': 'ti', 'legend': 'Experiment'},
            title='85126', time=2.175,
        )
        make_plot(
            {'runid': '85126T02', 'yname': 'q', 'legend': 'Prediction'},
            {'runid': '85126W05', 'yname': 'q', 'legend': 'Analysis'},
            extradata=(csv_to_plotdata('85126q.csv', '85126', 3.5, 'q', 'rho', legend='Experiment'),),
            title='85126', time=3.5, ymax=7, ymin=0,
        )
        make_plot(
            {'runid': '85126T02', 'yname': 'xke'},
            {'runid': '85126T02', 'yname': 'xki'},
            title='85126', time=2.175, ymin=0, ymax=4, xmax=0.91,
        )


        """85610"""
        make_plot(
            {'runid': '85610T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '85610W02', 'yname': 'te', 'legend': 'Experiment'},
            time=2.118, title='85610',
        )
        make_plot(
            {'runid': '85610T01', 'yname': 'ti', 'legend': 'Prediction'},
            # {'runid': '85610W02', 'yname': 'ti', 'legend': 'Experiment'},
            time=2.118, title='85610',
        )
        make_plot(
            {'runid': '85610T01', 'yname': 'q', 'legend': 'Prediction'},
            {'runid': '85610W02', 'yname': 'q', 'legend': 'Analysis'},
            extradata=(csv_to_plotdata('85610q.csv', '85610', 3.02, 'q', 'rho', legend='Experiment'),),
            time=3.02, title='85610', ymax=7, ymin=0,
        )
        make_plot(
            {'runid': '85610T01', 'yname': 'xke'},
            {'runid': '85610T01', 'yname': 'xki'},
            title='85610', time=2.118, ymin=0, ymax=7, xmax=0.91,
        )


        """85122"""
        make_plot(
            {'runid': '85122T04', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '85122W57', 'yname': 'te', 'legend': 'Experiment'},
            time=2.17, title='85122',
        )
        make_plot(
            {'runid': '85122T04', 'yname': 'ti', 'legend': 'Prediction'},
            # {'runid': '85122W57', 'yname': 'ti', 'legend': 'Experiment'},
            time=2.17,
            title='85122',
        )
        make_plot(
            {'runid': '85122T04', 'yname': 'q', 'legend': 'Prediction'},
            {'runid': '85122W57', 'yname': 'q', 'legend': 'Analysis'},
            extradata=(csv_to_plotdata('85122q.csv', '85122', 4.38, 'q', 'rho', legend='Experiment'),),
            time=4.38, title='85122', ymax=7, ymin=0,
        )
        make_plot(
            {'runid': '85122T04', 'yname': 'xke'},
            {'runid': '85122T04', 'yname': 'xki'},
            time=2.17, title='85122', ymin=0, ymax=4, xmax=0.91,
        )


        """85124"""
        make_plot(
            {'runid': '85124T02', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '85124W02', 'yname': 'te', 'legend': 'Experiment'},
            time=2.076, title='85124',
        )
        make_plot(
            {'runid': '85124T02', 'yname': 'ti', 'legend': 'Prediction'},
            # {'runid': '85124W02', 'yname': 'ti', 'legend': 'Experiment'},
            time=2.076, title='85124',
        )
        make_plot(
            {'runid': '85124T02', 'yname': 'q', 'legend': 'Prediction'},
            {'runid': '85124W02', 'yname': 'q', 'legend': 'Analysis'},
            extradata=(csv_to_plotdata('85124q.csv', '85124', 3.5, 'q', 'rho', legend='Experiment'),),
            time=3.5, title='85124', ymax=7, ymin=0,
        )
        make_plot(
            {'runid': '85124T02', 'yname': 'xke'},
            {'runid': '85124T02', 'yname': 'xki'},
            time=2.076, title='85124', ymin=0, ymax=4, xmax=0.91,
        )


        """90328"""
        make_plot(
            {'runid': '90328T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '90328W23', 'yname': 'te', 'legend': 'Experiment'},
            time=1.997, title='90328',
        )
        make_plot(
            {'runid': '90328T01', 'yname': 'ti', 'legend': 'Prediction'},
            # {'runid': '90328W23', 'yname': 'ti', 'legend': 'Experiment'},
            time=1.997, title='90328',
        )
        make_plot(
            {'runid': '90328T01', 'yname': 'q', 'legend': 'Prediction'},
            {'runid': '90328W23', 'yname': 'q', 'legend': 'Analysis'},
            extradata=(csv_to_plotdata('90328q.csv', '90328', 3.02, 'q', 'rho', legend='Experiment'),),
            time=3.02, title='90328', ymax=7, ymin=0,
        )
        make_plot(
            {'runid': '90328T01', 'yname': 'xke'},
            {'runid': '90328T01', 'yname': 'xki'},
            time=1.997, title='90328', ymin=0, ymax=6, xmax=0.91,
        )

        """80208"""
        # make_plot(
        #     {'runid': '80208T04', 'yname': 'te', 'legend': 'Prediction'},
        #     {'runid': '80208R02', 'yname': 'te', 'legend': 'Experiment'},
        #     time=4, title='80208',
        # )
        # make_plot(
        #     {'runid': '80208T04', 'yname': 'te', 'legend': 'Prediction'},
        #     {'runid': '80208R02', 'yname': 'te', 'legend': 'Experiment'},
        #     time=3.5, title='80208',
        # )
        # make_plot(
        #     {'runid': '80208T04', 'yname': 'ti', 'legend': 'Prediction'},
        #     # {'runid': '80208R02', 'yname': 'ti', 'legend': 'Experiment'},
        #     time=3.5, title='80208',
        # )
        # make_plot(
        #     {'runid': '80208T04', 'yname': 'ti', 'legend': 'Prediction'},
        #     # {'runid': '80208R02', 'yname': 'ti', 'legend': 'Experiment'},
        #     time=4, title='80208',
        # )
        # make_plot(
        #     {'runid': '80208T04', 'yname': 'q', 'legend': 'Prediction'},
        #     {'runid': '80208R02', 'yname': 'q', 'legend': 'Analysis'},
        #     extradata=csv_to_plotdata('80208q', '80208', 3.0, 'q', 'rho', legend='Experiment'),
        #     time=3.0, title='80208',
        # )
        make_plot(
            {'runid': '80208T04', 'yname': 'xke'},
            {'runid': '80208T04', 'yname': 'xki'},
            time=3.5, title='80208', ymin=0, ymax=5, xmax=0.91,
        )


    """KSTR PLOTS"""
    if plot_kstr:
        """Beta Plots"""
        make_plot(
            extradata=(
                csv_to_plotdata('16295BTN.dat', '16295', 4, 'betant', 'time', yunits='', legend=''),
                csv_to_plotdata('16295BTH.dat', '16295', 4, 'betanh', 'time', yunits='', legend=''),
            ),
            title='16295', time=0.5, timeplot=True,
        )
        make_plot(
            extradata=(
                csv_to_plotdata('16296BTN.dat', '16296', 4, 'betant', 'time', yunits='', legend=''),
                csv_to_plotdata('16296BTH.dat', '16296', 4, 'betanh', 'time', yunits='', legend=''),
            ),
            title='16296', time=0.5, timeplot=True,
        )
        make_plot(
            extradata=(
                csv_to_plotdata('16297BTN.dat', '16297', 4, 'betant', 'time', yunits='', legend=''),
                csv_to_plotdata('16297BTH.dat', '16297', 4, 'betanh', 'time', yunits='', legend=''),
            ),
            title='16297', time=0.5, timeplot=True,
        )
        make_plot(
            extradata=(
                csv_to_plotdata('16299BTN.dat', '16299', 4, 'betant', 'time', yunits='', legend=''),
                csv_to_plotdata('16299BTH.dat', '16299', 4, 'betanh', 'time', yunits='', legend=''),
            ),
            title='16299', time=0.5, timeplot=True, ymax=5,
        )
        make_plot(
            extradata=(
                csv_to_plotdata('16901BTN.dat', '16901', 4, 'betant', 'time', yunits='', legend=''),
                csv_to_plotdata('16901BTH.dat', '16901', 4, 'betanh', 'time', yunits='', legend=''),
            ),
            title='16901', time=0.5, timeplot=True, ymax=5,
        )
        """16901"""
        make_plot(
            {'runid': '16901T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '16901P03', 'yname': 'te', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='16901', xmin=1, xmax=9, ymin=0.8,
        )
        make_plot(
            {'runid': '16901T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '16901P03', 'yname': 'ti', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='16901', xmin=1, xmax=9, ymin=0.8,
        )
        make_plot(
            {'runid': '16901T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '16901P03', 'yname': 'te', 'legend': 'Experiment'},
            time=5.0, title='16901',
        )
        make_plot(
            {'runid': '16901T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '16901P03', 'yname': 'ti', 'legend': 'Experiment'},
            time=5.0, title='16901',
        )
        make_plot(
            {'runid': '16901T01', 'yname': 'xke', 'legend': 'Total'},
            {'runid': '16901T01', 'yname': 'xkemmm07', 'legend': 'MMM'},
            time=5.0, title='16901', xmax=0.8,
        )
        make_plot(
            {'runid': '16901T01', 'yname': 'xki', 'legend': 'Total'},
            {'runid': '16901T01', 'yname': 'xkimmm07', 'legend': 'MMM'},
            time=5.0, title='16901', xmax=0.8, ymax=5,
        )

        """16295, 16296, 16297"""
        make_plot(
            {'runid': '16295T10', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '16295P03', 'yname': 'te', 'legend': 'Experiment'},
            # time=0.5, timeplot=True, title='16295:T10,P03', xmin=2, ymax=1.6,
            time=0.5, timeplot=True, title='16295', xmin=2, ymax=1.6,
        )
        make_plot(
            {'runid': '16295T10', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '16295P03', 'yname': 'ti', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='16295', xmin=2,
        )
        make_plot(
            {'runid': '16296T10', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '16296P02', 'yname': 'te', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='16296', xmin=0.8, xmax=3, ymax=3,
        )
        make_plot(
            {'runid': '16296T10', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '16296P02', 'yname': 'ti', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='16296', xmin=0.8, xmax=3,
        )
        make_plot(
            {'runid': '16297T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '16297P02', 'yname': 'te', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='16297', xmin=0.8,
        )
        make_plot(
            {'runid': '16297T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '16297P02', 'yname': 'ti', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='16297', xmin=0.8,
        )

        """18399, 18400, 16325"""
        make_plot(
            {'runid': '18399T05', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18399P01', 'yname': 'te', 'legend': 'Experiment'},
            time=10.0, title='18399',
        )
        make_plot(
            {'runid': '18399T05', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18399P01', 'yname': 'ti', 'legend': 'Experiment'},
            time=10.0, title='18399',
        )
        make_plot(
            {'runid': '18400T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18400P01', 'yname': 'te', 'legend': 'Experiment'},
            time=11.5, title='18400',
        )
        make_plot(
            {'runid': '18400T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18400P01', 'yname': 'ti', 'legend': 'Experiment'},
            time=11.5, title='18400',
        )
        make_plot(
            {'runid': '16325T10', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '16325P01', 'yname': 'te', 'legend': 'Experiment'},
            time=10.0, title='16325',
        )
        make_plot(
            {'runid': '16325T10', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '16325P01', 'yname': 'ti', 'legend': 'Experiment'},
            time=10.0, title='16325',
        )

        """18404, 16949"""
        make_plot(
            {'runid': '18404T05', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18404P01', 'yname': 'te', 'legend': 'Experiment'},
            time=12.0, title='18404',
        )
        make_plot(
            {'runid': '18404T05', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18404P01', 'yname': 'ti', 'legend': 'Experiment'},
            time=12.0, title='18404',
        )
        make_plot(
            {'runid': '16949T02', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '16949P01', 'yname': 'te', 'legend': 'Experiment'},
            time=9.0, title='16949',
        )
        make_plot(
            {'runid': '16949T02', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '16949P01', 'yname': 'ti', 'legend': 'Experiment'},
            time=9.0, title='16949',
        )

        """18402"""
        make_plot(
            {'runid': '18402T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18402P01', 'yname': 'te', 'legend': 'Experiment'},
            time=11.5, title='18402',
        )
        make_plot(
            {'runid': '18402T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18402P01', 'yname': 'ti', 'legend': 'Experiment'},
            time=11.5, title='18402',
        )
        make_plot(
            {'runid': '18402T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18402P01', 'yname': 'te', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='18402', xmin=1, xmax=12,
        )
        make_plot(
            {'runid': '18402T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18402P01', 'yname': 'ti', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='18402', xmin=1, xmax=12,
        )
        make_plot(
            {'runid': '18402T01', 'yname': 'q',},
            {'runid': '18402T01', 'yname': 'xke',},
            {'runid': '18402T01', 'yname': 'xki',},
            time=10.0, title='18402',
        )
        make_plot(
            {'runid': '18402T01', 'yname': 'q',},
            {'runid': '18402T01', 'yname': 'xke',},
            {'runid': '18402T01', 'yname': 'xki',},
            time=0.5, timeplot=True, title='18402', xmin=1, xmax=12,
        )

        """18476"""
        make_plot(
            {'runid': '18476T02', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18476P01', 'yname': 'te', 'legend': 'Experiment'},
            time=7.0, title='18476',
        )
        make_plot(
            {'runid': '18476T02', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18476P01', 'yname': 'ti', 'legend': 'Experiment'},
            time=7.0, title='18476',
        )
        make_plot(
            {'runid': '18476T02', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18476P01', 'yname': 'te', 'legend': 'Experiment'},
            timeplot=True, time=0.4, title='18476', xmin=1, xmax=7,
        )
        make_plot(
            {'runid': '18476T02', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18476P01', 'yname': 'ti', 'legend': 'Experiment'},
            timeplot=True, time=0.4, title='18476', xmin=1, xmax=7,
        )
        make_plot(
            {'runid': '18476T02', 'yname': 'omega'},
            time=7.0, title='18476',
        )
        make_plot(
            {'runid': '18476T02', 'yname': 'xke', 'legend': 'Total'},
            {'runid': '18476T02', 'yname': 'xkemmm07', 'legend': 'MMM'},
            time=7.0, title='18476',
            ymin=0, xmax=0.8,
        )
        make_plot(
            {'runid': '18476T02', 'yname': 'xki', 'legend': 'Total'},
            {'runid': '18476T02', 'yname': 'xkimmm07', 'legend': 'MMM'},
            time=7.0, title='18476',
            ymin=0, ymax=6, xmax=0.8,
        )

        """18477"""
        make_plot(
            {'runid': '18477T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18477P02', 'yname': 'te', 'legend': 'Experiment'},
            time=6.0, title='18477',
        )
        make_plot(
            {'runid': '18477T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18477P02', 'yname': 'ti', 'legend': 'Experiment'},
            time=6.0, title='18477',
        )

        """18492"""
        make_plot(
            {'runid': '18492T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18492P02', 'yname': 'te', 'legend': 'Experiment'},
            time=6.7, title='18492',
        )
        make_plot(
            {'runid': '18492T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18492P02', 'yname': 'ti', 'legend': 'Experiment'},
            time=6.7, title='18492',
        )

        """18495"""
        make_plot(
            {'runid': '18495T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18495P02', 'yname': 'te', 'legend': 'Experiment'},
            time=6.1, title='18495',
        )
        make_plot(
            {'runid': '18495T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18495P02', 'yname': 'ti', 'legend': 'Experiment'},
            time=6.1, title='18495',
        )

        """18499"""
        make_plot(
            {'runid': '18499T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18499P02', 'yname': 'te', 'legend': 'Experiment'},
            time=6.0, title='18499',
        )
        make_plot(
            {'runid': '18499T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18499P02', 'yname': 'ti', 'legend': 'Experiment'},
            time=6.0, title='18499',
        )

        """18602"""
        make_plot(
            {'runid': '18602T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18602P03', 'yname': 'te', 'legend': 'Experiment'},
            time=6.9, title='18602',
        )
        make_plot(
            {'runid': '18602T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18602P03', 'yname': 'ti', 'legend': 'Experiment'},
            time=6.9, title='18602',
        )
        make_plot(
            {'runid': '18602T01', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '18602P03', 'yname': 'te', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='18602', xmin=4, xmax=8, ymin=1.0,
        )
        make_plot(
            {'runid': '18602T01', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '18602P03', 'yname': 'ti', 'legend': 'Experiment'},
            time=0.5, timeplot=True, title='18602', xmin=4, xmax=8, ymin=0.8,
        )

    """D3D PLOTS"""
    if plot_d3d:
        """118341"""
        make_plot(
            {'runid': '118341T54', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '118341V01', 'yname': 'te', 'legend': 'Experiment'},
            time=5.85, title='118341',
        )
        make_plot(
            {'runid': '118341T54', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '118341V01', 'yname': 'ti', 'legend': 'Experiment'},
            time=5.85, title='118341',
        )
        make_plot(
            {'runid': '118341T54', 'yname': 'ne', 'legend': 'Prediction'},
            {'runid': '118341V01', 'yname': 'ne', 'legend': 'Experiment'},
            time=5.85, title='118341',
        )
        make_plot(
            {'runid': '118341T54', 'yname': 'omega', 'legend': 'Prediction'},
            {'runid': '118341V01', 'yname': 'omega', 'legend': 'Experiment'},
            time=5.85, title='118341',
        )

        """144449"""
        make_plot(
            {'runid': '144449T54', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '144449T52', 'yname': 'te', 'legend': 'Experiment'},
            time=3.0, title='144449',
        )
        make_plot(
            {'runid': '144449T54', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '144449T52', 'yname': 'ti', 'legend': 'Experiment'},
            time=3.0, title='144449',
        )
        make_plot(
            {'runid': '144449T54', 'yname': 'ne', 'legend': 'Prediction'},
            {'runid': '144449T52', 'yname': 'ne', 'legend': 'Experiment'},
            time=3.0, title='144449',
        )

        """153283"""
        make_plot(
            {'runid': '153283T51', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '153283T52', 'yname': 'te', 'legend': 'Experiment'},
            time=3.3, title='153283',
        )
        make_plot(
            {'runid': '153283T51', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '153283T52', 'yname': 'ti', 'legend': 'Experiment'},
            time=3.3, title='153283',
        )
        make_plot(
            {'runid': '153283T51', 'yname': 'ne', 'legend': 'Prediction'},
            {'runid': '153283T52', 'yname': 'ne', 'legend': 'Experiment'},
            time=3.3, title='153283',
        )
        make_plot(
            {'runid': '153283T51', 'yname': 'omega', 'legend': 'Prediction'},
            {'runid': '153283T52', 'yname': 'omega', 'legend': 'Experiment'},
            time=3.3, title='153283',
        )

        """150840"""
        make_plot(
            {'runid': '150840T03', 'yname': 'te', 'legend': 'Prediction'},
            {'runid': '150840T02', 'yname': 'te', 'legend': 'Experiment'},
            time=3.1, title='150840',
        )
        make_plot(
            {'runid': '150840T03', 'yname': 'ti', 'legend': 'Prediction'},
            {'runid': '150840T02', 'yname': 'ti', 'legend': 'Experiment'},
            time=3.1, title='150840',
        )
        make_plot(
            {'runid': '150840T03', 'yname': 'ne', 'legend': 'Prediction'},
            {'runid': '150840T02', 'yname': 'ne', 'legend': 'Experiment'},
            time=3.1, title='150840',
        )
        make_plot(
            {'runid': '150840T03', 'yname': 'omega', 'legend': 'Prediction'},
            {'runid': '150840T02', 'yname': 'omega', 'legend': 'Experiment'},
            time=3.1, title='150840',
        )


    """JET PLOTS"""
    if plot_jet:
        ...
    
