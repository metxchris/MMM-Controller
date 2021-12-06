# Standard Packages
import sys
sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt

# Local Packages
from main import utils
from main.enums import ShotType
from main.options import Options

def init_variables(rho_value, var_to_plot):
    Options.instance.load_options(runid, scan_num)

    rho_path = utils.get_rho_path(runid, scan_num, Options.instance.var_to_scan)

    scan_range = Options.instance.scan_range

    rho_str = '{:.3f}'.format(rho_value)

    input_file = f'{rho_path}\\Input rho = {rho_str}.csv'
    output_file = f'{rho_path}\\Output rho = {rho_str}.csv'

    input_data = np.genfromtxt(input_file, delimiter=',', dtype=float)
    input_vars = np.array(np.genfromtxt(input_file, delimiter=',', dtype=float, names=True).dtype.names)

    output_data = np.genfromtxt(output_file, delimiter=',', dtype=float)
    output_vars = np.array(np.genfromtxt(output_file, delimiter=',', dtype=float, names=True).dtype.names)

    xvals = input_data[:, input_vars == Options.instance.var_to_scan]
    yvals = output_data[:, output_vars == var_to_plot]

    from plots.styles import standard as ps
    from plots.colors import mmm

    fig = plt.figure(figsize=(3.5,3))
    plt.subplots_adjust(left=0.15, right=0.9, bottom=0.15, top=0.9)

    plot_title = f'{var_to_plot} vs. {Options.instance.var_to_scan}' + r' ($\rho = $' + rho_str + ')'
    plt.plot(xvals, yvals)
    plt.xlabel(Options.instance.var_to_scan)
    plt.ylabel(var_to_plot)
    plt.title(plot_title)
    plt.show()


if __name__ == '__main__':
    runid = '129041A10'
    scan_num = 16
    rho_value = .64
    var_to_plot = 'xteETGM'

    init_variables(rho_value, var_to_plot)
