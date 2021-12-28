# Standard Packages
import sys; sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt

# Local Packages
import modules.datahelper as datahelper
import modules.options as options
from modules.enums import ShotType
from plotting.modules.styles import single as plotlayout
from plotting.modules.colors import mmm as plotcolors


def simple_plot(x1var, y1var, l1='', x2var=None, y2var=None, l2=''):

    plotlayout.init()
    plotcolors.init()

    opts = options.instance

    t_idx = opts.time_idx

    plt.figure(figsize=(3.5, 3))
    plt.subplots_adjust(left=0.15, right=0.9, bottom=0.15, top=0.9)

    if x1var.values is not None and y1var.values is not None:
        plt.plot(x1var.values[:, t_idx], y1var.values[:, t_idx], label=y1var.label + l1)
    if x2var is not None and y2var is not None:
        plt.plot(x2var.values[:, t_idx], y2var.values[:, t_idx], label=y2var.label + l2)

    xmin = min(x1var.values.min(), x2var.values.min()) if x2var is not None else x1var.values.min()
    xmax = max(x1var.values.max(), x2var.values.max()) if x2var is not None else x1var.values.max()

    plt.xlim(xmin, xmax)
    plt.xlabel(x1var.label)
    plt.ylabel(y1var.units_label)
    plt.legend()
    plt.title(f'{opts.runid}, t={opts.time_str}s')
    plt.show()


# Run this file directly to make a simple plot of variable profiles
if __name__ == '__main__':
    '''
    CDF Options:
    * Uncomment the line you wish to use
    * Edit enums.py to view or add additional ShotTypes
    '''
    cdf_name, shot_type, input_time = '120982A09', ShotType.NSTX, 0.50
    # cdf_name, shot_type, input_time = '129041A10', ShotType.NSTX, 0.5
    # cdf_name, shot_type, input_time = '132017T01', ShotType.D3D, 2.1
    # cdf_name, shot_type, input_time = '141552A01', ShotType.D3D, 2.1

    '''
    Input Options:
    * input_points is the number of points to use when making the MMM input file
    * Set input_points = None to match the number of points used in the CDF
    * Set var_to_scan = var_name (str) to run a scan of the specified variable
    * Set var_to_scan = None to skip the variable scan
    * E.g.: var_to_scan = 'te'
    '''

    # TODO: Add local dataclass to replace the need to set options
    options.instance.set(
        runid=cdf_name,
        shot_type=shot_type,
        input_time=input_time,
        input_points=None,
        uniform_rho=False,
        var_to_scan=None,
        scan_range=None)

    # Initialize variable objects and call simple_plot function
    mmm_vars, cdf_vars, raw_cdf_vars = datahelper.initialize_variables()

    simple_plot(cdf_vars.xb, cdf_vars.gnz, r' (CDF)', mmm_vars.xb, mmm_vars.gnz)
