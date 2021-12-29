# Standard Packages
import sys; sys.path.insert(0, '../')

# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt

# Local Packages
import modules.datahelper as datahelper
from modules.options import Options
from modules.enums import ShotType
from plotting.modules.styles import single as plotlayout
from plotting.modules.colors import mmm as plotcolors


def simple_plot(title, t1, t2, x1var, y1var, l1='', x2var=None, y2var=None, l2=''):

    plotlayout.init()
    plotcolors.init()

    plt.figure(figsize=(3.5, 3))
    plt.subplots_adjust(left=0.15, right=0.9, bottom=0.15, top=0.9)

    if x1var.values is not None and y1var.values is not None:
        plt.plot(x1var.values[:, t1], y1var.values[:, t1], label=y1var.label + l1)
    if x2var is not None and y2var is not None:
        plt.plot(x2var.values[:, t2], y2var.values[:, t2], label=y2var.label + l2)

    xmin = min(x1var.values.min(), x2var.values.min()) if x2var is not None else x1var.values.min()
    xmax = max(x1var.values.max(), x2var.values.max()) if x2var is not None else x1var.values.max()

    plt.xlim(xmin, xmax)
    plt.xlabel(x1var.label)
    plt.ylabel(y1var.units_label)
    plt.legend()
    plt.title(title)
    plt.show()


# Run this file directly to make a simple plot of variable profiles
if __name__ == '__main__':

    vars1, __, __ = datahelper.initialize_variables(Options(runid='120982A09', input_time=0.50))
    t1 = vars1.options.time_idx

    vars2, __, __ = datahelper.initialize_variables(Options(runid='129041A10', input_time=0.45))
    t2 = vars2.options.time_idx

    # title = f'{options.runid}, t={options.time_str}s'
    title = 'Electron Temperature'
    label1 = f' {vars1.options.runid} t={vars1.options.time_str}s'
    label2 = f' {vars2.options.runid} t={vars2.options.time_str}s'
    simple_plot(title, t1, t2, vars1.rho, vars1.te, label1, vars2.rho, vars2.te, label2)
