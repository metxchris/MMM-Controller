import sys
sys.path.insert(0, '../')
# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from cycler import cycler
from plots.styles import standard as ps


def plot(input_options, x1var, y1var, l1='', x2var=None, y2var=None, l2=''):
    t_idx = input_options.time_idx
    
    fig = plt.figure(figsize=(3.5,3))
    plt.subplots_adjust(left=0.15, right=0.9, bottom=0.15, top=0.9)

    plt.plot(x1var.values[:, t_idx], y1var.values[:, t_idx], label=y1var.label + l1)
    if x2var is not None and y2var is not None:
        plt.plot(x2var.values[:, t_idx], y2var.values[:, t_idx], label=y2var.label + l2)

    plt.xlim(0, 1)
    plt.xlabel(x1var.label)
    plt.ylabel(y1var.units)
    plt.legend()
    plt.title(f'{input_options.runid}, t={input_options.time}s')
    plt.show()

if __name__ == '__main__':
    pass
