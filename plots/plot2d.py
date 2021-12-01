# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

def plot(input_options, x1var, y1var, l1='', x2var=None, y2var=None, l2=''):
    plt.rcParams.update({
        'font.family': 'serif',
        'mathtext.fontset': 'cm',
        'font.serif': 'cmr10',
        'axes.unicode_minus': False # (unicode_minus does not work in Computer Modern font)
    })

    t_idx = input_options.time_idx
    
    fig = plt.figure(figsize=(5,3.75))

    plt.plot(x1var.values[:, t_idx], y1var.values[:, t_idx], 'b-', lw=2, label=y1var.label + l1)
    if x2var is not None and y2var is not None:
        plt.plot(x2var.values[:, t_idx], y2var.values[:, t_idx], 'r-.', lw=2.5, label=y2var.label + l2)

    plt.title('{0}, t={1}s'.format(input_options.runid, input_options.time), fontsize=14)
    plt.xlabel(x1var.label, labelpad=0, fontsize=12)
    plt.ylabel(y1var.units, labelpad=0, fontsize=12)
    plt.legend(borderpad=0, labelspacing=0, frameon=False, prop=fm.FontProperties(size=14))
    plt.show()

if __name__ == '__main__':
    pass
