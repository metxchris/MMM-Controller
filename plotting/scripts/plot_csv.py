#!/usr/bin/python3

"""Simple tool to plot data from CSV"""

# Standard Packages
import sys
sys.path.insert(0, '../')
sys.path.insert(0, '../../')
import io

# 3rd Party Packages
import numpy as np 
import matplotlib.pyplot as plt
import scipy.ndimage
from PyQt5.QtGui import QImage
from PyQt5.QtWidgets import QApplication

# Local Packages
from modules.options import Options
from plotting.modules.plotstyles import PlotStyles, StyleType
import modules.datahelper as datahelper
import modules.utils as utils


def load_csv(path, delimiter=""):
    data = np.genfromtxt(path, delimiter=delimiter)

    varx = []
    vary = []

    for i in range(data.shape[1]):
        if i % 2:
            vary.append(data[:, i])
        else:
            varx.append(data[:, i])

    print(varx[1], vary[1])

    return varx, vary

def plot_csv(varxs, varys, labels):

    for varx, vary, label in zip(varxs, varys, labels):
        plt.plot(varx, vary, label=label)

    plt.legend().set_draggable(state=True)
    


def on_press(event):
    fig, ax = plt.gcf(), plt.gca()

    if event.key == 'x':  # flip x-axis limits
        plt.xlim(plt.xlim()[::-1])
        fig.canvas.draw()

    if event.key == 'y':  # flip y-axis limits
        plt.ylim(plt.ylim()[::-1])
        fig.canvas.draw()

    if event.key == "ctrl+c":  # copy figure to clipboard
        save_format = plt.rcParams['savefig.format']
        plt.rcParams.update({'savefig.format': 'png'})
        with io.BytesIO() as buffer:
            fig.savefig(buffer)
            QApplication.clipboard().setImage(QImage.fromData(buffer.getvalue()))
            plt.rcParams.update({'savefig.format': save_format})


if __name__ == '__main__':

    PlotStyles(
        axes=StyleType.Axes.GRAY,
        lines=StyleType.Lines.MMM,
        layout=StyleType.Layout.SINGLE1B,
    )

    # options = Options(runid='138536A01', input_time=0.629, input_points=51)
    # mmm_vars, cdf_vars, __ = datahelper.initialize_variables(options)

    plt.figure()
    plt.gcf().canvas.mpl_connect('key_press_event', on_press)

    # csv_path = r"C:\Users\metxc\Documents\MMM-Explorer\plotting\scripts\data\fort.36"
    csv_path = r"C:\Users\metxc\Documents\MMM-Explorer\plotting\scripts\data\vei_dbm.csv"
    labels = ['Default', 'ETANC Formula']
    plot_csv(*load_csv(csv_path), labels)

    plt.title(r'Formula Comparison')
    plt.xlabel(r'$r$')
    plt.ylabel(r'$\nu_\mathrm{ei}$')

    plt.show()