# Standard Packages
import sys; sys.path.insert(0, '../'); sys.path.insert(0, '../../')
import logging
import io

# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Rectangle
from PyQt5.QtGui import QImage
from PyQt5.QtWidgets import QApplication

# Local Packages
import modules.options
# import modules.constants as constants
# import modules.datahelper as datahelper
import modules.utils as utils
from modules.enums import SaveType
from plotting.modules.plotstyles import PlotStyles, StyleType


_log = logging.getLogger(__name__)


RLIM = 10.5
LLIM = 0.5
SUBDIVISIONS = 10
LEVELMULT = 100
SLEVEL = SUBDIVISIONS * LEVELMULT  # Must be integer multiple of SUBDIVISIONS for rectangles to line up properly
ROUND = 8

FUNC_COLOR = (0.094, 0.353, 0.663)
FUNC_COLOR = '#000'
RECT_LINE_COLOR = '#aaa'
RECT_FILL_COLOR = '#ddd'
DIFF_EDGE_COLOR = '#E77474'
DIFF_FILL_COLOR = '#E2a8a8'

USE_LINEAR = False
# USE_LINEAR = True


class Function:
    def __init__(self) -> None:
        self.rlim = RLIM
        self.llim = LLIM
        self.slevel = SLEVEL  # how smooth will be the plotted curve
        self.subd = SUBDIVISIONS
        self.yoffset = 0.04

        self.rect_points = np.zeros((self.subd + 1))
        self.rect_width = np.zeros((self.subd))

        if USE_LINEAR:
            self.rect_points = np.round(np.linspace(self.llim, self.rlim, self.subd + 1), ROUND)
            self.curve_points = np.round(np.linspace(self.llim, self.rlim, self.slevel + 1), ROUND)
            
        else:
            inc = (self.rlim / max(self.llim, 1e-6))**(1 / max(self.subd, 1))
            self.rect_points = np.round(max(self.llim, 1e-6) * inc**(np.arange(0, self.subd + 1)), ROUND)

            inc = (self.rlim / max(self.llim, 1e-6))**(1 / max(self.slevel, 1))
            self.curve_points = np.round(max(self.llim, 1e-6) * inc**(np.arange(0, self.slevel + 1)), ROUND)

        self.rect_width[:] = self.rect_points[1:] - self.rect_points[:-1]
        # self.curve_points = np.linspace(self.llim, self.rlim, self.slevel + 1)
        self.rect_function = self.function(self.rect_points)
        self.curve_function = self.function(self.curve_points)


        self.fill_function = np.zeros_like(self.curve_function)
        self.fill_points = np.zeros_like(self.fill_function)
        self.fill_change = np.zeros_like(self.rect_function, dtype=int)

        j = 0
        for i in range(0, len(self.fill_function)):
            
            self.fill_function[i] = self.rect_function[j]

            idx = int(min(i, len(self.fill_function) - 1))
            jdx = int(min(j + 1, len(self.rect_function) - 1))
            if round(self.curve_function[idx], ROUND) <= round(self.rect_function[jdx], ROUND):
                self.fill_change[jdx] = idx + 1
                j = j + 1

        self.exact_area = -self.function(self.rlim) + self.function(self.llim) + self.yoffset * (self.rlim - self.llim)

        func_str = 'e^{{-x}}'
        self.func_str = rf'${func_str} + {self.yoffset}$'


    def function(self, points):
        return self.yoffset + np.exp(-points)

def on_press(event):
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


def clear_prev_patches(ax):
    patches = ax.patches
    while len(patches) != 0:
        for each_patch in patches:
            each_patch.remove()


def approximate_area(func_val, width):
    height_sum = 0
    for each in func_val:
        # print(each)
        height_sum += 0 if np.isnan(each) else each
    approx_area = height_sum*width
    return approx_area


def set_title(func, approx):
    title = fr'$A = ${func.exact_area:.3f}, $A_{{\rm rect}} = ${approx:.3f}, Error$ = ${(approx / func.exact_area - 1) * 100:.1f}%'
    plt.title(title, fontsize=8)

    title = 'Linear' if USE_LINEAR else 'Exponential'
    # title = fr'Rectangle function: {title}, Rectangle count: {func.subd}'
    title = fr'{func.subd} {title} Increments'
    plt.suptitle(title, fontsize=9)


def plot_subdivision(func, ax):
    # clear_prev_patches(ax)
    xCords, width = np.linspace(
        func.llim, func.rlim, func.subd+1, retstep=True)
    # print('xcords')
    # print(xCords)
    # print(func.rect_points)
    # print(width)
    # print(func.rect_width)


    xCords = func.rect_points
    width = func.rect_width
    # print(width[0])
    xCords = xCords[:-1]
    func_val = func.rect_function
    lw = 1.0 if func.subd <= 400 else 0.1
    j = 0
    for i, height in zip(xCords, func_val):

        ax_rect = Rectangle((i, 0), width[j], height,
                            ec=RECT_LINE_COLOR,
                            fc=RECT_FILL_COLOR,
                            alpha=1,
                            linewidth=lw)
        ax.add_patch(ax_rect)
        j = j + 1
    # approx_area = approximate_area(func_val[:-1], width[0])
    approx_area = np.sum(func.rect_function[:-1] * func.rect_width)
    # print('---')
    # print(func.rect_function[:-1])
    # print(func.rect_width)
    set_title(func, approx_area)


if __name__ == '__main__':
    utils.init_logging()

    # Initialize visual styles for the plot
    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.AIP2,
    )

    func = Function()
    fig, ax = plt.subplots()
    fig.canvas.mpl_connect('key_press_event', on_press)
    graph, = ax.plot(func.curve_points, func.curve_function, color=FUNC_COLOR, linewidth=1.5, label=fr'{func.func_str}')
    plt.legend()
    # plt.grid()
    ax.set_xlim(func.curve_points.min(), func.curve_points.max())
    ax.set_ylim(min(func.curve_function.min(), 0), func.curve_function.max())
    plt.xticks(np.unique(np.linspace(max(int(func.llim), 1), int(func.rlim), func.subd, dtype=int)))
    plt.subplots_adjust(top=.84, bottom=.12, right=.93, left=.12)

    # mng = plt.get_current_fig_manager()
    plot_subdivision(func, ax)
    # ax.set(xlabel='x', ylabel='y')

    plt.draw()

    # s0 = 1
    s = 1
    k0 = max(func.fill_change[s] - 1, 0)
    if s + 1 >= len(func.fill_change):
        k = len(func.curve_function)
    else:
        k = func.fill_change[s + 1]
    
    # print(s, k0, k)
    # print(func.curve_points[k0 - 1], func.curve_points[k0], func.rect_points[s])
    # ax.fill_between(
    #     func.curve_points[k0:k], 
    #     func.curve_function[k0:k], 
    #     np.hstack((func.fill_function[k0 + 1], func.fill_function[k0 + 1:k])), 
    #     alpha=1, color=DIFF_FILL_COLOR
    # )

    # print(func.curve_points[k-1], func.rect_points[s+1])

    s = 0
    for k in func.fill_change:
        if k < func.llim:
            continue
        # print(s)
        k0 = max(func.fill_change[s] - 1, 0)
        if s + 1 >= len(func.fill_change):
            k = len(func.curve_function)
        else:
            k = func.fill_change[s + 1]

        # print(func.curve_points[k - 1], func.curve_points[k - 2], func.rect_points[s + 1])
        ax.fill_between(
            func.curve_points[k0:k], 
            func.curve_function[k0:k], 
            np.hstack((func.fill_function[k0 + 1], func.fill_function[k0 + 1:k])), 
            alpha=1, color=DIFF_FILL_COLOR, hatch='/////', ec=DIFF_EDGE_COLOR, linewidth=1,
        )
        s = s + 1
        # print(k, len(func.curve_function))
        if k == len(func.curve_function):
            break
        # print ('***')
        # print(s, k, k-1)


    plt.show()
    # mng.window.showMaximized()
    # create_tboxes(func, graph, ax)