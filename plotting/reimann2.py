import sys; sys.path.insert(0, '../')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Rectangle
import io
from PyQt5.QtGui import QImage
from PyQt5.QtWidgets import QApplication
from plotting.modules.plotstyles import PlotStyles, StyleType
import modules.utils as utils

RLIM = 10  # Upper x-value
LLIM = 1  # Lower x-value
SUBDIVISIONS = 10  # Number of rectangles
LEVELMULT = 100  # Larger = smoother curve
ROUND = 6  # Rounding numbers
EXP_FUNC_MULT = 1  # Argument of exponential exp_function
EXP_FUNC_ARG = -1.1  # Argument of exponential exp_function
YOFFSET = 0.01

# FUNC_COLOR = (0.094, 0.353, 0.663)
FUNC_COLOR = '#000'
RECT_LINE_COLOR = '#777'
RECT_FILL_COLOR = '#ddd'
DIFF_FILL_COLOR = (0.933, 0.180, 0.184)

USE_LINEAR = False
USE_LINEAR = True


class Riemann:
    def __init__(self, **kwargs) -> None:
        self.rlim = RLIM
        self.llim = max(LLIM, 1)
        self.slevel = SUBDIVISIONS * LEVELMULT  # Must be integer multiple of SUBDIVISIONS for rectangles to line up properly
        self.subd = SUBDIVISIONS
        self.yoffset = YOFFSET
        self.modifier_min = None
        self.modifer_exp = None

        self.set(**kwargs)

        self.rect_points = np.zeros((self.subd))
        self.rect_width = np.zeros((self.subd))

        if USE_LINEAR:
            self.rect_points = np.round(np.linspace(self.llim, self.rlim, self.subd + 1), ROUND)
            self.curve_points = np.round(np.linspace(self.llim, self.rlim, self.slevel + 1), ROUND)

        else:
            inc = (self.rlim / self.llim)**(1 / max(self.subd, 1))
            self.rect_points = np.round(max(self.llim, 1) * inc**(np.arange(0, self.subd + 1)), ROUND)

            if self.modifier_min is not None and self.modifer_exp is not None:

                modifier = (np.linspace(self.modifier_min, 1, self.rect_points.size))**self.modifer_exp
                self.rect_points = (self.rect_points - self.rect_points[0]) * modifier + self.rect_points[0]
                # print(self.rect_points)

            inc = (self.rlim / self.llim)**(1 / max(self.slevel, 1))
            self.curve_points = np.round(max(self.llim, 1) * inc**(np.arange(0, self.slevel + 1)), ROUND)

        # self.yoffset = np.round(self.exp_function(self.curve_points).max() / 20, OFFSET_ROUND)

        self.rect_width[:] = self.rect_points[1:] - self.rect_points[:-1]
        self.rect_function = self.exp_function(self.rect_points)
        self.approx_area = np.sum(self.rect_function[:-1] * self.rect_width)
        self.curve_function = self.exp_function(self.curve_points)

        self.fill_function = np.zeros_like(self.curve_function)
        self.fill_points = np.zeros_like(self.fill_function)
        self.fill_change = np.zeros_like(self.rect_function, dtype=int)



        # Exact formula, assuming no changes have been made
        # exact = EXP_FUNC_MULT * (self.exp_function(self.rlim) - self.exp_function(self.llim)) / EXP_FUNC_ARG + self.yoffset * (self.rlim - self.llim)
        
        # Exact area using trapezoid integration rule
        self.exact_area = np.trapz(self.exp_function(self.curve_points), x=self.curve_points)
        self.error = (self.approx_area / self.exact_area - 1) * 100  # Percent error formula

        self.func_str = rf'${{{EXP_FUNC_MULT}}}\,e^{{{EXP_FUNC_ARG}x}} + {self.yoffset}$'
        self.func_str = self.func_str.replace('e^{-1x}', 'e^{-x}').replace('{1}\,e', 'e')

        j = 0
        for i in range(0, len(self.fill_function)):
            self.fill_function[i] = self.rect_function[j]

            idx = int(min(i, len(self.fill_function) - 1))
            jdx = int(min(j + 1, len(self.rect_function) - 1))
            if round(self.curve_function[idx], ROUND) <= round(self.rect_function[jdx], ROUND):
                self.fill_change[jdx] = idx + 1
                j = min(j + 1, len(self.rect_function) - 1)

    def exp_function(self, points):
        return self.yoffset + EXP_FUNC_MULT * np.exp(EXP_FUNC_ARG * points)

    def set(self, **kwargs):
        '''Sets specified control values'''
        for key, value in kwargs.items():
            if not hasattr(self, key):
                raise ValueError(f'Invalid control specified: {key}')
            setattr(self, key, value)

    def save_to_csv(self, file_name=''):
        """
        Save plotted data to a CSV

        Data is saved in the order it is generated, so the first CSV column
        will be the first variable defined in FigData, etc.  Filenames
        are chosen as sequentially increasing integers.

        Raises:
        * NameError: If a file name can not be chosen
        * FileNotFoundError: If the file cannot be found after saving it
        """

        show_save_message = True

        if not file_name:  # Automatically generate the save name using a unique number
            show_save_message = True
            save_name_digits = 4
            save_dir = f'{utils.get_plotting_singles_path()}\\misc'
            utils.create_directory(save_dir)
            saved_files = utils.get_files_in_dir(save_dir, '*.csv')
            for save_number in range(1, 10**save_name_digits):
                file_name = f'{save_dir}\\{save_number:0>{save_name_digits}d}'
                if f'{file_name}.csv' not in saved_files:
                    break

            if not file_name:
                raise NameError('The filename for the CSV could not be set\n'
                                '\tMake sure Python has file reading permissions,'
                                ' and try deleting old CSVs from the CSV folder\n')
        
        title = 'Linear' if USE_LINEAR else 'Exponential'
        file_name = f'{file_name}_{SUBDIVISIONS}_{title}'
        save_dir = f'{utils.get_plotting_singles_path()}\\misc'
        file_name_curve = f'{save_dir}\\{file_name}_curve.csv'
        file_name_rect = f'{save_dir}\\{file_name}_rect.csv'

        output_curve = np.full((func.curve_function.shape[0], 2), np.nan, dtype=float)
        output_curve[:, 0] = func.curve_points
        output_curve[:, 1] = func.curve_function

        output_rect = np.full((func.rect_function.shape[0], 2), np.nan, dtype=float)
        output_rect[:, 0] = func.rect_points
        output_rect[:, 1] = func.rect_function

        prec, col_pad = 4, 7
        col_len = prec + col_pad
        fmt_str = f'%{col_len}.{prec}e'
        header = 'x, y'

        np.savetxt(file_name_curve, output_curve, header=header, fmt=fmt_str, delimiter=',')
        np.savetxt(file_name_rect, output_rect, header=header, fmt=fmt_str, delimiter=',')

        if not utils.check_exists(file_name_curve):
            raise FileNotFoundError('Failed to save plot data to a CSV in the CSV folder\n'
                                    '\tMake sure Python has file writing permissions')

        if not utils.check_exists(file_name_rect):
            raise FileNotFoundError('Failed to save plot data to a CSV in the CSV folder\n'
                                    '\tMake sure Python has file writing permissions')

        if show_save_message:
            print(
                f'Plot Data Saved:'
                f'\n\t{file_name_curve}'
                f'\n\t{file_name_rect}\n'
            )

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

    if event.key == 'alt+s':  # save plot lines to csv
        func.save_to_csv('reimann2')

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


def set_title(func):
    # title = f'Exact area: {func.exact_area:.3f}, Rectangle area: {approx:.3f}, Error: {(approx / func.exact_area - 1) * 100:.1f}%'
    # plt.suptitle(title, fontsize=8)

    # title = 'Linear' if USE_LINEAR else 'Exponential'
    # title = fr'Rectangle exp_function: {title}, Rectangle count: {func.subd}'
    # plt.title(title, fontsize=8)

    title = 'Linear' if USE_LINEAR else 'Exponential'
    title = f'{func.subd} {title} Segments, Error: {func.error:.1f}%'
    plt.title(title, fontsize=8)


def plot_subdivision(func, ax):
    # clear_prev_patches(ax)
    xCords, width = np.linspace(func.llim, func.rlim, func.subd + 1, retstep=True)
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
    lw = plt.rcParams['patch.linewidth'] if func.subd <= 400 else 0.1
    j = 0
    for i, height in zip(xCords, func_val):

        ax_rect = Rectangle((i, 0), width[j], height,
                            ec=RECT_LINE_COLOR,
                            fc=RECT_FILL_COLOR,
                            alpha=1,
                            linewidth=lw)
        ax.add_patch(ax_rect)
        j = j + 1

    set_title(func)


if __name__ == '__main__':

    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.AIP3,
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
        'patch.linewidth': 0.8,  # Width of rectangle lines
        'legend.handlelength': 1.5,  # length of lines
        'axes.linewidth': 0.5,
        'lines.linewidth': 1.0,
        'figure.subplot.left': 0.15,
        'legend.fontsize': 9,
    })

    mmin = np.linspace(0.6, 1, 10)
    mexp = np.linspace(1, 4, 4)
    print(mexp)

    best_error = 100
    best_min = -1
    best_exp = -1
    best_points = np.zeros_like((11, 1))

    for m in mmin:
        for n in mexp:
            func = Riemann(modifier_min=m, modifer_exp=n)
            if func.error < best_error:
                best_error = func.error
                best_min = m
                best_exp = n
                best_points = func.rect_points

             
    print(f'min: {best_min:.4f}, exp: {best_exp:.4f}, Error: {best_error:.4f}%')
    func = Riemann(modifier_min=best_min, modifer_exp=best_exp)   

    func = Riemann()   
    print(f'Exact area: {func.exact_area:.4f}, Rectangle area: {func.approx_area:.4f}, Error: {func.error:.4f}%')
    fig, ax = plt.subplots()
    fig.canvas.mpl_connect('key_press_event', on_press)
    graph, = ax.plot(func.curve_points, func.curve_function, color=FUNC_COLOR, linewidth=1.0, label=fr'f(x) = {func.func_str}')
    # ax.set(ylabel=fr'{func.func_str}')
    ax.set(xlabel=r'$x$')
    # ax.yaxis.set_ticks_position('right')
    ax.tick_params(bottom=False)  # Disable bottom ticks to not confuse with rectangle lines
    plt.legend()
    # plt.grid()
    ax.set_xlim(func.curve_points.min(), func.curve_points.max() )
    ax.set_ylim(min(func.curve_function.min(), 0), func.curve_function.max())
    # plt.subplots_adjust(top=.88, bottom=.1, right=.93, left=.1)

    # mng = plt.get_current_fig_manager()
    plot_subdivision(func, ax)

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
            alpha=0.5, color=DIFF_FILL_COLOR, interpolate=False,
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