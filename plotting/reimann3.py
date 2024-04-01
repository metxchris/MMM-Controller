import sys; sys.path.insert(0, '../')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Rectangle
import io
from PyQt5.QtGui import QImage
from PyQt5.QtWidgets import QApplication
from plotting.modules.plotstyles import PlotStyles, StyleType
import modules.utils as utils
import matplotlib.transforms as mtransforms
import types


RLIM = 10  # Upper x-value
LLIM = 1  # Lower x-value
SUBDIVISIONS = 10  # Number of rectangles
LEVELMULT = 50  # Larger = smoother curve
ROUND = 16  # Rounding numbers
EXP_FUNC_MULT = 1  # Multiplier of exponential exp_function
EXP_FUNC_ARG = -1  # Argument of exponential exp_function
YOFFSET = 0.01

# FUNC_COLOR = (0.094, 0.353, 0.663)
FUNC_COLOR = '#333'
RECT_LINE_COLOR = '#444'
RECT_FILL_COLOR = '#ccc'
DIFF_FILL_COLOR = np.array([239, 87, 88]) / 255  # Red


class Riemann:
    def __init__(self, **kwargs) -> None:
        self.rlim = RLIM
        self.llim = max(LLIM, 1)
        self.slevel = SUBDIVISIONS * LEVELMULT  # Must be integer multiple of SUBDIVISIONS for rectangles to line up properly
        self.subd = SUBDIVISIONS
        self.yoffset = YOFFSET
        self.modifier_min = None
        self.modifier_exp = None
        self.modifier_exp2 = 1
        self.segments = 'Linear'

        self.set(**kwargs)
        self.segments = self.segments.capitalize()

        self.rect_points = np.zeros((self.subd))
        self.rect_width = np.zeros((self.subd))

        if self.segments == 'Linear':
            self.rect_points = np.round(np.linspace(self.llim, self.rlim, self.subd + 1), ROUND)
            self.curve_points = np.round(np.linspace(self.llim, self.rlim, self.slevel + 1), ROUND)

        else:
            inc = (self.rlim / self.llim)**(1 / (max(self.subd, 1)**self.modifier_exp2))
            self.rect_points = np.round(max(self.llim, 1) * inc**(np.arange(0, self.subd + 1)**self.modifier_exp2), ROUND)

            if self.modifier_min is not None and self.modifier_exp is not None:

                modifier = (np.linspace(self.modifier_min, 1, self.rect_points.size))**self.modifier_exp
                self.rect_points = (self.rect_points - self.rect_points[0]) * modifier + self.rect_points[0]
                # print(self.rect_points)


            inc = (self.rlim / self.llim)**(1 / max(self.slevel, 1))
            self.curve_points = np.round(max(self.llim, 1) * inc**(np.arange(0, self.slevel + 1)), ROUND)

        # self.yoffset = np.round(self.exp_function(self.curve_points).max() / 20, OFFSET_ROUND)

        self.rect_width[:] = self.rect_points[1:] - self.rect_points[:-1]
        self.rect_function = self.exp_function(self.rect_points)
        self.curve_function = self.exp_function(self.curve_points)

        self.fill_function = np.zeros_like(self.curve_function)
        self.fill_points = np.zeros_like(self.fill_function)
        self.fill_change = np.zeros_like(self.rect_function, dtype=int)



        # Exact formula, assuming no changes have been made
        # exact = EXP_FUNC_MULT * (self.exp_function(self.rlim) - self.exp_function(self.llim)) / EXP_FUNC_ARG + self.yoffset * (self.rlim - self.llim)
        
        # Areas
        self.approx_area = np.sum(self.rect_function[:-1] * self.rect_width)
        self.exact_area = np.trapz(self.exp_function(self.curve_points), x=self.curve_points)
        self.error = (self.approx_area / self.exact_area - 1) * 100  # Percent error formula

        self.func_str = rf'${{{EXP_FUNC_MULT}}}\,e^{{{EXP_FUNC_ARG}x}} + {self.yoffset}$'
        self.func_str = self.func_str.replace('e^{-1x}', 'e^{-x}').replace('{1}\,e', 'e')

        j = 1
        fill_len = len(self.fill_function)
        for i in range(0, fill_len):
            self.fill_function[i] = self.rect_function[j - 1]

            if self.curve_function[i] <= self.rect_function[j] and self.curve_function[i] < self.rect_function[j - 1]:
                self.fill_change[j] = i + 1
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
        
        file_name = f'{file_name}_{SUBDIVISIONS}_{self.segments}'
        save_dir = f'{utils.get_plotting_singles_path()}\\misc'

        file_name_curve = f'{save_dir}\\{file_name}_curve.csv'
        file_name_rect = f'{save_dir}\\{file_name}_rect.csv'

        output_curve = np.vstack((func.curve_points, func.curve_function)).T
        output_rect = np.vstack((func.rect_points, func.rect_function)).T

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
        func.save_to_csv('reimann3')

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

def get_title(func):
    ncount = r'$n_{{\rm XXX}} = $'.replace('XXX', f'{func.segments[0]}')
    return f'{ncount}{func.subd}, Error: {func.error:.1f}%'

def plot_fill_rects(func, ax):

    gray_fill = np.linspace(0.7, 0.9, len(func.fill_change))
    gray_line = np.linspace(0.2, 0.5, len(func.fill_change))
    error_fill_exp = np.linspace(0.9, 0.7, len(func.fill_change))

    for i, k in enumerate(func.fill_change[1:]):
 
        k0 = max(func.fill_change[i] - 1, 0)

        # Fill in error rectangles above curve
        ax.fill_between(
            func.curve_points[k0:k], 
            func.curve_function[k0:k], 
            np.hstack((func.fill_function[k0 + 1], func.fill_function[k0 + 1:k])), 
            ec=np.array(DIFF_FILL_COLOR)**10,
            # ec=np.full(3, gray_line[i]),
            fc=np.array(DIFF_FILL_COLOR)**error_fill_exp[i],
            alpha=1,
            interpolate=False,
        )

        # Fill in rectangles below curve
        ax.fill_between(
            func.curve_points[k0:k], 
            func.curve_function[k0:k], 
            np.zeros_like(func.curve_function[k0:k]), 
            alpha=1,
            ec=np.full(3, gray_line[i]),
            fc=np.full(3, gray_fill[i]),
            interpolate=False,
        )

def monkey_patch(axis, func):
    axis._update_offset_text_position = types.MethodType(func, axis)

def x_update_offset_text_position(self, bboxes, bboxes2):
    x, y = self.offsetText.get_position()

    if self.offset_text_position == 'bottom':
        if bboxes:
            bbox = mtransforms.Bbox.union(bboxes)
        else:
            bbox = self.axes.bbox
        y = bbox.ymin - self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        self.offsetText.set(va='top', ha='right')

    else:
        if bboxes2:
            bbox = mtransforms.Bbox.union(bboxes2)
        else:
            bbox = self.axes.bbox
        y = bbox.ymax + self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        self.offsetText.set(va='bottom', ha='left')

    self.offsetText.set_position((x, y))

def y_update_offset_text_position(self, bboxes, bboxes2):
    x, y = self.offsetText.get_position()

    if self.offset_text_position == 'left':
        # y in axes coords, x in display coords
        self.offsetText.set_transform(mtransforms.blended_transform_factory(
                self.axes.transAxes, mtransforms.IdentityTransform()))

        top = self.axes.bbox.ymax

        x = -0.13
        # x = 0 - self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        y = top + self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        
    else:
        # x & y in display coords
        self.offsetText.set_transform(mtransforms.IdentityTransform())

        # Northwest of upper-right corner of right-hand extent of tick labels
        if bboxes2:
            bbox = mtransforms.Bbox.union(bboxes2)
        else:
            bbox = self.axes.bbox
        top, right = bbox.ymax, bbox.xmax
        x = right + self.OFFSETTEXTPAD * self.figure.dpi / 72.0
        y = top + self.OFFSETTEXTPAD * self.figure.dpi / 72.0

    self.offsetText.set_position((x, y))



if __name__ == '__main__':

    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.AIP3,
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
        'patch.linewidth': 0.5,  # Width of rectangle lines
        'legend.handlelength': 1.5,  # length of lines
        'axes.linewidth': 0.5,
        'lines.linewidth': 1.2,
        'figure.subplot.left': 0.15,
        'legend.fontsize': 9,
    })

    modifier_min = np.linspace(0.5, 1, 300)
    modifier_exp = np.linspace(1, 4, 4)
    modifier_exp2 = np.linspace(1, 4, 100)

    best_error = 100
    best_min = -1
    best_exp = -1
    best_points = np.zeros_like((11, 1))
    best_func = None

    segments = 'Linear'
    # segments = 'Exp.'

    func = Riemann(segments=segments)
    print(f'Exact area: {func.exact_area:.4f}, Rectangle area: {func.approx_area:.4f}, Error: {func.error:.4f}%')

    # if segments != 'Linear':
        # for m in modifier_min:
        #     for n in modifier_exp:
        #         func = Riemann(modifier_min=m, modifier_exp=n, segments='Exp.')
        #         if func.error < best_error:
        #             best_error = func.error
        #             best_min = m
        #             best_exp = n
        #             best_points = func.rect_points
        #             best_func = func

        # for n in modifier_exp2:
        #     func = Riemann(modifier_exp2=n, segments='Exp.')
        #     if func.error < best_error:
        #         best_error = func.error
        #         best_exp = n
        #         best_points = func.rect_points
        #         best_func = func

    if best_func:
        print(f'\tmin: {best_min:.4f}, exp: {best_exp:.4f}, Error: {best_error:.4f}%')

    func = best_func or Riemann(segments=segments)
    print(f'Exact area: {func.exact_area:.4f}, Rectangle area: {func.approx_area:.4f}, Error: {func.error:.4f}%')

    fig, ax = plt.subplots()
    fig.canvas.mpl_connect('key_press_event', on_press)
    ax.plot(func.curve_points, func.curve_function, color=FUNC_COLOR, linewidth=plt.rcParams['lines.linewidth'], label=fr'{func.func_str}')
    monkey_patch(ax.xaxis, x_update_offset_text_position)
    monkey_patch(ax.yaxis, y_update_offset_text_position)
    ax.set(xlabel=r'$x$')
    plt.title(get_title(func), fontsize=8)
    ax.tick_params(bottom=False, left=False)  # Disable bottom ticks to not confuse with rectangle lines
    plt.legend()
    ax.set_xlim(func.curve_points.min(), func.curve_points.max() )
    ax.set_ylim(min(func.curve_function.min(), 0), func.curve_function.max()*1.05)

    plot_fill_rects(func, ax)




    plt.show()
