import types

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.transforms as mtransforms

def main():
    x = np.linspace(0, 1e8)
    y = np.exp(x**0.2)

    fig, ax = plt.subplots()
    ax.plot(x, y)

    monkey_patch(ax.xaxis, x_update_offset_text_position)
    monkey_patch(ax.yaxis, y_update_offset_text_position)

    ax.xaxis.tick_top()
    ax.yaxis.tick_right()
    ax.xaxis.offset_text_position = "top"
    ax.yaxis.offset_text_position = "right"

    plt.show()

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

main()