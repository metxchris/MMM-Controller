# 3rd Party Packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

def plot2d(x1, y1, x2=None, y2=None):
    plt.figure()
    plt.plot(x1,y1, 'b-', lw=2)
    if x2 is not None and y2 is not None:
        plt.plot(x2, y2, 'r-.', lw=2.5, mew=1.5, fillstyle='none', markersize=4)
    plt.show()

if __name__ == '__main__':
    pass
