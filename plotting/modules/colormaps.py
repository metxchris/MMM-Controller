# 3rd Party Packages
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


colormaps = {}


def _init_colormaps():
    colormap_key = 'magma_positive'
    cmap = cm.get_cmap('magma_r', 60)
    colors = np.array(cmap(np.arange(0, cmap.N))[:-10])
    cmap = LinearSegmentedColormap.from_list(colormap_key, colors, N=256)
    cmap.set_under([1, 1, 1, 0])
    colormaps[colormap_key] = cmap

    colormap_key = 'magma_both'
    cmap = cm.get_cmap('binary_r', 30)
    cmap2 = cm.get_cmap('magma_r', 30)
    colors = np.vstack((
        cmap(np.arange(0, cmap.N))[5:],
        # np.array([1, 1, 1, 1]),
        cmap2(np.arange(0, cmap2.N))[:-5]
    ))
    cmap = LinearSegmentedColormap.from_list(colormap_key, colors, N=512)
    colormaps[colormap_key] = cmap

    colormap_key = 'magma_negative'
    cmap = cm.get_cmap('binary', 200)
    colors = np.vstack((np.array([1, 1, 1, 0]), cmap(np.arange(0, cmap.N))[180:]))
    cmap = LinearSegmentedColormap.from_list('colormap_key', colors, N=256)
    colormaps[colormap_key] = cmap

    colormap_key = 'binary'
    cmap = cm.get_cmap('binary', 200)
    colors = np.array(cmap(np.arange(0, cmap.N)))[180:]
    cmap = LinearSegmentedColormap.from_list(colormap_key, colors, N=256)
    colormaps[colormap_key] = cmap


def get_colormaps():
    if not colormaps:
        _init_colormaps()

    return colormaps
