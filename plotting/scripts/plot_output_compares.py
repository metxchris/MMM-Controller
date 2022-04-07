#!/usr/bin/python3

"""Saved output comparison plots for the ETGM paper"""

# Standard Packages
import sys; sys.path.insert(0, '../'); sys.path.insert(0, '../../')
import logging

# 3rd Party Packages
import matplotlib.pyplot as plt

# Local Packages
import modules.utils as utils
from plotting.modules.plotstyles import PlotStyles, StyleType
from plotting.plot_variables import AllPlotData, PlotDataCsv, main


_log = logging.getLogger(__name__)


def plot_profiles(profile_list, all_data_list, saveall=True):
    for profiles, all_data in zip(profile_list, all_data_list):
        for p in profiles:
            all_data.set(*p)
            main(all_data, savefig=saveall, savedata=saveall)


def get_compared_profiles(ptype, r, scan_nums, legends, savename_append, title=''):

    if ptype == 'max':
        var_names = [
            'gmaETGM', 'omgETGM', 'xteETGM', 'xte2ETGM',
            'kyrhosETGM', 'gaveETGM', 'phi2ETGM', 'Apara2ETGM', 'satETGM',
        ]
        all_data_max.title_override = title
        all_data_max.savename_append = savename_append
    else:
        var_names = ['xteETGM', 'xte2ETGM']
        all_data_sum.title_override = title
        all_data_sum.summed_modes = True
        all_data_sum.savename_append = f'sum_{savename_append}'

    profiles = []

    for var_name in var_names:
        profiles.append(
            [PlotDataCsv(r, scan_num, var_name, legend=legend) for scan_num, legend in zip(scan_nums, legends)]
        )

    return profiles


if __name__ == '__main__':
    """Run this module directly to plot variable data stored in CDF or CSV files"""

    utils.init_logging()

    # Initialize visual styles for the plot
    PlotStyles(
        axes=StyleType.Axes.WHITE,
        lines=StyleType.Lines.RHO_MMM,
        layout=StyleType.Layout.SINGLE1B,
    )

    plt.rcParams.update({
        'savefig.format': 'pdf',  # Common save formats: png, pdf, eps
    })

    all_data_max = AllPlotData(
        replace_offset_text=False,
        allow_title_runid=False,
        allow_title_time=False,
        allow_title_factor=True,
        allow_title_rho=True,
        invert_y_axis=False,
        invert_x_axis=False,
        nomralize_y_axis=False,
        nomralize_x_axis=False,
        savename_append='',
        title_override=' ',
        ylabel_override='',
        xlabel_override='',
    )

    all_data_sum = AllPlotData(
        replace_offset_text=False,
        allow_title_runid=False,
        allow_title_time=False,
        allow_title_factor=True,
        allow_title_rho=True,
        invert_y_axis=False,
        invert_x_axis=False,
        nomralize_y_axis=False,
        nomralize_x_axis=False,
        savename_append='',
        title_override=' ',
        ylabel_override='',
        xlabel_override='',
    )

    r = '138536A01'  # Discharge
    nmax, nsum = 1787, 1788  # Scan numbers of default profiles

    # Plot legend labels
    ahyd = r'$a_\mathrm{hyd}$'
    kpara2 = r'$\langle k^2_\parallel\rangle\propto$'
    kxky = r'$k_\mathrm{x}/k_\mathrm{y}$'
    kyrhos = r'$k_y\rho_\mathrm{s}$'
    kyrhosmin = r'$(k_y\rho_\mathrm{s})_\mathrm{min}$'
    wexb = r'$\omega_{E \times B}$'
    bp = r'$\beta^\prime$'

    # COMPARED PROFILES
    pmax = get_compared_profiles('max', r, [nmax], [''], 'base')
    psum = get_compared_profiles('sum', r, [nsum], [''], 'base')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1741, 1745], [fr'{kpara2}$1/3(qR)^2$', fr'{kpara2}$1/5.3(qR)^2$', fr'{kpara2}$1/(qR)^2$'], 'kpara2')
    # psum = get_compared_profiles('sum', r, [nsum, 1742, 1746], [fr'{kpara2}$1/3(qR)^2$', fr'{kpara2}$1/5.3(qR)^2$', fr'{kpara2}$1/(qR)^2$'], 'kpara2')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1696], [r'using $\hat{s}_{\nabla\rho}$', r'using $s$'], 's')
    # psum = get_compared_profiles('sum', r, [nsum, 1697], [r'using $\hat{s}_{\nabla\rho}$', r'using $s$'], 's')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1700], ['Default', r'$\beta_{\mathrm{e}} = 0,\,\alpha_\mathrm{m}\,\, \mathrm{unchanged}$'], 'betae0')
    # psum = get_compared_profiles('sum', r, [nsum, 1701], ['Default', r'$\beta_{\mathrm{e}} = 0,\,\alpha_\mathrm{m}\,\, \mathrm{unchanged}$'], 'betae0')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1670], ['Default', r'$g_\mathrm{Bu}=1$'], 'gbu')
    # psum = get_compared_profiles('sum', r, [nsum, 1671], ['Default', r'$g_\mathrm{Bu}=1$'], 'gbu')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1674, 1676], [r'using $B_\mathrm{unit}$', r'using $B_\mathrm{unit}, \overline{G}=1$', r'using $B_\phi,\, \overline{G}=1$'], 'btorgave1')
    # psum = get_compared_profiles('sum', r, [nsum, 1675, 1677], [r'using $B_\mathrm{unit}$', r'using $B_\mathrm{unit}, \overline{G}=1$', r'using $B_\phi,\, \overline{G}=1$'], 'btorgave1')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1672], [r'using $B_\mathrm{unit}$', r'using $B_\phi$'], 'btor')
    # psum = get_compared_profiles('sum', r, [nsum, 1673], [r'using $B_\mathrm{unit}$', r'using $B_\phi$'], 'btor')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1666], ['Collisional', 'Collisionless'], 'cl')
    # psum = get_compared_profiles('sum', r, [nsum, 1667], ['Collisional', 'Collisionless'], 'cl')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1668], ['Electromagnetic', 'Electrostatic'], 'es')
    # psum = get_compared_profiles('sum', r, [nsum, 1669], ['Electromagnetic', 'Electrostatic'], 'es')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1678, 1680], [fr'{kxky}$=0.5$', fr'{kxky}$=0.4$', fr'{kxky}$=0.3$'], 'kxky')
    # psum = get_compared_profiles('sum', r, [nsum, 1679, 1681], [fr'{kxky}$=0.5$', fr'{kxky}$=0.4$', fr'{kxky}$=0.3$'], 'kxky')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1682, 1684], [fr'{kyrhos} scan', fr'{kyrhos}$=5$', fr'{kyrhos}$=10$'], 'kyrhos')
    # psum = get_compared_profiles('sum', r, [nsum, 1683, 1685], [fr'{kyrhos} scan', fr'{kyrhos}$=5$', fr'{kyrhos}$=10$'], 'kyrhos')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # legend = [r'$n=25$', r'$n=50$', r'$n=100$']
    # pmax = get_compared_profiles('max', r, [1816, 1817, 1818], legend, 'kyrhoscount')
    # psum = get_compared_profiles('sum', r, [1830, 1831, 1832], legend, 'kyrhoscount')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # legend = [r'Converged', r'$n=500$, Linear', r'$n=50,\,\,$ Exponential']
    # pmax = get_compared_profiles('max', r, [1822, 1826, 1828], legend, 'kyrhosinc')
    # psum = get_compared_profiles('sum', r, [1823, 1827, 1829], legend, 'kyrhosinc')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1686], [fr'{wexb} off', fr'{wexb} on'], 'wexb')
    # psum = get_compared_profiles('sum', r, [nsum, 1687], [fr'{wexb} off', fr'{wexb} on'], 'wexb')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [1833, 1835], [fr'{wexb} off', fr'{wexb} on'], 'wexbkyrhos', title=fr'Minimum {kyrhos}$=0.1$')
    # psum = get_compared_profiles('sum', r, [1834, 1836], [fr'{wexb} off', fr'{wexb} on'], 'wexbkyrhos', title=fr'Minimum {kyrhos}$=0.1$')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1688], ['Default', r'$g_{\mathrm{ne}}=0$'], 'gne0')
    # psum = get_compared_profiles('sum', r, [nsum, 1689], ['Default', r'$g_{\mathrm{ne}}=0$'], 'gne0')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1690], ['Default', r'$g_{\mathrm{Te}}=0$'], 'gte0')
    # psum = get_compared_profiles('sum', r, [nsum, 1691], ['Default', r'$g_{\mathrm{Te}}=0$'], 'gte0')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1692], ['Default', r'$g_{\mathrm{ne}}, g_{\mathrm{Te}}=0$'], 'gnegte0')
    # psum = get_compared_profiles('sum', r, [nsum, 1693], ['Default', r'$g_{\mathrm{ne}}, g_{\mathrm{Te}}=0$'], 'gnegte0')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1694], [r'Saturation$=2$', r'Saturation$=1$'], 'sat')
    # psum = get_compared_profiles('sum', r, [nsum, 1695], [r'Saturation$=2$', r'Saturation$=1$'], 'sat')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1712, 1714], [fr'{kyrhosmin}$ = 1.00$', fr'{kyrhosmin}$ = 0.50$', fr'{kyrhosmin}$ = 0.25$'], 'kyrhosmin')
    # psum = get_compared_profiles('sum', r, [nsum, 1713, 1715], [fr'{kyrhosmin}$ = 1.00$', fr'{kyrhosmin}$ = 0.50$', fr'{kyrhosmin}$ = 0.25$'], 'kyrhosmin')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1708, 1710], [fr'{ahyd}$=2$ (Default)', fr'{ahyd}$=1$', fr'{ahyd}$=3$'], 'ahyd')
    # psum = get_compared_profiles('sum', r, [nsum, 1709, 1711], [fr'{ahyd}$=2$ (Default)', fr'{ahyd}$=1$', fr'{ahyd}$=3$'], 'ahyd')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [1853, 1852, 1854], [fr'{ahyd}$=2$ (Default)', fr'{ahyd}$=1$', fr'{ahyd}$=3$'], 'ahyd2')
    # psum = get_compared_profiles('sum', r, [1855, 1856, 1857], [fr'{ahyd}$=2$ (Default)', fr'{ahyd}$=1$', fr'{ahyd}$=3$'], 'ahyd2')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # pmax = get_compared_profiles('max', r, [nmax, 1754], [fr'Default', fr'Geometry Disabled'], 'geo')
    # psum = get_compared_profiles('sum', r, [nsum, 1755], [fr'Default', fr'Geometry Disabled'], 'geo')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # legend = [fr'$1\,${bp}', fr'$0\,${bp}', fr'$2\,${bp}']
    # pmax = get_compared_profiles('max', r, [nmax, 1789, 1791], legend, 'betaep')
    # psum = get_compared_profiles('sum', r, [nsum, 1790, 1792], legend, 'betaep')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # legend = [r'kpc$=1$', r'kpc$ = c \langle k_\parallel \rangle / (\omega_\mathrm{De} / \overline{G})$', r'kpc$ = c \langle k_\parallel \rangle / \omega_\mathrm{De}$']
    # pmax = get_compared_profiles('max', r, [nmax, 1811, 1809], legend, 'kpc')
    # psum = get_compared_profiles('sum', r, [nsum, 1812, 1810], legend, 'kpc')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # Change plot boundary to not include first variable values
    # legend = [r'using $|\hat{\phi}|$', r'using $|\hat{A}_{\!\parallel}\!|$', r'using $|\hat{\phi}|+|\hat{A}_{\!\parallel}\!|$']
    # pmax = get_compared_profiles('max', r, [nmax, 1702, 1704], legend, 'phi')
    # psum = get_compared_profiles('sum', r, [nsum, 1703, 1705], legend, 'phi')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    """
    Testing isotropic mass change
    """
    # pmax = get_compared_profiles('max', r, [1708, 1847], [fr'{ahyd}$=1 $ (old)', fr'{ahyd}$=1 $ (new)'], 'ahyd')
    # psum = get_compared_profiles('sum', r, [1709, 1849], [fr'{ahyd}$=1 $ (old)', fr'{ahyd}$=1 $ (new)'], 'ahyd')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=0)
    # pmax = get_compared_profiles('max', r, [1710, 1850], [fr'{ahyd}$=3 $ (old)', fr'{ahyd}$=3 $ (new)'], 'ahyd')
    # psum = get_compared_profiles('sum', r, [1711, 1851], [fr'{ahyd}$=3 $ (old)', fr'{ahyd}$=3 $ (new)'], 'ahyd')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=0)
