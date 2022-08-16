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
        print(f'Plotting profiles for {all_data.savename_append}')
        for p in profiles:
            all_data.set(*p)
            main(all_data, savefig=saveall, savedata=saveall)


def get_compared_profiles(ptype, r, scan_nums, legends, savename_append, title=' '):

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
        layout=StyleType.Layout.AIP,
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
        ylabel_override='',
        xlabel_override='',
    )

    r = '138536A01'  # Discharge
    nmax, nsum = 1787, 1788  # Scan numbers of default profiles

    # Plot legend labels
    ah = r'$a_\mathrm{hyd}$'
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

    legend = [fr'{kpara2}$1/3(qR)^2$', fr'{kpara2}$1/5.3(qR)^2$', fr'{kpara2}$1/(qR)^2$']
    pmax = get_compared_profiles('max', r, [nmax, 1741, 1745], legend, 'kpara2')
    psum = get_compared_profiles('sum', r, [nsum, 1742, 1746], legend, 'kpara2')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [r'using $\hat{s}_{\nabla\rho}$', r'using $s$']
    pmax = get_compared_profiles('max', r, [nmax, 1696], legend, 's')
    psum = get_compared_profiles('sum', r, [nsum, 1697], legend, 's')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = ['Default', r'$\beta_{\mathrm{e}} = 0,\,\alpha_\mathrm{m}\,\, \mathrm{unchanged}$']
    pmax = get_compared_profiles('max', r, [nmax, 1700], legend, 'betae0')
    psum = get_compared_profiles('sum', r, [nsum, 1701], legend, 'betae0')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = ['Default', r'$g_\mathrm{Bu}=1$']
    pmax = get_compared_profiles('max', r, [nmax, 1670], legend, 'gbu')
    psum = get_compared_profiles('sum', r, [nsum, 1671], legend, 'gbu')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [r'using $B_\mathrm{unit}$', r'using $B_\mathrm{unit}, \overline{G}=1$', r'using $B_\phi,\, \overline{G}=1$']
    pmax = get_compared_profiles('max', r, [nmax, 1674, 1676], legend, 'btorgave1')
    psum = get_compared_profiles('sum', r, [nsum, 1675, 1677], legend, 'btorgave1')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [r'using $B_\mathrm{unit}$', r'using $B_\phi$']
    pmax = get_compared_profiles('max', r, [nmax, 1672], legend, 'btor')
    psum = get_compared_profiles('sum', r, [nsum, 1673], legend, 'btor')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = ['Collisional', 'Collisionless']
    pmax = get_compared_profiles('max', r, [nmax, 1666], legend, 'cl')
    psum = get_compared_profiles('sum', r, [nsum, 1667], legend, 'cl')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = ['EM', 'ES']
    pmax = get_compared_profiles('max', r, [nmax, 1668], legend, 'es')
    psum = get_compared_profiles('sum', r, [nsum, 1669], legend, 'es')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [fr'{kxky}$=0.5$', fr'{kxky}$=0.4$', fr'{kxky}$=0.3$']
    pmax = get_compared_profiles('max', r, [nmax, 1678, 1680], legend, 'kxky')
    psum = get_compared_profiles('sum', r, [nsum, 1679, 1681], legend, 'kxky')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [fr'{kyrhos} scan', fr'{kyrhos}$=5$', fr'{kyrhos}$=10$']
    pmax = get_compared_profiles('max', r, [nmax, 1682, 1684], legend, 'kyrhos')
    psum = get_compared_profiles('sum', r, [nsum, 1683, 1685], legend, 'kyrhos')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [r'$n=25$', r'$n=50$', r'$n=100$']
    pmax = get_compared_profiles('max', r, [1816, 1817, 1818], legend, 'kyrhoscount')
    psum = get_compared_profiles('sum', r, [1830, 1831, 1832], legend, 'kyrhoscount')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [r'Converged', r'$n=500$, Linear', r'$n=50,\,\,$ Exponential']
    pmax = get_compared_profiles('max', r, [1822, 1826, 1828], legend, 'kyrhosinc')
    psum = get_compared_profiles('sum', r, [1823, 1827, 1829], legend, 'kyrhosinc')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [fr'{wexb} off', fr'{wexb} on']
    pmax = get_compared_profiles('max', r, [nmax, 1686], legend, 'wexb')
    psum = get_compared_profiles('sum', r, [nsum, 1687], legend, 'wexb')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [fr'{wexb} off', fr'{wexb} on']
    pmax = get_compared_profiles('max', r, [1833, 1835], legend, 'wexbkyrhos', title=fr'Minimum {kyrhos}$=0.1$')
    psum = get_compared_profiles('sum', r, [1834, 1836], legend, 'wexbkyrhos', title=fr'Minimum {kyrhos}$=0.1$')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = ['Default', r'$g_{\mathrm{ne}}=0$']
    pmax = get_compared_profiles('max', r, [nmax, 1688], legend, 'gne0')
    psum = get_compared_profiles('sum', r, [nsum, 1689], legend, 'gne0')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = ['Default', r'$g_{\mathrm{Te}}=0$']
    pmax = get_compared_profiles('max', r, [nmax, 1690], legend, 'gte0')
    psum = get_compared_profiles('sum', r, [nsum, 1691], legend, 'gte0')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = ['Default', r'$g_{\mathrm{ne}}, g_{\mathrm{Te}}=0$']
    pmax = get_compared_profiles('max', r, [nmax, 1692], legend, 'gnegte0')
    psum = get_compared_profiles('sum', r, [nsum, 1693], legend, 'gnegte0')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [r'Saturation$=2$', r'Saturation$=1$']
    pmax = get_compared_profiles('max', r, [nmax, 1694], legend, 'sat')
    psum = get_compared_profiles('sum', r, [nsum, 1695], legend, 'sat')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [fr'{kyrhosmin}$ = 1.00$', fr'{kyrhosmin}$ = 0.50$', fr'{kyrhosmin}$ = 0.25$']
    pmax = get_compared_profiles('max', r, [nmax, 1712, 1714], legend, 'kyrhosmin')
    psum = get_compared_profiles('sum', r, [nsum, 1713, 1715], legend, 'kyrhosmin')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [fr'{ah}$=1$', fr'{ah}$=2$', fr'{ah}$=3$']
    pmax = get_compared_profiles('max', r, [1708, nmax, 1710], legend, 'ah')
    psum = get_compared_profiles('sum', r, [1709, nsum, 1711], legend, 'ah')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [fr'{ah}$=1$', fr'{ah}$=2$', fr'{ah}$=3$']
    pmax = get_compared_profiles('max', r, [1852, 1853, 1854], legend, 'ah2')
    psum = get_compared_profiles('sum', r, [1856, 1855, 1857], legend, 'ah2')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [fr'Default', fr'Geometry Disabled']
    pmax = get_compared_profiles('max', r, [nmax, 1754], legend, 'geo')
    psum = get_compared_profiles('sum', r, [nsum, 1755], legend, 'geo')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [fr'$0\,${bp}', fr'$1\,${bp}', fr'$2\,${bp}']
    pmax = get_compared_profiles('max', r, [1789, nmax, 1791], legend, 'betaep')
    psum = get_compared_profiles('sum', r, [1790, nsum, 1792], legend, 'betaep')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    legend = [r'kpc$=1$', r'kpc$ = c \langle k_\parallel \rangle / (\omega_\mathrm{De} / \overline{G})$', r'kpc$ = c \langle k_\parallel \rangle / \omega_\mathrm{De}$']
    pmax = get_compared_profiles('max', r, [nmax, 1811, 1809], legend, 'kpc')
    psum = get_compared_profiles('sum', r, [nsum, 1812, 1810], legend, 'kpc')
    plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    # Change plot boundary to not include first variable values
    # legend = [r'using $|\hat{\phi}|$', r'using $|\hat{A}_{\!\parallel}\!|$', r'using $|\hat{\phi}|+|\hat{A}_{\!\parallel}\!|$']
    # pmax = get_compared_profiles('max', r, [nmax, 1702, 1704], legend, 'phi')
    # psum = get_compared_profiles('sum', r, [nsum, 1703, 1705], legend, 'phi')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=1)

    """
    Testing isotropic mass change
    """
    # pmax = get_compared_profiles('max', r, [1708, 1847], [fr'{ah}$=1 $ (old)', fr'{ah}$=1 $ (new)'], 'ah')
    # psum = get_compared_profiles('sum', r, [1709, 1849], [fr'{ah}$=1 $ (old)', fr'{ah}$=1 $ (new)'], 'ah')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=0)
    # pmax = get_compared_profiles('max', r, [1710, 1850], [fr'{ah}$=3 $ (old)', fr'{ah}$=3 $ (new)'], 'ah')
    # psum = get_compared_profiles('sum', r, [1711, 1851], [fr'{ah}$=3 $ (old)', fr'{ah}$=3 $ (new)'], 'ah')
    # plot_profiles([pmax, psum], [all_data_max, all_data_sum], saveall=0)
