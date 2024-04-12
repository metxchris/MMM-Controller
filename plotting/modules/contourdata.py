#!/usr/bin/python3

# Standard Packages
import logging

# 3rd Party Packages
import scipy.ndimage
import numpy as np

# Local Packages
import modules.utils as utils
import modules.datahelper as datahelper
import modules.constants as constants
import modules.options
from modules.variables import InputVariables, OutputVariables


_log = logging.getLogger(__name__)


class PlotOptions:
    """ Store options for the contour plot

    Initialization parameters determine plot settings and must be specified as
    keyword arguments.
    * savenameend (str): Name to end file save name with
    """

    def __init__(self, **kwargs):
        self.ymin: float | None = None
        self.ymax: float | None = None
        self.xmin: float | None = None
        self.xmax: float | None = None
        self.zmin_diff: float | None = None
        self.zmax_diff: float | None = None
        self.time_min: float = 0
        self.time_max: float = 1
        self.time_count: int = 100
        self.savefig: bool = False
        self.savedata: bool = False
        self.smoothing: bool = False
        self.showfig: bool = True
        self.raw: bool = False
        self.plotidentical: bool = False
        self.difftype: str = ''
        self.saveappend: str = ''

        self._set_kwargs(kwargs)
        self._validate_members()

    def _set_kwargs(self, kwargs):
        """Set member values using keyword arguments"""
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                _log.error(f'\n\t"{key}" is not a valid parameter for PlotOptions')

    def _validate_members(self):
        if self.difftype:
            difftypes = ['diff', 'absdiff', 'absrel', 'rel', 'ratio']
            if self.difftype not in difftypes:
                raise ValueError(
                    f'difftype \'{self.difftype}\' not in list of accepted types:'
                    f'\n\t{difftypes}'
                )


class ContourData():

    def __init__(self, options, plot_options):
        self.options = options
        self.plot_options = plot_options
        self.var_to_plot = None
        self.xvar = None
        self.yvar = None
        self.zvar = None
        self.X: np.ndarray | None = None
        self.Y: np.ndarray | None = None
        self.Z: np.ndarray | None = None

    def _get_variable(self, name, obj_list):
        """Returns variable object from Variables class"""
        var = None
        for obj in obj_list:
            if hasattr(obj, name):
                var = getattr(obj, name)

        return var

    def init_meshgrid(self, x, y):

        if self.plot_options is not None:
            idx = np.ones_like(x, dtype=bool)
            if self.plot_options.xmin is not None:
                idx *= x >= self.plot_options.xmin
            if self.plot_options.xmax is not None:
                idx *= x <= self.plot_options.xmax
            x = x[idx]

            idx = np.ones_like(y, dtype=bool)
            if self.plot_options.ymin is not None:
                idx *= y >= self.plot_options.ymin
            if self.plot_options.ymax is not None:
                idx *= y <= self.plot_options.ymax
            y = y[idx]

        self.X, self.Y = np.meshgrid(x, y)

    def verify_z(self):
        if np.isnan(self.Z).any():
            # This can happen due to calculation errors or or when loading
            # data using old variable definitions, so a warning is thrown but
            # no exception is raised
            _log.warning(f'\n\tnan values found in {self.var_to_plot}, continuing to next variable...')
            return True  # Error flag

        if self.use_log():
            if self.Z.max() <= 0:
                _log.warning(f'\n\tcan\'t use log for {self.var_to_plot} when no positive values exist...')
                return True  # Error flag
            self.Z = np.ma.masked_where(self.Z <= 0, self.Z)

        if self.plot_options.smoothing:
            self.smooth_data(self.var_to_plot)

    def smooth_data(self):
        """Smooth contour data using a Gaussian filter"""
        Zmax, Zmin = self.Z.max(), self.Z.min()  # Cache bounds before smoothing

        # TODO: attach all of this to individual variables
        medium_smoothing = ['xteETGM', 'xte2ETGM', 'xdiETGM']
        light_smoothing = [
            'omgETGM', 'kyrhoeETGM', 'kyrhosETGM',
            'omegateETGM', 'walfvenunit', 'omegadETGM',
            'omegasETGM', 'omegasetaETGM', 'omegadiffETGM',
            'gammadiffETGM'
        ]

        sigma = (0, 0)
        if self.options.var_to_scan == 'time':
            sigma = (0, 0)
        elif self.var_to_plot in medium_smoothing:
            sigma = (2, 2)
        elif self.var_to_plot in light_smoothing:
            sigma = (1, 1)

        self.Z = scipy.ndimage.gaussian_filter(self.Z, sigma=sigma)
        self.Z = np.minimum(np.maximum(self.Z, Zmin), Zmax)  # restore original bounds

    def get_xlabel(self):
        """Get the xlabel of the plot"""
        return self.xvar.label

    def get_ylabel(self):
        """Get the ylabel of the plot"""
        ylabel = self.yvar.label
        var_to_scan = self.options.var_to_scan
        if var_to_scan == 'gne' and self.options.use_gneabs:
            ylabel = f'$|${ylabel}$|$'

        if var_to_scan and 'time' not in var_to_scan and 'kyrho' not in var_to_scan:  # Most yvariables will be plotted as multipliers
            ylabel = f'{ylabel} (multipliers)'
        elif self.yvar.units_label:
            ylabel = f'{ylabel} ({self.yvar.units_label})'
        return ylabel

    def get_title(self):
        """Get the plot title"""
        title = f'{self.zvar.label} ({self.zvar.units_label})' if self.zvar.units_label else f'{self.zvar.label}'
        if self.options.use_gnezero:
            title = fr'{title} [$g_\mathrm{{ne}} = 0$]'
        if self.options.use_gtezero:
            title = fr'{title} [$g_\mathrm{{Te}} = 0$]'
        if self.options.use_gneabs:
            title = fr'{title} with $|g_\mathrm{{ne}}|$'
        if hasattr(self, 'controls') and self.controls is not None:
            # Could use some cleanup
            if self.controls.etgm_sum_modes.values and 'ETGM' in self.var_to_plot and '\chi' in self.zvar.label:
                title = fr'$_{{^\sum}}${title}'
        if self.var_to_plot in ['nR8TOMSQZ', 'nCubic', 'nSolver']:
            title = f'{title} ({round(np.average(self.Z[:, 1:]), 1)})'
        if self.var_to_plot in ['nWarning', 'nError']:
            title = f'{title} ({int(np.sum(self.Z))})'
        return title

    def get_savename(self):
        """Get the name to save the file as"""
        savename = self.options.adjustment_name or self.options.var_to_scan or 'cdf'
        if self.options.use_gnezero:
            savename = f'{savename}_gne0'
        if self.options.use_gtezero:
            savename = f'{savename}_gte0'
        if self.options.use_gneabs:
            savename = f'{savename}_gneabs'
        if hasattr(self, 'controls') and self.controls is not None:
            if self.controls.etgm_exbs.values:
                savename = f'{savename}_exbs1'
            if self.controls.etgm_sum_modes.values and 'xte' in self.var_to_plot:
                savename = f'{savename}_sum'
        if self.plot_options.difftype:
            savename = f'{savename}_{self.plot_options.difftype}'
        if self.options.scan_num:
            savename = f'{savename}_{self.options.scan_num}'
        if hasattr(self, 'options2') and self.options2.scan_num:
            savename = f'{savename}_{self.options2.scan_num}'
        if self.plot_options.saveappend:
            savename = f'{savename}_{self.plot_options.saveappend}'
        return f'{savename}_{self.var_to_plot}'

    def save_to_csv(self, savename):
        """
        Save plotted data to a CSV

        Data is saved into three CSV's, one for each dimension X, Y, and Z.  When
        using normalized data for X and Y, the CSV's could be simplified and
        saved as single arrays.  However, the full matrices are being saved
        anyways if generalized values for X and Y are ever used.

        Raises:
        * FileNotFoundError: If the file cannot be found after saving it
        """

        prec, col_pad = 4, 6
        col_len = prec + col_pad
        fmt_str = f'%{col_len}.{prec}e'

        # Deep copy np arrays
        X = np.array(self.X)
        Y = np.array(self.Y)
        Z = np.array(self.Z)

        # Crop X, Y data if values are uniform per row or column, respectively
        if (X.max(0) == X.min(0)).all():
            X = X[0, :]
        if (Y.max(1) == Y.min(1)).all():
            Y = Y[:, 0]

        data_to_save = [X, Y, Z]
        names_to_save = ['X', 'Y', 'Z']

        print(f'\nPlot Data Saved:')
        for data, name in zip(data_to_save, names_to_save):
            var_savename = f'{savename}_{name}.csv'
            np.savetxt(var_savename, data, fmt=fmt_str, delimiter=',')
            if not utils.check_exists(var_savename):
                raise FileNotFoundError(f'Failed to save {var_savename}\n'
                                        '\tMake sure Python has file writing permissions')
            print(f'\t{var_savename}')

    def use_log(self):
        """Check whether the axis should be logarithmic"""
        return self.var_to_plot in ['errorEPM']


class ContourDataMMM(ContourData):

    def __init__(self, options, plot_options):
        super().__init__(options, plot_options)  # Init parent class

        self.input_vars_rho = None
        self.output_vars_rho = None
        self.controls_rho = None
        self.input_vars = None
        self.output_vars = None
        self.controls = None
        self.rho_strs: list | None = None

        i, o, c = datahelper.get_all_rho_data(self.options)
        self.input_vars_rho = i
        self.output_vars_rho = o
        self.controls_rho = c

        i, o, c = datahelper.get_data_objects(self.options)
        self.input_vars = i
        self.output_vars = o
        self.controls = c

        self.xvar = self.input_vars.rho
        self.yvar = self._get_variable(self.options.var_to_scan, [self.input_vars, self.output_vars, self.controls])

        x = np.array(list(self.output_vars_rho.keys()), dtype=float)
        y = self.options.scan_range
        self.init_meshgrid(x, y)

        self.rho_strs = [f'{val:{constants.RHO_VALUE_FMT}}' for val in self.X[0, :]]

    def set_z(self, var_to_plot):
        self.var_to_plot = var_to_plot

        if hasattr(self.output_vars, var_to_plot):
            self.zvar = getattr(self.output_vars, var_to_plot)
        elif hasattr(self.input_vars, var_to_plot):
            self.zvar = getattr(self.input_vars, var_to_plot)
        elif hasattr(self.controls, var_to_plot):
            self.zvar = getattr(self.controls, var_to_plot)
        if not self.zvar or not hasattr(self.zvar, 'values') or not isinstance(self.zvar.values, np.ndarray):
            # Contours need to be np.ndarray. This can happen when using the
            # same variable plotting list for many different scan types, so
            # no exception is raised
            _log.warning(f'\n\t{var_to_plot} is not an np.ndarray, continuing to next variable...')
            return True  # Error flag

        Z_dict = None
        if hasattr(self.output_vars, var_to_plot):
            Z_dict = self.output_vars_rho
        elif hasattr(self.input_vars, var_to_plot):
            Z_dict = self.input_vars_rho

        if Z_dict is None:
            _log.warning(f'\n\tFailed to load dictionary csv for {var_to_plot}, continuing to next variable...')
            return True  # Error flag

        y = self.options.scan_range
        idxs = np.ones_like(y, dtype=bool)
        if self.plot_options.ymin is not None:
            idxs *= y >= self.plot_options.ymin
        if self.plot_options.ymax is not None:
            idxs *= y <= self.plot_options.ymax

        self.Z = np.zeros_like(self.X)
        for i, rho_str in enumerate(self.rho_strs):
            self.Z[:, i] = getattr(Z_dict[rho_str], var_to_plot).values[idxs]

        self.verify_z()

    def verify_z(self):
        super().verify_z()

        # Clamp Z between min and max defined on variable
        self.Z = np.minimum(np.maximum(self.Z, self.zvar.contour_min), self.zvar.contour_max)


class ContourDataCDF(ContourData):

    def __init__(self, options, plot_options):
        super().__init__(options, plot_options)  # Init parent class

        mmm_vars, __, raw_vars = datahelper.initialize_variables(options)

        self.mmm_vars = mmm_vars if not plot_options.raw else raw_vars
        self.xvar = mmm_vars.rho if not plot_options.raw else raw_vars.xb
        self.yvar = self.mmm_vars.time

        x = self.xvar.values[:, 0]
        y = self.yvar.values

        self.init_meshgrid(x, y)
        self.X0, self.Y0 = np.meshgrid(x, y)

    def set_z(self, var_to_plot):
        self.var_to_plot = var_to_plot

        self.zvar = getattr(self.mmm_vars, var_to_plot)
        if not self.zvar or not hasattr(self.zvar, 'values') or not isinstance(self.zvar.values, np.ndarray):
            # Contours need to be np.ndarray. This can happen when using the
            # same variable plotting list for many different scan types, so
            # no exception is raised
            _log.warning(f'\n\t{var_to_plot} is not an np.ndarray (check if raw)...')
            return True  # Error flag

        idx = np.ones_like(self.X0, dtype=bool)
        if self.plot_options.xmin is not None:
            idx *= self.X0 >= self.plot_options.xmin
        if self.plot_options.xmax is not None:
            idx *= self.X0 <= self.plot_options.xmax
        if self.plot_options.ymin is not None:
            idx *= self.Y0 >= self.plot_options.ymin
        if self.plot_options.ymax is not None:
            idx *= self.Y0 <= self.plot_options.ymax

        self.Z = self.zvar.values.T
        self.Z = self.Z[idx].reshape(self.X.shape)

        self.verify_z()

    def verify_z(self):
        super().verify_z()

        # Clamp Z between min and max defined on variable
        self.Z = np.minimum(np.maximum(self.Z, self.zvar.contour_min), self.zvar.contour_max)


class ContourDataDiff(ContourData):

    def __init__(self, options, options2, plot_options):
        super().__init__(options, plot_options)  # Init parent class

        self.options2 = options2

        self.isidentical: bool = False  # If both vars are identical
        self.input_vars_rho = None
        self.output_vars_rho = None
        self.controls_rho = None
        self.input_vars_rho2 = None
        self.output_vars_rho2 = None
        self.controls_rho2 = None
        self.input_vars = None
        self.output_vars = None
        self.controls = None
        self.input_vars2 = None
        self.output_vars2 = None
        self.controls2 = None
        self.zvar = None
        self.zvar2 = None
        self.rho_strs: list | None = None

        i, o, c = datahelper.get_all_rho_data(self.options)
        self.input_vars_rho = i
        self.output_vars_rho = o
        self.controls_rho = c

        i, o, c = datahelper.get_all_rho_data(self.options2)
        self.input_vars_rho2 = i
        self.output_vars_rho2 = o
        self.controls_rho2 = c

        i, o, c = datahelper.get_data_objects(self.options)
        self.input_vars = i
        self.output_vars = o
        self.controls = c

        i, o, c = datahelper.get_data_objects(self.options2)
        self.input_vars2 = i
        self.output_vars2 = o
        self.controls2 = c

        self.rho_strs = self.input_vars_rho.keys()
        if (self.rho_strs != self.input_vars_rho2.keys()):
            raise ValueError('Rho values must match from each data set')

        self.xvar = self.input_vars.rho
        self.yvar = self._get_variable(self.options.var_to_scan, [self.input_vars, self.output_vars, self.controls])

        x = np.array(list(self.rho_strs), dtype=float)
        y = self.options.scan_range
        self.init_meshgrid(x, y)

        self.rho_strs = [f'{val:{constants.RHO_VALUE_FMT}}' for val in self.X[0, :]]

    def set_z(self, var_to_plot):
        self.var_to_plot = var_to_plot

        self.zvar = self._get_variable(var_to_plot, [self.input_vars, self.output_vars, self.controls])
        if not self.zvar or not hasattr(self.zvar, 'values') or not isinstance(self.zvar.values, np.ndarray):
            # Contours need to be np.ndarray. This can happen when using the
            # same variable plotting list for many different scan types, so
            # no exception is raised
            _log.warning(f'\n\t{var_to_plot} is not an np.ndarray, continuing to next variable...')
            return True  # Error flag

        self.zvar2 = self._get_variable(var_to_plot, [self.input_vars2, self.output_vars2, self.controls2])
        if not self.zvar2 or not hasattr(self.zvar2, 'values') or not isinstance(self.zvar2.values, np.ndarray):
            # Contours need to be np.ndarray. This can happen when using the
            # same variable plotting list for many different scan types, so
            # no exception is raised
            _log.warning(f'\n\t{var_to_plot} is not an np.ndarray, continuing to next variable...')
            return True  # Error flag

        Z_dict = None
        if hasattr(self.output_vars, var_to_plot):
            Z_dict = self.output_vars_rho
        elif hasattr(self.input_vars, var_to_plot):
            Z_dict = self.input_vars_rho
        if Z_dict is None:
            _log.warning(f'\n\tFailed to load dictionary csv for {var_to_plot}, continuing to next variable...')
            return True  # Error flag

        Z_dict2 = None
        if hasattr(self.output_vars2, var_to_plot):
            Z_dict2 = self.output_vars_rho2
        elif hasattr(self.input_vars2, var_to_plot):
            Z_dict2 = self.input_vars_rho2
        if Z_dict2 is None:
            _log.warning(f'\n\tFailed to load dictionary csv for {var_to_plot}, continuing to next variable...')
            return True  # Error flag

        Z1 = np.zeros_like(self.X)
        Z2 = np.zeros_like(self.X)

        y = self.options.scan_range
        idxs = np.ones_like(y, dtype=bool)
        if self.plot_options.ymin is not None:
            idxs *= y >= self.plot_options.ymin
        if self.plot_options.ymax is not None:
            idxs *= y <= self.plot_options.ymax
        for i, rho_str in enumerate(self.rho_strs):
            Z1[:, i] = getattr(Z_dict[rho_str], self.var_to_plot).values[idxs]
            Z2[:, i] = getattr(Z_dict2[rho_str], self.var_to_plot).values[idxs]

        if (Z1 == Z2).all():
            _log.warning(f'\n\tIdentical data sets found for {var_to_plot}')
            if not self.plot_options.plotidentical:
                return True  # Error flag
            else:
                self.isidentical = True

        # Take difference
        Zdiff = Z1 - Z2
        Zsum = np.absolute(Z1 + Z2)
        Zsum[Zsum == 0] = 1  # Divide by zero errors

        if self.plot_options.difftype == 'diff':
            self.Z = Zdiff
        elif self.plot_options.difftype == 'absdiff':
            self.Z = np.absolute(Zdiff)
        elif self.plot_options.difftype == 'rel':
            self.Z = Zdiff / Zsum
        elif self.plot_options.difftype == 'absrel':
            self.Z = np.absolute(Zdiff / Zsum)
        elif self.plot_options.difftype == 'ratio':
            Z1[Z1 == 0] = 1
            Z2[Z2 == 0] = 1
            self.Z = Z1 / Z2
        else:
            _log.error(f'\n\tNo formula defined for difference type: {self.plot_options.difftype}')
            return True  # Error flag

        if self.plot_options.zmax_diff:
            self.Z = np.minimum(self.Z, self.zvar.contour_max)
        if self.plot_options.zmin_diff:
            self.Z = np.maximum(self.Z, self.zvar.contour_min)

        # self.Z -= 1
        # self.Z[self.Z > 2] = 2
        # self.Z[self.Z < -2] = -2

        self.verify_z()

    def get_title(self):
        title = f'{self.zvar.label} ({self.zvar.units_label})' if self.zvar.units_label else f'{self.zvar.label}'
        print(self.var_to_plot)
        if self.var_to_plot in ['nR8TOMSQZ', 'nCubic', 'nSolver']:
            title = f'{title} ({round(np.average(self.Z[:, 1:]), 1)})'
        elif self.var_to_plot in ['nWarning', 'nError']:
            title = f'{title} ({int(np.sum(self.Z))})'

        if self.plot_options.difftype == 'diff':
            title = fr'{title} [Diff]'
        elif self.plot_options.difftype == 'absdiff':
            title = fr'{title} [Abs Diff]'
        elif self.plot_options.difftype == 'rel':
            title = fr'{title} [Rel]'
        elif self.plot_options.difftype == 'absrel':
            title = fr'{title} [Abs Rel]'
        elif self.plot_options.difftype == 'ratio':
            title = fr'{title} [Ratio]'

        return title


def verify_vars_to_plot(vars_to_plot):
    '''
    Variables are verified to be members of OutputVariables before the plotting loop begins

    Parameters:
    * vars_to_plot (list):  List of output variables to plot

    Raises:
    * NameError: If the variable to plot is not found in OutputVariables
    '''

    options = modules.options.Options()
    output_vars = OutputVariables(options)
    input_vars = InputVariables(options)
    for var_to_plot in vars_to_plot:
        if var_to_plot == 'var_to_scan':
            continue
        if not hasattr(output_vars, var_to_plot) and not hasattr(input_vars, var_to_plot):
            raise NameError(f'{var_to_plot} not found in Variables classes')
