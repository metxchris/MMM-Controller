"""Stores options needed to run modules in MMM Explorer

Some options are set at runtime by the user, and others may be set while
various code executes.  Options are saved using the pickle module, which
serializes the entire object in binary, and is unreadable by the user.
Consequently, we also save a CSV of all options values stored in the pickle
file, which allows the user to see what options values were stored for each
scan; this CSV is otherwise not used by MMM Explorer.

The Options class is coupled to both the InputControls and Variables classes,
as options values are needed frequently when working with controls or
variables objects.  Specifically, the Options class must be instantiated
before instantiating either InputControls or Variables classes.  When loading
data, options should be loaded first, and then the loaded options object
should be used to load controls and variables objects.

As a consequence of this coupling, the options module does not import either
the controls or variables modules.  Any options related methods requiring
either of these modules are implemented in those modules themselves.

Example Usage:
    # Set options values
    options = modules.options.Options(
        runid='120968A02',
        shot_type=ShotType.NSTX,
        input_time=0.5,
        input_points=51,
    )

    # Save options (to pickle file)
    options.save()

    # Load options (from pickle file)
    options = modules.options.Options()
    options.load('120968A02', 1)
"""

# Standard Packages
import pickle
import inspect
import logging

# 3rd Party Packages
import numpy as np

# Local Packages
import modules.utils as utils
import modules.datahelper as datahelper
import modules.constants as constants
from modules.enums import ShotType, ScanType


_log = logging.getLogger(__name__)

_adjustment_name_to_var_dict = {
    'zeff': 'nz',
    'nuei_alphaconst': 'nuei',
    'nuei_lareunitconst': 'nuei',
    'gne_alphaconst': 'gne',
}


class Options:
    '''
    Stores options for MMM Controller

    Properties and Members:
    * adjustment_name (string): the name of the adjustment to make
    * apply_smoothing (bool): kill-switch to disable smoothing of all variables
    * ignore_exceptions (bool): Exceptions for nonphysical values or NaN values are ignored when True
    * input_points (int): the amount of radial points each variable is interpolated to when sent to MMM
    * input_time (float): the time to check the CDF for values
    * input_time_range (np.ndarray[float]): the input range of time values to use in a time scan
    * normalize_time_range (bool): treat input time range as a range of normalized time values
    * runid (str): the Runid in the CDF, usually also the name of the CDF
    * scan_factor_str (str): the string of the scan factor, rounded for better visual presentation
    * scan_num (int): the number identifying where data is stored within the ./output/runid/ directory
    * scan_range_idxs (list[int]): indices corresponding to scan range values (used for time scans)
    * scan_range (np.ndarray[float]): the range of factors to multiply the var_to_scan by
    * scan_type (ScanType): the type of the scan
    * shot_type (ShotType): the shot type of the CDF
    * temperature_profiles (bool): replace temperature variables with experimental profiles
    * time_str (str): the string of the measurement time, rounded for better visual presentation
    * time_idx (int): the index of the CDF time value that is closest to input_time
    * use_gnezero (bool): set gne equal to zero (sets gne to a small number to avoid division by 0)
    * use_gtezero (bool): set gte equal to zero (sets gte to a small number to avoid division by 0)
    * use_gtizero (bool): set gti equal to zero (sets gte to a small number to avoid division by 0)
    * use_gneabs (bool): take absolute value of gne
    * use_gnethreshold (bool): used when calculating the threshold for gne
    * use_gtethreshold (bool): used when calculating the threshold for gte
    * use_etgm_btor (bool): replaces bu, gbu values with btor, gbtor values in the input file
    * use_experimental_profiles (bool): replaces te, ti, q, with CDF variables TEPRO, TIPRO, QPRO (if available)
    * var_to_scan (str): the variable to scan; syntax must match a member of InputVariables
    '''

    def __init__(self, **kwargs):
        # Private members (each has a property)
        self._adjustment_name = None
        self._input_points = None
        self._runid = None
        self._scan_range = None
        self._time_str = None
        self._var_to_scan = None
        # Public members
        self.apply_smoothing = False
        self.ignore_exceptions = False
        self.input_time = None
        self.input_time_range = None
        self.normalize_time_range = True
        self.scan_num = None
        self.scan_range_idxs = None
        self.scan_type = ScanType.NONE
        self.shot_type = ShotType.NONE
        self.temperature_profiles = False
        # self.time_str = None
        self.time_idx = None
        self.use_gnezero = False
        self.use_gtezero = False
        self.use_gtizero = False
        self.use_gneabs = False
        self.use_gnethreshold = False
        self.use_gtethreshold = False
        self.use_etgm_btor = False
        self.use_experimental_profiles = False

        self.set(**kwargs)

    # Properties
    @property
    def adjustment_name(self):
        return self._adjustment_name

    @adjustment_name.setter
    def adjustment_name(self, adjustment_name):
        self._adjustment_name = adjustment_name

        if adjustment_name in _adjustment_name_to_var_dict:
            self.var_to_scan = _adjustment_name_to_var_dict[adjustment_name]
        else:
            self.var_to_scan = adjustment_name

    @property
    def input_points(self):
        return self._input_points

    @input_points.setter
    def input_points(self, points):
        if points:
            self._input_points = max(points, 5)

    @property
    def runid(self):
        return self._runid

    @runid.setter
    def runid(self, runid):
        if not isinstance(runid, str):
            runid = str(runid)
        self._runid = runid.strip()

    @property
    def scan_range(self):
        return self._scan_range

    @scan_range.setter
    def scan_range(self, scan_range):
        if scan_range is not None:
            if not isinstance(scan_range, np.ndarray):
                raise TypeError(f'scan_range must be {np.ndarray} or {None} and not {type(scan_range)}')
            too_small = np.absolute(scan_range) < constants.ABSMIN_SCAN_FACTOR_VALUE
            if too_small.any():
                value_signs = np.sign(scan_range[too_small])
                value_signs[value_signs == 0] = 1  # np.sign(0) = 0, so set these to +1
                scan_range[too_small] = constants.ABSMIN_SCAN_FACTOR_VALUE * value_signs
        self._scan_range = scan_range

    @property
    def time_str(self):
        return self._time_str

    @time_str.setter
    def time_str(self, time_value):
        if isinstance(time_value, str):
            self._time_str = time_value
        else:
            self._time_str = f'{time_value:{constants.TIME_VALUE_FMT}}'

    @property
    def var_to_scan(self):
        return self._var_to_scan

    @var_to_scan.setter
    def var_to_scan(self, var_to_scan):
        self._var_to_scan = var_to_scan
        self.scan_type = datahelper.get_scan_type(var_to_scan)

    # Methods
    def get_keys(self):
        '''Returns (list): Names of properties and public members'''
        kvps = self.get_key_value_pairs()
        return [kvp[0] for kvp in kvps]

    def get_key_value_pairs(self):
        '''Returns (list): Names and values of properties and public members'''
        attributes = inspect.getmembers(self, lambda m: not inspect.isroutine(m))
        return [a for a in attributes if not a[0].startswith('_')]

    def set(self, **kwargs):
        '''Sets specified options values'''
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                _log.warning(f'\n\tOptions class does not contain member {key}')

    def load(self, runid, scan_num):
        '''
        Loads Options object from a pickle file

        Parameters:
        * runid (str): The run id to load options from
        * scan_num (int): The scan number to load options from

        Returns:
        * self (Options)
        '''

        pickle_path = utils.get_options_path(runid, scan_num)

        with open(pickle_path, 'rb') as handle:
            loaded_options = pickle.load(handle)

        # Setting options values one-by-one will not break any existing references to Options
        options_to_set = loaded_options.get_keys()
        for option in options_to_set:
            setattr(self, option, getattr(loaded_options, option))

        return self

    def save(self):
        '''
        Saves Options object to a pickle file

        A CSV is also saved to make saved options viewable by the user, but is otherwise not used
        '''

        pickle_path = utils.get_options_path(self.runid, self.scan_num)

        with open(pickle_path, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # Options are also written to a CSV to make viewing their values easier
        csv_path = pickle_path.replace('.pickle', '.csv')
        options_values = self.get_key_value_pairs()
        with open(csv_path, 'w') as handle:
            handle.write('# Note: This file is for inspection of values, and is not otherwise used by MMM Explorer\n')
            for o in options_values:
                handle.write(f'{o[0]}, {o[1]}\n')

        _log.info(f'\n\tSaved: {pickle_path}\n')

    def set_measurement_time(self, time_values):
        '''Find the index of the measurement time closest to the input_time, then store that value and its index'''
        if self.input_time is not None:
            self.time_idx = np.argmin(np.abs(time_values - self.input_time))
            self.time_str = time_values[self.time_idx]

    def set_measurement_rho(self, rho_values):
        '''Find the index of the measurement time closest to the input_time, then store that value and its index'''
        if self.input_time is not None:
            self.time_idx = np.argmin(np.abs(time_values - self.input_time))
            self.time_str = time_values[self.time_idx]

    def set_time_ranges(self, time_values):
        '''
        Set the time range and time range indices using time values from the CDF

        This is only needed when conducting a time scan

        Parameters:
        * time (np.ndarray(float)): Time values from the CDF
        '''

        if self.scan_range is not None:

            times = time_values
            if self.normalize_time_range:
                # Compare normalized input time range against normalized time values
                times = (time_values - time_values.min()) / (time_values.max() - time_values.min())

            rounded_range = np.round(self.scan_range, constants.TIME_VALUE_SIGFIGS)
            times_tile = np.tile(np.round(times, constants.TIME_VALUE_SIGFIGS), (rounded_range.shape[0], 1))
            rounded_range_tile = np.tile(rounded_range, (1, 1)).T
            unique_idxs = np.unique(np.argmin(np.abs(times_tile - rounded_range_tile), axis=1))

            scan_range_idxs = []
            scan_range_time_strs = []
            for i in unique_idxs:
                time_idx = np.argmin(np.abs(time_values - time_values[i]))
                time_str = f'{time_values[time_idx]:{constants.TIME_VALUE_FMT}}'
                if time_str not in scan_range_time_strs:
                    scan_range_time_strs.append(time_str)
                    scan_range_idxs.append(i)

            self._scan_range = np.array(scan_range_time_strs).astype(float)
            self.scan_range_idxs = scan_range_idxs

    def find_scan_factor(self, scan_factor):
        '''
        Finds the value in scan range closest to the supplied scan factor

        Parameters:
        * scan_factor (float): The scan factor to find

        Returns:
        * (str | None): Value in scan_range closest to the specified
          scan_factor, and None if scan_factor is None

        Raises:
        * ValueError: If scan_range is None
        '''

        if scan_factor is None:
            return scan_factor
        elif self.scan_range is None:
            raise ValueError('Cannot find scan_factor value when scan_range is None')
        return self.scan_range[np.argmin(np.abs(self.scan_range - scan_factor))]
